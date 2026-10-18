"""src/run_vasttools_turbo.py against a mock vast-tools and a mock bowtie.

The real tools are compared by RealTurboEquivalenceTests, which runs only when
VASTTOOLS_TEST_DIR (a vast-tools checkout), VASTTOOLS_TEST_DB (a VASTDB folder
holding e.g. Hsa/), VASTTOOLS_TEST_SPECIES (e.g. hg19) and VASTTOOLS_TEST_FASTQ
(a gzipped FASTQ with a few hundred thousand reads) are set and bowtie is on PATH.
"""

import gzip
import json
import os
import shutil
import stat
import sys
import tarfile
import tempfile
import textwrap
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))

from src import run_vasttools_turbo as turbo  # noqa: E402

# mock bowtie: one output line per input read; the Nth call (one per mate) calls
# round(fraction * reads) reads "-" strand, fractions from MOCK_FRACS ("0.5,0.5")
BOWTIE = textwrap.dedent('''\
    #!/usr/bin/env python3
    import os, sys
    from pathlib import Path
    counter = Path(os.environ["MOCK_BT_COUNT"])
    call = int(counter.read_text()) if counter.exists() else 0
    counter.write_text(str(call + 1))
    fractions = [float(x) for x in os.environ["MOCK_FRACS"].split(",")]
    fraction = fractions[min(call, len(fractions) - 1)]
    mapped = int(os.environ.get("MOCK_BT_MAPPED", "100000"))
    lines = sys.stdin.read().splitlines()
    reads = len(lines) // 4
    minus = round(fraction * min(reads, mapped))
    for i in range(min(reads, mapped)):
        print("\\t".join(["r%d" % i, "-" if i < minus else "+", "mRNA", "0", "ACGT"]))
''')

# mock vast-tools: align writes the to_combine files (info says -SS when no --ns);
# merge sums nothing but writes every group file and the group info like MergeOutputs.pl
VAST = textwrap.dedent('''\
    #!/usr/bin/env python3
    import json, os, sys
    from pathlib import Path
    args = sys.argv[1:]
    command = args[0]
    def value(flag):
        return args[args.index(flag) + 1]
    with open(os.environ["MOCK_CALLS"], "a") as stream:
        stream.write(json.dumps(args) + "\\n")
    skip = os.environ.get("MOCK_SKIP", "")
    suffixes = (".eej2", ".exskX", ".micX", ".MULTI3X", ".IR2", ".IR.summary_v2.txt")
    if command == "align":
        out = Path(value("-o")) / "to_combine"
        out.mkdir(parents=True)
        name = value("--name")
        for suffix in suffixes:
            if suffix != skip:
                (out / (name + suffix)).write_text("" if os.environ.get("MOCK_EMPTY") == suffix else "x\\t1\\n")
        stranded = "--ns" not in args
        tail = "\\t--norc\\t-SS\\tdone" if stranded else "\\t\\tdone"
        (out / (name + ".info")).write_text("single\\tR1\\tHsa\\t--s" + tail)
    elif command == "merge":
        groups = [line.split("\\t") for line in Path(value("--groups")).read_text().splitlines()]
        out = Path(value("-o")) / "to_combine"
        for sub, sample in groups:
            assert (out / (sub + ".eej2")).exists(), "merge needs plain files"
            assert "--sp" in args, "merge needs --sp"
        infos = [(out / (sub + ".info")).read_text() for sub, _ in groups]
        ss = all(i.split("\\t")[-2] == "-SS" for i in infos)
        for sample in {g[1] for g in groups}:
            for suffix in suffixes:
                (out / (sample + suffix)).write_text("merged\\t%d\\n" % len(groups))
            (out / (sample + ".info")).write_text("--norc\\t-SS\\tdone" if ss else "\\t\\tdone")
        if "--move_to_PARTS" in args:
            (out / "PARTS").mkdir(exist_ok=True)
            for sub, _ in groups:
                for path in out.glob(sub + ".*"):
                    path.rename(out / "PARTS" / path.name)
''')


def fastq(path, reads, mate=1):
    with gzip.open(path, "wt") as stream:
        for i in range(reads):
            stream.write("@r{}/{}\nACGTACGTAC\n+\nIIIIIIIIII\n".format(i, mate))


def names(path):
    with gzip.open(path, "rt") as stream:
        return [line.strip()[1:] for n, line in enumerate(stream) if n % 4 == 0]


class Base(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp(prefix="turbo-test-"))
        self.addCleanup(shutil.rmtree, self.dir, True)
        self.bin = self.dir / "bin"
        self.bin.mkdir()
        for name, text in (("bowtie", BOWTIE), ("vast-tools", VAST)):
            path = self.bin / name
            path.write_text(text)
            path.chmod(path.stat().st_mode | stat.S_IEXEC)
        self.calls = self.dir / "calls.jsonl"
        self.old = dict(os.environ)
        os.environ.update(PATH=str(self.bin) + os.pathsep + os.environ["PATH"], MOCK_CALLS=str(self.calls),
                          MOCK_BT_COUNT=str(self.dir / "bt_count"), MOCK_FRACS="0.5")
        for key in ("MOCK_SKIP", "MOCK_EMPTY", "MOCK_BT_MAPPED"):
            os.environ.pop(key, None)
        self.addCleanup(self.restore)
        self.db = self.dir / "db"
        (self.db / "Hsa" / "EXPRESSION").mkdir(parents=True)

    def restore(self):
        os.environ.clear()
        os.environ.update(self.old)

    def vast_calls(self, command):
        if not self.calls.exists():
            return []
        return [c for c in map(json.loads, self.calls.read_text().splitlines()) if c[0] == command]

    def reads(self, count=3000, mates=1):
        paths = []
        for mate in range(1, mates + 1):
            path = self.dir / "R{}.fastq.gz".format(mate)
            fastq(path, count, mate)
            paths.append(path)
        return paths

    def split(self, chunks=3, count=3000, mates=1, **kwargs):
        reads = self.reads(count, mates)
        out = self.dir / "work"
        record = turbo.split([str(p) for p in reads], self.db, "hg19", chunks, out, "run1", **kwargs)
        return out, record

    def align_chunks(self, out, run="run1"):
        state = json.loads((out / "split.json").read_text())
        tars, metas = [], []
        for k in range(state["chunks"]):
            files = sorted(str(p) for p in (out / "chunks").glob("c{}_R*.fastq.gz".format(k)))
            tar, meta = self.dir / "chunk_{}.tar".format(k), self.dir / "chunk_{}.json".format(k)
            turbo.chunk(self.bin, self.db, "hg19", out / "split.json", k, files, run, tar, meta, threads=2)
            tars.append(str(tar))
            metas.append(str(meta))
        return tars, metas


class IdentityTests(unittest.TestCase):
    """The turbo module is independent of the native one and keyed by what changes its counts."""

    def ids(self, **config):
        from src.umbrella_modules import module_reference_id
        return {tool: module_reference_id("ref", tool, config) for tool in ("vasttools", "vasttools_turbo")}

    def test_turbo_and_native_have_separate_identities(self):
        ids = self.ids()
        self.assertNotEqual(ids["vasttools"], ids["vasttools_turbo"])

    def test_native_identity_is_unchanged_by_turbo_settings(self):
        self.assertEqual(self.ids()["vasttools"], self.ids(vasttools_turbo_chunks=3, vasttools_turbo_threads=9)["vasttools"])

    def test_chunks_change_the_turbo_identity_but_threads_do_not(self):
        base = self.ids()["vasttools_turbo"]
        self.assertEqual(base, self.ids(vasttools_turbo_threads=2)["vasttools_turbo"])
        self.assertNotEqual(base, self.ids(vasttools_turbo_chunks=4)["vasttools_turbo"])
        self.assertEqual(base, self.ids(vasttools_turbo_chunks=8)["vasttools_turbo"])

    def test_species_and_expr_still_count(self):
        base = self.ids()["vasttools_turbo"]
        self.assertNotEqual(base, self.ids(vasttools_species="mm10")["vasttools_turbo"])
        self.assertNotEqual(base, self.ids(vasttools_expr=True)["vasttools_turbo"])

    def test_the_module_is_selectable_and_opt_in(self):
        from src.umbrella_modules import parse_selection
        self.assertEqual(parse_selection({"umbrella_modules": ["vasttools_turbo", "vasttools"]})[0],
                         ["vasttools", "vasttools_turbo"])
        from src import umbrella_execution
        self.assertIn("vasttools_turbo", umbrella_execution.COMPLEMENT)
        self.assertNotIn("vasttools_turbo", umbrella_execution.DEFAULT_COMPLEMENT)
        self.assertEqual(umbrella_execution.tool_for("umbrella_vasttools_turbo_chunk", ""), "vasttools_turbo")
        self.assertEqual(umbrella_execution.tool_for("umbrella_vasttools_align", ""), "vasttools")


class StrandTests(unittest.TestCase):
    def test_decision_follows_the_native_thresholds(self):
        cases = [([0.5], False, []), ([0.35], False, []), ([0.65], False, []), ([0.3499], True, [False]),
                 ([0.6501], True, [True]), ([0.05], True, [False]), ([0.99], True, [True]),
                 ([0.5, 0.5], False, []), ([0.5, 0.66], True, [False, True]),
                 ([0.9, 0.1], True, [True, False]), ([0.4, 0.2], True, [False, False])]
        for fractions, stranded, rc in cases:
            with self.subTest(fractions=fractions):
                decision = turbo.strand_decision(fractions)
                self.assertEqual(decision["strand_specific"], stranded)
                self.assertEqual(decision["rc"], rc)

    def test_flags_force_every_mate(self):
        self.assertEqual(turbo.strand_flags({"strand_specific": False, "rc": []}), ["--ns"])
        self.assertEqual(turbo.strand_flags({"strand_specific": True, "rc": [True]}), ["--rc1"])
        self.assertEqual(turbo.strand_flags({"strand_specific": True, "rc": [False, True]}),
                         ["--nrc1", "--rc2"])

    def test_info_strandedness_is_read_like_merge_does(self):
        self.assertTrue(turbo._info_is_stranded("single\tR1\tHsa\t--s\t--norc\t-SS\tdone"))
        self.assertFalse(turbo._info_is_stranded("single\tR1\tHsa\t--ns\t\t\tdone"))
        self.assertFalse(turbo._info_is_stranded("\t\tdone"))


class SplitTests(Base):
    def test_every_record_lands_in_exactly_one_chunk_in_balance(self):
        out, record = self.split(chunks=4, count=4001)
        seen = []
        for k in range(4):
            seen += names(out / "chunks" / "c{}_R1.fastq.gz".format(k))
        self.assertEqual(sorted(seen), sorted("r{}/1".format(i) for i in range(4001)))
        self.assertEqual(sum(record["reads_per_chunk"]), 4001)
        self.assertLessEqual(max(record["reads_per_chunk"]) - min(record["reads_per_chunk"]), 1)
        self.assertEqual(json.loads((out / "split.json").read_text())["reads"], 4001)

    def test_mates_stay_together(self):
        out, record = self.split(chunks=3, count=3000, mates=2)
        for k in range(3):
            one = [n.split("/")[0] for n in names(out / "chunks" / "c{}_R1.fastq.gz".format(k))]
            two = [n.split("/")[0] for n in names(out / "chunks" / "c{}_R2.fastq.gz".format(k))]
            self.assertEqual(one, two)
        self.assertEqual(record["mates"], 2)
        self.assertEqual(len(record["probe_minus_fraction"]), 2)

    def test_unequal_mates_are_refused(self):
        a, b = self.dir / "a.fastq.gz", self.dir / "b.fastq.gz"
        fastq(a, 3000)
        fastq(b, 2999, 2)
        with self.assertRaisesRegex(RuntimeError, "different numbers of reads"):
            turbo.split([str(a), str(b)], self.db, "hg19", 3, self.dir / "w", "run1")

    def test_too_few_reads_for_the_chunks_are_refused(self):
        with self.assertRaisesRegex(RuntimeError, "fewer than 1000 per chunk"):
            self.split(chunks=8, count=3000)

    def test_a_probe_that_maps_too_few_reads_is_refused(self):
        os.environ["MOCK_BT_MAPPED"] = "100"
        with self.assertRaisesRegex(RuntimeError, "mapped to the mRNA library"):
            self.split()

    def test_decision_is_recorded(self):
        os.environ["MOCK_FRACS"] = "0.9"
        _, record = self.split()
        self.assertTrue(record["strand_specific"])
        self.assertEqual(record["rc"], [True])
        self.assertEqual(record["probe_minus_fraction"], [0.9])

    def test_the_probe_uses_the_native_bowtie_command_and_index(self):
        # the probe is the first 2M lines of each mate, mapped as RunDBS_1.pl does
        self.assertEqual(turbo.PROBE_LINES, 2000000)
        self.assertEqual(turbo.NOT_SPECIFIC, (0.35, 0.65))
        self.assertEqual(turbo.MIN_PROBE_READS, 500)


class ChunkTests(Base):
    def test_every_chunk_gets_the_runs_strandedness(self):
        os.environ["MOCK_FRACS"] = "0.9,0.1"
        out, _ = self.split(chunks=3, mates=2)
        self.align_chunks(out)
        calls = self.vast_calls("align")
        self.assertEqual(len(calls), 3)
        for call in calls:
            self.assertIn("--rc1", call)
            self.assertIn("--nrc2", call)
            self.assertNotIn("--ns", call)
            self.assertEqual(call[call.index("--IR_version") + 1], "2")

    def test_unspecific_runs_force_ns_and_keep_the_chunk_names(self):
        out, _ = self.split(chunks=3)
        self.align_chunks(out)
        calls = self.vast_calls("align")
        self.assertTrue(all("--ns" in c for c in calls))
        self.assertEqual(sorted(c[c.index("--name") + 1] for c in calls), ["run1_c0", "run1_c1", "run1_c2"])

    def test_a_chunk_with_the_wrong_number_of_reads_is_refused(self):
        out, _ = self.split(chunks=3)
        fastq(out / "chunks" / "c1_R1.fastq.gz", 5)
        with self.assertRaisesRegex(RuntimeError, "split recorded"):
            self.align_chunks(out)

    def test_a_chunk_missing_a_file_is_refused(self):
        out, _ = self.split(chunks=3)
        os.environ["MOCK_SKIP"] = ".micX"
        with self.assertRaisesRegex(RuntimeError, "left no .micX"):
            self.align_chunks(out)

    def test_expr_is_passed_on(self):
        out, _ = self.split(chunks=3)
        files = [str(out / "chunks" / "c0_R1.fastq.gz")]
        turbo.chunk(self.bin, self.db, "hg19", out / "split.json", 0, files, "run1", self.dir / "t.tar",
                    self.dir / "t.json", expr=True)
        self.assertIn("--expr", self.vast_calls("align")[0])


class MergeTests(Base):
    def merged(self, **kwargs):
        out, record = self.split(chunks=3)
        tars, metas = self.align_chunks(out)
        target = self.dir / "cache"
        turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars, metas, "run1", target,
                    identity={"cache_id": "abc"}, **kwargs)
        return out, tars, metas, target

    def test_merge_writes_the_native_cache_layout(self):
        _, _, _, target = self.merged()
        with tarfile.open(target / "to_combine.tar") as archive:
            names_ = sorted(archive.getnames())
        for suffix in (".eej2", ".exskX", ".micX", ".MULTI3X", ".info", ".IR2", ".IR.summary_v2.txt"):
            self.assertIn("to_combine/run1" + suffix + ".gz", names_)
        self.assertFalse([n for n in names_ if "run1_c" in n], "chunk files must not be in the cache")
        cache = json.loads((target / "cache.json").read_text())
        self.assertEqual((cache["sample"], cache["tool"], cache["chunks"], cache["reads"]),
                         ("run1", "vasttools_turbo", 3, 3000))
        self.assertEqual(cache["cache_id"], "abc")
        self.assertEqual(sorted(cache["members"]), sorted(names_))

    def test_merge_pools_every_chunk_into_the_sample(self):
        self.merged()
        call = self.vast_calls("merge")[0]
        self.assertIn("--move_to_PARTS", call)
        self.assertEqual(call[call.index("--sp") + 1], "hg19")
        self.assertNotIn("--expr", call)

    def test_the_cache_extracts_the_way_compare_reads_it(self):
        _, _, _, target = self.merged()
        with tarfile.open(target / "to_combine.tar") as archive:
            member = archive.getmember("to_combine/run1.eej2.gz")
            text = gzip.GzipFile(fileobj=archive.extractfile(member)).read().decode()
        self.assertEqual(text, "merged\t3\n")

    def test_expr_merges_the_expression_files(self):
        out, _ = self.split(chunks=3)
        tars, metas = self.align_chunks(out)
        turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars, metas, "run1", self.dir / "c", expr=True)
        self.assertIn("--expr", self.vast_calls("merge")[0])

    def test_a_missing_chunk_is_refused(self):
        out, tars, metas, _ = self.merged()
        with self.assertRaisesRegex(RuntimeError, "3 chunks were split but 2"):
            turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars[:2], metas[:2], "run1", self.dir / "c2")

    def test_read_counts_that_do_not_add_up_are_refused(self):
        out, tars, metas, _ = self.merged()
        meta = json.loads(Path(metas[1]).read_text())
        meta["reads"] += 1
        Path(metas[1]).write_text(json.dumps(meta))
        with self.assertRaisesRegex(RuntimeError, "do not add up"):
            turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars, metas, "run1", self.dir / "c2")

    def test_a_chunk_with_the_wrong_strandedness_is_refused(self):
        out, tars, metas, _ = self.merged()
        meta = json.loads(Path(metas[2]).read_text())
        meta["strand_specific"] = True
        Path(metas[2]).write_text(json.dumps(meta))
        with self.assertRaisesRegex(RuntimeError, "chunk 2 was aligned strand-specific"):
            turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars, metas, "run1", self.dir / "c2")

    def test_an_empty_merged_eej2_is_refused(self):
        out, tars, metas, _ = self.merged()
        script = self.bin / "vast-tools"
        script.write_text(VAST.replace('"merged\\t%d\\n" % len(groups)', '""'))
        with self.assertRaisesRegex(RuntimeError, "empty"):
            turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars, metas, "run1", self.dir / "c2")

    def test_intron_retention_version_one_is_refused(self):
        out, tars, metas, _ = self.merged()
        with self.assertRaises(SystemExit):
            turbo.merge(self.bin, self.db, "hg19", out / "split.json", tars, metas, "run1", self.dir / "c2",
                        ir_version=1)


class ComplexityTests(Base):
    HEAD = "\t".join("c%d" % i for i in range(26)) + "\n"

    def exsk(self, path, rows):
        lines = [self.HEAD]
        for event, value in rows:
            fields = ["x"] * 26
            fields[3], fields[18] = event, value
            lines.append("\t".join(fields) + "\n")
        Path(path).write_text("".join(lines))

    def test_the_annotation_that_merge_overwrites_is_restored(self):
        merged, a, b = self.dir / "m", self.dir / "a", self.dir / "b"
        self.exsk(merged, [("E1", "S"), ("E2", "S"), ("E3", "S")])
        self.exsk(a, [("E1", "S"), ("E2", "C"), ("E3", "S")])
        self.exsk(b, [("E1", "S"), ("E2", "C"), ("E3", "S")])
        self.assertEqual(turbo.restore_exsk_complexity(merged, [a, b]), 1)
        rows = [line.split("\t") for line in merged.read_text().splitlines()[1:]]
        self.assertEqual([(r[3], r[18]) for r in rows], [("E1", "S"), ("E2", "C"), ("E3", "S")])
        self.assertEqual(merged.read_text().splitlines()[0], self.HEAD.strip())

    def test_chunks_that_disagree_are_refused(self):
        merged, a, b = self.dir / "m", self.dir / "a", self.dir / "b"
        self.exsk(merged, [("E1", "S")])
        self.exsk(a, [("E1", "S")])
        self.exsk(b, [("E1", "C")])
        with self.assertRaisesRegex(RuntimeError, "disagree on the complexity of E1"):
            turbo.restore_exsk_complexity(merged, [a, b])

    def test_short_lines_and_unknown_events_are_left_alone(self):
        merged, a = self.dir / "m", self.dir / "a"
        merged.write_text("h\nshort\n")
        self.exsk(a, [("E1", "C")])
        self.assertEqual(turbo.restore_exsk_complexity(merged, [a]), 0)
        self.assertEqual(merged.read_text(), "h\nshort\n")


class CliTests(Base):
    def test_split_chunk_merge_through_main(self):
        reads = self.reads(3000)
        out = self.dir / "w"
        base = ["--dbdir", str(self.db), "--species", "hg19", "--run-id", "run1"]
        turbo.main(["split"] + base + ["--reads", str(reads[0]), "--chunks", "3", "--out-dir", str(out)])
        tars, metas = [], []
        for k in range(3):
            tar, meta = str(self.dir / "k{}.tar".format(k)), str(self.dir / "k{}.json".format(k))
            turbo.main(["chunk", "--vast-tools", str(self.bin)] + base + [
                "--split-json", str(out / "split.json"), "--index", str(k), "--reads",
                str(out / "chunks" / "c{}_R1.fastq.gz".format(k)), "--out-tar", tar, "--out-meta", meta])
            tars.append(tar)
            metas.append(meta)
        turbo.main(["merge", "--vast-tools", str(self.bin)] + base + [
            "--split-json", str(out / "split.json"), "--tars"] + tars + ["--metas"] + metas + [
            "--out-dir", str(self.dir / "cache"), "--identity", json.dumps({"cache_id": "z"})])
        self.assertTrue((self.dir / "cache" / "to_combine.tar").exists())
        self.assertEqual(json.loads((self.dir / "cache" / "cache.json").read_text())["cache_id"], "z")


@unittest.skipUnless(all(os.environ.get(k) for k in ("VASTTOOLS_TEST_DIR", "VASTTOOLS_TEST_DB",
                                                      "VASTTOOLS_TEST_SPECIES", "VASTTOOLS_TEST_FASTQ")),
                     "set VASTTOOLS_TEST_DIR/_DB/_SPECIES/_FASTQ to compare against the real vast-tools")
class RealTurboEquivalenceTests(unittest.TestCase):
    """Native align against split + chunk + merge on real data (the FASTQ as single-end,
    and its two halves as paired-end), with 2 and 3 chunks.

    Single-end: the same events, identical raw EEJ counts, every number within vast-tools'
    own rounding (0.06), identical labels. Paired-end: Trim.pl numbers its 50-mers with one
    counter over the whole run, which fixes the order of a pair's fragments after sorting
    and so which of them vast-tools' adjacent-line logic counts once; chunks renumber, so a
    few events move by a read (measured: up to 3.6% of EEJ rows, by at most 2 reads). Totals stay within 0.1% on the EEJ counts."""

    CHUNKS = (2, 3)
    KEYS = {".eej2": (0, 1), ".micX": (1,), ".exskX": (3,), ".MULTI3X": (3,), ".IR2": (0,),
            ".IR.summary_v2.txt": (0,)}

    @classmethod
    def setUpClass(cls):
        from src import run_vasttools as vt
        cls.dir = Path(tempfile.mkdtemp(prefix="turbo-real-"))
        cls.species = os.environ["VASTTOOLS_TEST_SPECIES"]
        cls.db = os.environ["VASTTOOLS_TEST_DB"]
        cls.tool = os.environ["VASTTOOLS_TEST_DIR"]
        single = os.environ["VASTTOOLS_TEST_FASTQ"]
        with gzip.open(single, "rt") as stream:
            records = [line for line in stream]
        half = len(records) // 8 * 4
        for mate, part in ((1, records[:half]), (2, records[half:2 * half])):
            with gzip.open(cls.dir / "pe_{}.fq.gz".format(mate), "wt") as sink:
                sink.writelines(part)
        cls.layouts = {"se": [single], "pe": [str(cls.dir / "pe_1.fq.gz"), str(cls.dir / "pe_2.fq.gz")]}
        cls.native = {}
        for name, reads in cls.layouts.items():
            vt.align(cls.tool, cls.db, cls.species, reads, "run1", cls.dir / ("native_" + name), threads=4)
            cls.native[name] = cls.dir / ("native_" + name) / "to_combine.tar"

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.dir, True)

    @staticmethod
    def table(archive_path, suffix):
        with tarfile.open(archive_path) as archive:
            member = archive.getmember("to_combine/run1" + suffix + ".gz")
            return gzip.GzipFile(fileobj=archive.extractfile(member)).read().decode().splitlines()

    def turbo_cache(self, name, chunks):
        out = self.dir / "turbo_{}_{}".format(name, chunks)
        turbo.split(self.layouts[name], self.db, self.species, chunks, out / "w", "run1", threads=4)
        tars, metas = [], []
        for k in range(chunks):
            tar, meta = out / "k{}.tar".format(k), out / "k{}.json".format(k)
            files = sorted(str(p) for p in (out / "w" / "chunks").glob("c{}_R*.fastq.gz".format(k)))
            turbo.chunk(self.tool, self.db, self.species, out / "w" / "split.json", k, files, "run1", tar, meta,
                        threads=4)
            tars.append(str(tar))
            metas.append(str(meta))
        turbo.merge(self.tool, self.db, self.species, out / "w" / "split.json", tars, metas, "run1", out / "cache")
        return out / "cache" / "to_combine.tar"

    @staticmethod
    def numbers(field):
        try:
            return [float(part) for part in field.split("=")]
        except ValueError:
            return None

    def rows(self, lines, suffix):
        key = self.KEYS[suffix]
        return {tuple(f[i] for i in key): f for f in (l.split("\t") for l in lines) if len(f) > max(key)}

    def test_native_and_turbo_agree(self):
        for name in self.layouts:
            for chunks in self.CHUNKS:
                cache = self.turbo_cache(name, chunks)
                for suffix in self.KEYS:
                    with self.subTest(layout=name, chunks=chunks, table=suffix):
                        a = self.rows(self.table(self.native[name], suffix), suffix)
                        b = self.rows(self.table(cache, suffix), suffix)
                        if name == "se" or suffix not in (".eej2", ".IR2", ".IR.summary_v2.txt"):
                            self.assertEqual(set(a), set(b))
                        else:
                            self.assertGreaterEqual(len(set(a) & set(b)), 0.99 * len(a))
                        labels = moved = 0
                        for key in set(a) & set(b):
                            for x, y in zip(a[key], b[key]):
                                if x == y:
                                    continue
                                nx, ny = self.numbers(x), self.numbers(y)
                                if nx is None or ny is None or len(nx) != len(ny):
                                    labels += 1
                                elif max(abs(p - q) for p, q in zip(nx, ny)) > 0.06:
                                    moved += 1
                        if name == "se":
                            self.assertEqual((labels, moved), (0, 0))
                        else:
                            self.assertLessEqual(labels + moved, 0.06 * len(a), "too many paired-end differences")
                if True:
                    with self.subTest(layout=name, chunks=chunks, table="eej2 total"):
                        total = lambda path: sum(int(l.split("\t")[2]) for l in self.table(path, ".eej2"))
                        self.assertLessEqual(abs(total(cache) - total(self.native[name])),
                                             0.001 * total(self.native[name]))


if __name__ == "__main__":
    unittest.main()
