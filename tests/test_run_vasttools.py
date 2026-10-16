"""src/run_vasttools.py: parsing, statistics, and align/compare against a mock vast-tools.

The real tool is exercised by RealVastToolsTests when VASTTOOLS_TEST_DIR (a
vast-tools checkout), VASTTOOLS_TEST_DB (a VASTDB folder holding e.g. Dme/) and
VASTTOOLS_TEST_SPECIES (e.g. dm6) are set and bowtie is on PATH.
"""

import gzip
import io
import json
import os
import random
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

from src import run_vasttools as vt  # noqa: E402

# A mock `vast-tools`: align writes to_combine files; merge, combine, tidy and
# diff write small tables from MOCK_TABLES (JSON: inclusion text, tidy rows, diff rows).
MOCK = textwrap.dedent('''\
    #!/usr/bin/env python3
    import json, os, sys
    from pathlib import Path
    args = sys.argv[1:]
    command = args[0]
    def value(flag):
        return args[args.index(flag) + 1]
    calls = Path(os.environ["MOCK_CALLS"])
    with open(calls, "a") as stream:
        stream.write(json.dumps(args) + "\\n")
    tables = json.loads(Path(os.environ["MOCK_TABLES"]).read_text()) if "MOCK_TABLES" in os.environ else {}
    if command == "align":
        out = Path(value("-o")) / "to_combine"
        out.mkdir(parents=True)
        name = value("--name")
        skip = os.environ.get("MOCK_SKIP", "")
        for suffix in (".eej2", ".exskX", ".micX", ".MULTI3X", ".IR2", ".IR.summary_v2.txt"):
            if suffix != skip:
                (out / (name + suffix)).write_text("x\\t1\\n")
        (out / (name + ".info")).write_text("paired\\tR1\\tR2\\tHs2\\tunstranded\\n")
        (Path(value("-o")) / (name + "-50-e.fa.gz")).write_text("reads")
    elif command == "merge":
        groups = [line.split("\\t") for line in Path(value("--groups")).read_text().splitlines()]
        out = Path(value("-o")) / "to_combine"
        for sub, sample in groups:
            assert (out / (sub + ".eej2")).exists(), "merge needs plain files"
            (out / (sample + ".eej2")).write_text("merged")
    elif command == "combine":
        (Path(value("-o")) / "INCLUSION_LEVELS_FULL-hg38-6.tab").write_text(tables["inclusion"])
    elif command == "tidy":
        Path(value("-outFile")).write_text(tables["tidy"])
    elif command == "diff":
        (Path(value("-o")) / (value("-d") + ".tab")).write_text(tables["diff"])
''')

HEADER = "GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX"


def inclusion(samples, rows):
    lines = [HEADER + "".join("\t{0}\t{0}-Q".format(s) for s in samples)]
    for gene, event, coord, length, full, complexity, psis in rows:
        cells = "".join("\t{}\tOK,OK,5=10,na,na@5,10".format(p) for p in psis)
        lines.append("\t".join([gene, event, coord, length, full, complexity]) + cells)
    return "\n".join(lines) + "\n"


class ParsingTests(unittest.TestCase):
    def test_sample_names_follow_vast_tools_rules(self):
        self.assertEqual(vt.sample_name("SRR123"), "SRR123")
        self.assertEqual(vt.sample_name("lib-1.a"), "lib_1_a")
        self.assertEqual(vt.sample_name("1abc"), "s_1abc")

    def test_strand_from_full_coordinates(self):
        # GFRA1-like microexon on minus: C1donor > C2acceptor
        self.assertEqual(vt.strand_of("HsaEX0027358", "chr10:116216000,116211631-116211645,116200000", "MIC"), "-")
        self.assertEqual(vt.strand_of("X", "chr2L:6045921,6045839+6045845-6045852,6045526", "MIC"), "-")
        self.assertEqual(vt.strand_of("HsaEX1", "chrX:100,200-300,400+410", "S"), "+")
        self.assertEqual(vt.strand_of("HsaALTD1-1/2", "chr1:200-300+250,400", "Alt5"), "+")
        self.assertEqual(vt.strand_of("HsaALTA1-1/2", "chr1:500,200-300+250", "Alt3"), "-")
        self.assertEqual(vt.strand_of("HsaINT1", "chr1:100-200=300-400:-", "IR"), "-")
        self.assertEqual(vt.strand_of("HsaEX1", "", "S"), "")

    def test_event_classes_come_from_the_complex_column(self):
        self.assertEqual(vt.event_class("HsaEX1", "MIC"), "microexon")
        self.assertEqual(vt.event_class("FBgn0038108_CASSETTE1", "S"), "cassette_exon")
        self.assertEqual(vt.event_class("HsaEX6000001", "ANN"), "cassette_exon")
        self.assertEqual(vt.event_class("HsaINT1", "IR"), "intron_retention")
        self.assertEqual(vt.event_class("HsaALTA1-1/2", "Alt3"), "alt_3ss")
        self.assertEqual(vt.event_class("HsaALTD1-1/2", "Alt5"), "alt_5ss")

    def test_mann_whitney_and_bh(self):
        self.assertAlmostEqual(vt.mann_whitney([1, 2, 3], [4, 5, 6]), 0.1)
        self.assertEqual(vt.mann_whitney([5, 5, 5], [5, 5, 5]), 1.0)
        self.assertIsNone(vt.mann_whitney([1], [2, 3]))
        # normal approximation for large groups: well separated groups are significant
        p = vt.mann_whitney(list(range(30)), list(range(100, 130)), exact_limit=0)
        self.assertLess(p, 1e-6)
        self.assertEqual(vt.benjamini_hochberg([0.01, None, 0.04, 0.03]), [0.03, None, 0.04, 0.04])


class MockToolTests(unittest.TestCase):
    def setUp(self):
        self.temp = Path(tempfile.mkdtemp())
        self.tool = self.temp / "vast-tools-src"
        self.tool.mkdir()
        script = self.tool / "vast-tools"
        script.write_text(MOCK)
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        os.environ["MOCK_CALLS"] = str(self.temp / "calls.jsonl")
        os.environ.pop("MOCK_SKIP", None)

    def tearDown(self):
        shutil.rmtree(self.temp)
        for key in ("MOCK_CALLS", "MOCK_TABLES", "MOCK_SKIP"):
            os.environ.pop(key, None)

    def calls(self):
        return [json.loads(line) for line in (self.temp / "calls.jsonl").read_text().splitlines()]

    def align(self, run_id, out, reads=("R1.fastq.gz", "R2.fastq.gz")):
        vt.align(self.tool, self.temp / "db", "hg38", [self.temp / r for r in reads], run_id, out,
                 threads=8, tmp_root=self.temp / "tmp", extra={"cache_id": "abc"})

    def test_align_keeps_only_gzipped_to_combine_files(self):
        out = self.temp / "cache"
        self.align("lib-1", out)
        args = self.calls()[0]
        self.assertEqual(args[:3], ["align", str((self.temp / "R1.fastq.gz").resolve()),
                                    str((self.temp / "R2.fastq.gz").resolve())])
        for flag, value in (("-sp", "hg38"), ("--name", "lib_1"), ("-c", "8"), ("--IR_version", "2")):
            self.assertEqual(args[args.index(flag) + 1], value)
        self.assertNotIn("--expr", args)
        with tarfile.open(out / "to_combine.tar") as archive:
            names = sorted(archive.getnames())
            self.assertIn("to_combine/lib_1.eej2.gz", names)
            self.assertTrue(all(n.endswith(".gz") for n in names))
            self.assertEqual(gzip.decompress(archive.extractfile("to_combine/lib_1.micX.gz").read()), b"x\t1\n")
        record = json.loads((out / "cache.json").read_text())
        self.assertEqual((record["sample"], record["run_id"], record["cache_id"]), ("lib_1", "lib-1", "abc"))
        self.assertIn("unstranded", record["info"])
        self.assertEqual(sorted(p.name for p in out.iterdir()), ["align.log", "cache.json", "to_combine.tar"])
        self.assertEqual(list((self.temp / "tmp").iterdir()), [])   # split reads are gone

    def test_align_refuses_missing_outputs(self):
        os.environ["MOCK_SKIP"] = ".micX"
        with self.assertRaisesRegex(RuntimeError, r"left no \.micX"):
            self.align("SRR1", self.temp / "cache")
        self.assertFalse((self.temp / "cache" / "to_combine.tar").exists())

    def write_tables(self, samples, rows, tidy_events, diff_rows):
        tidy = "EVENT\t" + "\t".join(samples) + "\n" + "".join(
            "{}\t{}\n".format(event, "\t".join(str(p) for p in psis))
            for event, psis in tidy_events)
        diff = "GENE\tEVENT\tA\tB\tE[dPsi]\tMV[dPsi]_at_0.95\n" + "".join(
            "{}\t{}\t{}\t{}\t{}\t{}\n".format(*row) for row in diff_rows)
        tables = self.temp / "tables.json"
        tables.write_text(json.dumps({"inclusion": inclusion(samples, rows), "tidy": tidy, "diff": diff}))
        os.environ["MOCK_TABLES"] = str(tables)

    def preflight(self, a_units, b_units, collapse=None, supported=True):
        replicates = {"a": sorted({(collapse or {}).get(u, u) for u in a_units}),
                      "b": sorted({(collapse or {}).get(u, u) for u in b_units})}
        path = self.temp / "preflight.json"
        path.write_text(json.dumps({
            "project_id": "P", "comparison_id": "case_vs_ctrl", "inference_supported": supported,
            "reasons": [] if supported else ["fewer than two biological replicates"],
            "runs": {"a": a_units, "b": b_units}, "replicates": replicates,
            "collapse": {u: (collapse or {}).get(u, u) for u in a_units + b_units}}))
        return path

    def caches(self, units):
        paths = {}
        for unit in units:
            out = self.temp / "caches" / unit
            self.align(unit, out)
            paths[unit] = str(out)
        return paths

    def test_compare_merges_technical_runs_and_maps_microexons(self):
        units_a, units_b = ["a1", "a2x", "a2y", "a3"], ["b1", "b2", "b3"]
        collapse = {"a2x": "a2", "a2y": "a2"}
        caches = self.caches(units_a + units_b)
        samples = ["a1", "rep_a2", "a3", "b1", "b2", "b3"]
        rows = [
            ("GFRA1", "HsaEX0027358", "chr10:116211631-116211645", "15",
             "chr10:116216000,116211631-116211645,116200000", "MIC", [80, 85, 82, 20, 25, 22]),
            ("GENE2", "HsaEX0000002", "chr1:1001-1010", "10", "chr1:900,1001-1010,1200", "MIC",
             [50, 50, 50, 52, 49, 51]),
            ("GENE3", "HsaEX0000003", "chr1:5001-5100", "100", "chr1:4000,5001-5100,6000", "S",
             ["NA", "NA", "NA", "NA", "NA", "NA"]),
            ("GENE4", "HsaINT0000004", "chr2:100-200", "101", "chr2:50-99=201-300:+", "IR",
             [10, 12, 11, 40, 45, 41]),
        ]
        self.write_tables(samples, rows,
                          tidy_events=[("HsaEX0027358", [80, 85, 82, 20, 25, 22]),
                                       ("HsaEX0000002", [50, 50, 50, 52, 49, 51])],
                          diff_rows=[("GFRA1", "HsaEX0027358", 0.82, 0.22, 0.60, 0.41),
                                     ("GENE2", "HsaEX0000002", 0.50, 0.51, -0.01, 0.0),
                                     ("GENE4", "HsaINT0000004", 0.11, 0.42, -0.31, 0.2)])
        me_db = self.temp / "me_db.bed"
        me_db.write_text("chr10\t116211630\t116211645\tGFRA1\t0\t-\n"     # matched, minus strand
                         "chr1\t1000\t1010\tME2\t0\t+\n"                     # matched, not changing
                         "chr1\t1000\t1010\tME2opp\t0\t-\n"                  # same exon, other strand
                         "chr3\t10\t20\tnovel\t0\t+\n")                      # not in VASTDB
        out = self.temp / "out"
        vt.compare(self.preflight(units_a, units_b, collapse), caches, self.tool, self.temp / "db", "hg38",
                   me_db, out, threads=4, ids={"tool_version": "v", "analysis_id": "id1"})
        commands = [c[0] for c in self.calls() if c[0] != "align"]
        self.assertEqual(commands, ["merge", "combine", "tidy", "diff"])
        merge = next(c for c in self.calls() if c[0] == "merge")
        self.assertIn("--move_to_PARTS", merge)
        self.assertEqual(merge[merge.index("--sp") + 1], "hg38")
        diff = next(c for c in self.calls() if c[0] == "diff")
        self.assertEqual(diff[diff.index("-a") + 1], "a1,rep_a2,a3")
        self.assertEqual(diff[diff.index("-b") + 1], "b1,b2,b3")
        for flag, value in (("-r", "0.95"), ("-m", "0.1"), ("-e", "10"), ("-S", "2")):
            self.assertEqual(diff[diff.index(flag) + 1], value)
        self.assertIn("--noPDF", diff)
        tidy = next(c for c in self.calls() if c[0] == "tidy")
        for flag in ("--noVLOW", "--p_IR", "--noB3"):
            self.assertIn(flag, tidy)
        self.assertEqual(tidy[tidy.index("--min_N") + 1], "2")
        status = json.loads((out / "status.json").read_text())
        self.assertEqual((status["status"], status["method"], status["called"]), ("ok", "diff", 1))
        self.assertEqual(status["merged_replicates"], ["rep_a2"])
        with gzip.open(out / "results.tsv.gz", "rt") as stream:
            results = {row["feature_id"]: row for row in __import__("csv").DictReader(stream, delimiter="\t")}
        gfra1 = results["HsaEX0027358"]
        self.assertEqual((gfra1["status"], gfra1["effect"], gfra1["native_statistic"]), ("called", "0.6", "0.41"))
        self.assertEqual((gfra1["chrom"], gfra1["start"], gfra1["end"], gfra1["strand"]),
                         ("chr10", "116211630", "116211645", "-"))
        self.assertEqual(gfra1["event_class"], "microexon")
        self.assertEqual(results["HsaEX0000002"]["status"], "tested")
        self.assertEqual(results["HsaINT0000004"]["status"], "filtered")   # diff estimate, but no tidy pass
        self.assertNotIn("HsaEX0000003", results)                            # never quantified
        with gzip.open(out / "microexons.tsv.gz", "rt") as stream:
            microexons = {row["microexon_id"]: row for row in __import__("csv").DictReader(stream, delimiter="\t")}
        self.assertEqual(microexons["chr10_-_116211630_116211645"]["vasttools_status"], "called")
        self.assertEqual(microexons["chr10_-_116211630_116211645"]["vastdb_event"], "HsaEX0027358")
        self.assertEqual(microexons["chr1_+_1000_1010"]["vasttools_status"], "tested")
        self.assertEqual(microexons["chr1_-_1000_1010"]["vasttools_status"], "not_in_vastdb")
        self.assertEqual(microexons["chr3_+_10_20"]["vasttools_status"], "not_in_vastdb")
        self.assertEqual(status["microexons"], {"called": 1, "tested": 1, "not_in_vastdb": 2})
        for name in vt.TABLES:
            self.assertTrue((out / name).exists(), name)
        self.assertEqual(sorted(p.name for p in out.iterdir() if p.is_dir()), [])   # private folder removed

    def test_many_replicates_switch_to_mann_whitney(self):
        units_a = ["a{}".format(i) for i in range(6)]
        units_b = ["b{}".format(i) for i in range(6)]
        caches = self.caches(units_a + units_b)
        samples = units_a + units_b
        psis = [80, 82, 84, 86, 88, 90, 10, 12, 14, 16, 18, 20]
        self.write_tables(samples, [("G", "HsaEX1", "chr1:11-20", "10", "chr1:1,11-20,40", "MIC", psis)],
                          tidy_events=[("HsaEX1", psis)], diff_rows=[("G", "HsaEX1", 0.85, 0.15, 0.7, 0.02)])
        me_db = self.temp / "me.bed"
        me_db.write_text("chr1\t10\t20\tme\t0\t+\n")
        out = self.temp / "out"
        vt.compare(self.preflight(units_a, units_b), caches, self.tool, self.temp / "db", "hg38", me_db, out)
        status = json.loads((out / "status.json").read_text())
        self.assertEqual((status["method"], status["called"]), ("mann_whitney", 1))
        with gzip.open(out / "results.tsv.gz", "rt") as stream:
            row = next(__import__("csv").DictReader(stream, delimiter="\t"))
        self.assertEqual(row["status"], "called")      # MV 0.02 would fail diff; MWU calls it
        self.assertAlmostEqual(float(row["effect"]), 0.7)
        self.assertLess(float(row["p_value"]), 0.01)

    def test_unsupported_design_runs_nothing(self):
        caches = self.caches(["a1", "b1"])
        out = self.temp / "out"
        vt.compare(self.preflight(["a1"], ["b1"], supported=False), caches, self.tool, self.temp / "db",
                   "hg38", self.temp / "none.bed", out)
        self.assertEqual([c[0] for c in self.calls()], ["align", "align"])
        status = json.loads((out / "status.json").read_text())
        self.assertEqual(status["status"], "unsupported")
        for name in vt.TABLES:
            with gzip.open(out / name, "rt") as stream:
                self.assertTrue(stream.read().startswith("# inference not supported"))


class InstallDbTests(unittest.TestCase):
    def test_install_checks_layout_and_records_the_library(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            archive = temp / "library.tar.gz"
            with tarfile.open(archive, "w:gz") as library:
                for name, text in (("Hs2/FILES/x.fa", ">x\nA\n"), ("Hs2/TEMPLATES/t.txt", "t"),
                                   ("Hs2/VASTDB_VERSION", "vastdb.hs2.23.06.20\n")):
                    data = text.encode()
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    library.addfile(info, io.BytesIO(data))
            vt.install_db(archive, temp, "hg38", "https://x/lib.tar.gz", temp / "installed.json")
            record = json.loads((temp / "installed.json").read_text())
            self.assertEqual((record["species_key"], record["vastdb_version"]), ("Hs2", "vastdb.hs2.23.06.20"))
            self.assertEqual(len(record["archive_sha256"]), 64)
            self.assertFalse(archive.exists())
            self.assertTrue((temp / "Hs2" / "FILES" / "x.fa").exists())

    def test_wrong_species_is_refused_before_extraction(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            archive = temp / "library.tar.gz"
            with tarfile.open(archive, "w:gz") as library:
                info = tarfile.TarInfo("Mm2/FILES/x")
                library.addfile(info, io.BytesIO(b""))
            with self.assertRaisesRegex(SystemExit, "has no Hs2/ folder"):
                vt.install_db(archive, temp / "db", "hg38", "", temp / "installed.json")
            self.assertTrue(archive.exists())
            self.assertFalse((temp / "installed.json").exists())


@unittest.skipUnless(os.environ.get("VASTTOOLS_TEST_DIR") and os.environ.get("VASTTOOLS_TEST_DB")
                     and shutil.which("bowtie"), "real vast-tools run: set VASTTOOLS_TEST_DIR/_DB/_SPECIES")
class RealVastToolsTests(unittest.TestCase):
    """align x7 (one replicate split in two runs) and compare with the real tool,
    on reads sampled from the library's own microexon and exon junctions."""

    def test_align_and_compare(self):
        tool, db = os.environ["VASTTOOLS_TEST_DIR"], os.environ["VASTTOOLS_TEST_DB"]
        species = os.environ.get("VASTTOOLS_TEST_SPECIES", "dm6")
        key = vt.SPECIES_KEY[species]
        files = Path(db) / key / "FILES"
        inc, exc = {}, {}
        for library in (files / "{}_MIC-50.fasta".format(key), files / "EXSK-50.fasta"):
            name = None
            for line in open(library):
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:]
                elif ".INC." in name:
                    inc.setdefault(name.split(".")[0], []).append(line)
                elif ".EXC." in name:
                    exc.setdefault(name.split(".")[0], []).append(line)
        events = sorted(e for e in inc if e in exc)
        rng = random.Random(1)
        changed = set(rng.sample(events, min(40, len(events))))
        complement = str.maketrans("ACGTN", "TGCAN")
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            caches = {}
            plan = {"a1": 0.85, "a2x": 0.85, "a2y": 0.85, "a3": 0.85, "b1": 0.15, "b2": 0.15, "b3": 0.15}
            for unit, psi in plan.items():
                folder = temp / "reads" / unit
                folder.mkdir(parents=True)
                with gzip.open(folder / "R1.fastq.gz", "wt") as r1, gzip.open(folder / "R2.fastq.gz", "wt") as r2:
                    n = 0
                    for event in events:
                        for _ in range(20 if unit.startswith("a2") else 40):
                            share = psi if event in changed else 0.5
                            seq = rng.choice(inc[event] if rng.random() < share else exc[event])
                            if len(seq) < 60:
                                continue
                            n += 1
                            mate = seq.translate(complement)[::-1]
                            r1.write("@r{0}/1\n{1}\n+\n{2}\n".format(n, seq, "I" * len(seq)))
                            r2.write("@r{0}/2\n{1}\n+\n{2}\n".format(n, mate, "I" * len(mate)))
                vt.align(tool, db, species, [folder / "R1.fastq.gz", folder / "R2.fastq.gz"], unit,
                         temp / "cache" / unit, threads=4, tmp_root=temp / "tmp")
                caches[unit] = str(temp / "cache" / unit)
            preflight = temp / "preflight.json"
            preflight.write_text(json.dumps({
                "project_id": "P", "comparison_id": "c", "inference_supported": True, "reasons": [],
                "runs": {"a": ["a1", "a2x", "a2y", "a3"], "b": ["b1", "b2", "b3"]},
                "replicates": {"a": ["a1", "a2", "a3"], "b": ["b1", "b2", "b3"]},
                "collapse": {"a1": "a1", "a2x": "a2", "a2y": "a2", "a3": "a3",
                             "b1": "b1", "b2": "b2", "b3": "b3"}}))
            me_db = temp / "me.bed"
            me_db.write_text("chr2L\t1000\t1012\tfake\t0\t+\n")
            out = temp / "out"
            vt.compare(preflight, caches, tool, db, species, me_db, out, threads=4)
            status = json.loads((out / "status.json").read_text())
            self.assertEqual(status["status"], "ok")
            self.assertEqual(status["merged_replicates"], ["rep_a2"])
            self.assertGreater(status["called"], 10)
            self.assertEqual(status["microexons"], {"not_in_vastdb": 1})


if __name__ == "__main__":
    unittest.main()
