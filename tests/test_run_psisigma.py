"""src/run_psisigma.py: annotation, intron-retention universe, per-run caches, comparisons.

Equivalence with native PSI-Sigma (NativeEquivalenceTests, NativeComparisonTests)
runs when PSISIGMA_TEST_DIR points at the pinned PSI-Sigma scripts and samtools,
perl and pysam are available; NativeComparisonTests also needs PSISIGMA_TEST_PERL5LIB
(the PDL, PDL-Stats, PDL-GSL and Statistics-Multtest build).
"""

import gzip
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "tests"))

from src import run_psisigma as ps  # noqa: E402
import psisigma_fixture as fx  # noqa: E402

NATIVE = os.environ.get("PSISIGMA_TEST_DIR")
TOOLS = bool(NATIVE and shutil.which("samtools") and shutil.which("perl")
             and importlib.util.find_spec("pysam"))


class AnnotationTests(unittest.TestCase):
    def test_exon_only_gtf_gains_gene_and_transcript_lines_sorted_for_psisigma(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            source = temp / "in.gtf"
            source.write_text(fx.gtf(gene_lines=False, transcript_lines=False)
                              + 'chr1_KI270706v1_random\tt\texon\t1\t50\t.\t+\t.\tgene_id "X"; transcript_id "TX";\n'
                              + 'chr1.alt\tt\texon\t1\t50\t.\t+\t.\tgene_id "Y"; transcript_id "TY";\n')
            summary = ps.prepare_gtf(source, temp / "out.gtf")
            self.assertEqual(summary["synthesized"], {"transcript": 3, "gene": 1})
            self.assertEqual(summary["dropped_contigs"], 2)
            lines = [line.split("\t") for line in (temp / "out.gtf").read_text().splitlines()]
            self.assertEqual({line[0] for line in lines}, {"chr1"})
            # gene, then transcripts, before the exons that start at the same base
            features = [line[2] for line in lines if line[3] == "101"]
            self.assertEqual(features[0], "gene")
            self.assertLess(max(i for i, f in enumerate(features) if f == "transcript"),
                            min(i for i, f in enumerate(features) if f == "exon"))
            t2 = next(line for line in lines if line[2] == "transcript" and 'transcript_id "T2"' in line[8])
            self.assertEqual((t2[3], t2[4]), ("101", "2200"))
            self.assertIn('transcript_type "protein_coding"', t2[8])
            ps.write_mapping(temp / "out.gtf", temp / "mapping.txt")
            self.assertIn("T1\tGENE1\t+", (temp / "mapping.txt").read_text())

    def test_existing_gene_and_transcript_lines_are_kept(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            (temp / "in.gtf").write_text(fx.gtf())
            summary = ps.prepare_gtf(temp / "in.gtf", temp / "out.gtf")
            self.assertEqual(summary["synthesized"], {"transcript": 0, "gene": 0})
            self.assertEqual(sum(1 for line in (temp / "out.gtf").read_text().splitlines()
                                 if line.split("\t")[2] == "transcript"), 3)

    def test_ir_positions_follow_psisigma_ir(self):
        with tempfile.TemporaryDirectory() as temp:
            database = Path(temp) / "x.db"
            # chr, i1s, i1e, i2s, i2e, tes, tee, anno, as, ae, name, gene
            database.write_text("chr1\t301\t301\t301\t1000\t0\t0\tx\t301\t1000\tchr1_101_2200_R_T2_1\tGENE1\n"
                                "chr1\t301\t301\t301\t1000\t0\t0\tx\t301\t1000\tchr1_101_2200_R_T1_1\tGENE1\n"
                                "chr1\t1\t2\t3\t4\t0\t0\tx\t1\t4\tchr1_1_4_W_T1_1\tGENE1\n"
                                "chr1\t5\t5\t5\t12\t0\t0\tx\t5\t12\tchr1_5_12_R_T3_1\tGENE1\n")
            universe = ps.ir_positions(database)
            self.assertEqual(universe[("chr1", "301", "1000")], [302, 439, 578, 717, 999])
            # a short intron: coincident positions collapse, as in PSI-Sigma's hash keys
            self.assertEqual(universe[("chr1", "5", "12")], [5, 6, 7, 11])
            self.assertEqual(len(universe), 2)


class CountingTests(unittest.TestCase):
    def test_read_rules_match_psisigma_ir(self):
        M, I, D, N, S, H = 0, 1, 2, 3, 4, 5
        self.assertTrue(ps.qualifies(99, 1, [(M, 50)], True))
        self.assertFalse(ps.qualifies(99, 2, [(M, 50)], True))          # multi-mapper
        self.assertFalse(ps.qualifies(355, 1, [(M, 50)], True))         # secondary
        self.assertFalse(ps.qualifies(97, 1, [(M, 50)], True))          # not a proper pair
        self.assertTrue(ps.qualifies(0, 1, [(M, 50)], False))           # single-end
        for op in (N, S, H):
            self.assertFalse(ps.qualifies(99, 1, [(M, 20), (op, 5), (M, 25)], True))
        self.assertTrue(ps.qualifies(99, 1, [(M, 20), (D, 2), (M, 28)], True))
        self.assertTrue(ps.qualifies(2048 | 99, 1, [(M, 50)], True))   # supplementary is not excluded
        self.assertEqual(list(ps.blocks(100, [(M, 20), (D, 2), (M, 10), (I, 3), (M, 5)])),
                         [(100, 120), (122, 132), (132, 137)])
        self.assertEqual(list(ps.blocks(100, [(M, 10), (7, 4), (M, 5)])), [(100, 110), (110, 115)])

    def test_position_counter_and_median(self):
        counter = ps.PositionCounter({("chr1", "1", "100"): [10, 20, 30], ("chr2", "1", "9"): [5]})
        counter.add("chr1", [(5, 21), (20, 31), (31, 40)])
        counter.add("chr1", [(10, 11)])
        counter.add("chr3", [(1, 100)])
        self.assertEqual(counter.counts(), {("chr1", 10): 2, ("chr1", 20): 2, ("chr1", 30): 1})
        self.assertEqual(ps.median_value([0, 9, 7, 6, 5]), 6)
        self.assertEqual(ps.median_value([1, 2, 4, 5]), 3)
        self.assertEqual(ps.median_value([1, 2]), 1)        # 1.5 truncated, as sprintf("%d")

    def test_canonical_motifs_on_both_strands(self):
        genome = {"chr1": "AAGTCCCCCCAGAA" + "AACTCCCCCCACAA" + "AACACCCCCCTTAA"}
        fetch = lambda chrom, start, end: genome[chrom][start:end]
        self.assertTrue(ps.canonical(fetch, "chr1", 3, 12))       # GT...AG
        self.assertTrue(ps.canonical(fetch, "chr1", 17, 26))      # CT...AC: GT-AG on the minus strand
        self.assertFalse(ps.canonical(fetch, "chr1", 31, 40))     # CA...TT

    def test_filter_keeps_sorted_seven_columns(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            (temp / "raw").write_text("chr1\t9\t20\t1\t0\t1\t0\nchr1\t3\t8\t2\t1\t3\t1\nbad\n")
            result = ps.filter_canonical(temp / "raw", temp / "out.gz", enabled=False)
            self.assertEqual(result["kept"], 2)
            with gzip.open(temp / "out.gz", "rt") as stream:
                self.assertEqual(stream.read().splitlines()[0], "chr1\t3\t8\t2\t1\t3\t1")


def write_cache(folder, run, junctions, positions, universe_path, extra=None):
    folder.mkdir(parents=True)
    with gzip.open(folder / "run.SJ.tab.gz", "wt") as stream:
        stream.write("".join("\t".join(map(str, row)) + "\n" for row in junctions))
    with gzip.open(folder / "ir_positions.tsv.gz", "wt") as stream:
        stream.write("".join("{}\t{}\t{}\n".format(*row) for row in positions))
    record = dict(extra or {}, run_id=run, canonical_only=True, universe_sha256=ps._sha256(universe_path),
                  sha256={name: ps._sha256(folder / name) for name in ("run.SJ.tab.gz", "ir_positions.tsv.gz")})
    (folder / "cache.json").write_text(json.dumps(record))
    (folder / "qc.json").write_text("{}")


class ComparisonTests(unittest.TestCase):
    def test_replicates_sum_runs_like_a_merged_bam(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            universe = {("chr1", "301", "1000"): [302, 439, 578, 717, 999]}
            ps.write_universe(universe, temp / "u.gz")
            write_cache(temp / "x", "x", [("chr1", 301, 1000, 4, 1, 5, 2)], [("chr1", 439, 3), ("chr1", 578, 2)],
                        temp / "u.gz")
            write_cache(temp / "y", "y", [("chr1", 301, 1000, 6, 0, 7, 0), ("chr1", 1061, 2000, 2, 0, 2, 0)],
                        [("chr1", 439, 4), ("chr1", 717, 5), ("chr1", 999, 1)], temp / "u.gz")
            junctions, introns = ps.replicate_tables(["x", "y"], {"x": temp / "x", "y": temp / "y"}, universe)
            self.assertEqual(junctions, ["chr1\t301\t1000\t10\t1\t12\t2\n", "chr1\t1061\t2000\t2\t0\t2\t0\n"])
            # positions 302,439,578,717,999 -> 0,7,2,5,1 -> median 2
            self.assertEqual(introns, ["chr1\t301\t1000\t2\n"])

    def test_caches_from_another_reference_or_setting_are_refused(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            ps.write_universe({("chr1", "1", "9"): [2]}, temp / "u.gz")
            write_cache(temp / "c1", "r", [], [], temp / "u.gz", {"reference_id": "ref1", "module_reference_id": "m1",
                                                               "cache_id": "c1"})
            record = json.loads((temp / "c1" / "cache.json").read_text())
            ids = {"reference_id": "ref1", "module_reference_id": "m1"}
            digest = ps._sha256(temp / "u.gz")
            ps.check_cache(record, temp / "c1", "r", ids, {"canonical_only": True}, digest)
            for key, value, settings in (("reference_id", "ref2", {}), ("module_reference_id", "m2", {}),
                                         ("cache_id", "other", {}), (None, None, {"canonical_only": False})):
                changed = dict(record, **({key: value} if key else {}))
                with self.subTest(key=key or "canonical_only"), self.assertRaises(ValueError):
                    ps.check_cache(changed, temp / "c1", "r", ids, dict({"canonical_only": True}, **settings), digest)
            with self.assertRaisesRegex(ValueError, "universe_sha256"):
                ps.check_cache(record, temp / "c1", "r", ids, {"canonical_only": True}, "0" * 64)
            with gzip.open(temp / "c1" / "run.SJ.tab.gz", "wt") as stream:
                stream.write("tampered\n")
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                ps.check_cache(record, temp / "c1", "r", ids, {"canonical_only": True}, digest)

    def test_unsupported_designs_run_nothing(self):
        with tempfile.TemporaryDirectory() as temp:
            preflight = {"inference_supported": True, "reasons": [], "project_id": "p", "comparison_id": "c",
                         "replicates": {"a": ["a1"], "b": ["b1", "b2"]}, "collapse": {}}
            with patch.object(ps.subprocess, "run") as native:
                ps.compare(preflight, {}, temp, temp, Path(temp) / "out")
                native.assert_not_called()
            status = json.loads((Path(temp) / "out" / "status.json").read_text())
            self.assertEqual(status["status"], "unsupported")
            self.assertIn("two biological replicates", status["reasons"][0])

    def test_normalized_effect_is_a_minus_b(self):
        with tempfile.TemporaryDirectory() as temp:
            table = Path(temp) / "x.sorted.txt"
            table.write_text("\t".join(["Event Region", "Gene Symbol", "Target Exon", "Event Type", "N", "T",
                                        "Exon Type", "Reference Transcript", "ΔPSI (%)", "T-test p-value",
                                        "FDR (BH)", "N Values", "T Values", "Database ID"]) + "\n"
                             + "1:300-2001\tGENE1\t1:1001-1060\tSES\t2\t2\tKnown\tT1\t55.5\t0.001\t0.004\t"
                               "20|25\t80|75\tchr1_300_2001_W_T1_1\n"
                             + "1:300-2001\tGENE1\t1:301-1000\tIR\t2\t2\tKnown\tT2\tnan\tnan\tnan\t1|1\t1|1\t"
                               "chr1_301_1000_R_T2_1\n")
            rows = ps.normalize(table, {}, {}, {"comparison_id": "c"}, {"T1": "+", "T2": "+"})
            first, second = rows
            self.assertEqual((first["effect"], first["p_value"], first["q_value"]), ("0.555", "0.001", "0.004"))
            self.assertEqual((first["chrom"], first["start"], first["end"], first["strand"]),
                             ("chr1", 1000, 1060, "+"))
            self.assertEqual((first["n_a"], first["n_b"], first["status"]), ("2", "2", "tested"))
            self.assertEqual((second["status"], second["effect"]), ("untested", ""))
            self.assertIn("Welch", first["statistic_type"])


class InstallerTests(unittest.TestCase):
    def test_scripts_that_are_not_the_pinned_release_are_refused(self):
        import io
        import tarfile
        from src import install_psisigma as installer
        with tempfile.TemporaryDirectory() as temp:
            archive = Path(temp) / "psisigma.tar.gz"
            with tarfile.open(archive, "w:gz") as source:
                for name in list(installer.SCRIPTS) + ["LICENSE"]:
                    data = b"not the pinned script\n"
                    info = tarfile.TarInfo("PSI-Sigma-x/" + name)
                    info.size = len(data)
                    source.addfile(info, io.BytesIO(data))
            with self.assertRaisesRegex(SystemExit, "not the pinned v2.3"):
                installer.install_scripts(Path(temp) / "root", archive)

    def test_cpan_tarballs_are_checked(self):
        from src import install_psisigma as installer
        with tempfile.TemporaryDirectory() as temp:
            for name, _, _ in installer.CPAN:
                (Path(temp) / name).write_text("tampered")
            with patch.object(installer.subprocess, "run") as cpanm, \
                    self.assertRaisesRegex(SystemExit, "pinned SHA-256"):
                installer.install_perl(Path(temp) / "root", tarballs=temp)
            cpanm.assert_not_called()


def samtools(*args):
    subprocess.run(["samtools"] + [str(a) for a in args], check=True)


def write_fixture(folder, samples=(("a1", .8, 1), ("a2", .75, 2), ("b1", .2, 3), ("b2", .25, 4))):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "genome.fa").write_text(">chr1\n" + fx.genome() + "\n")
    (folder / "ann.gtf").write_text(fx.gtf(gene_lines=False, transcript_lines=False))
    for name, psi, seed in samples:
        text = fx.sam(psi, seed)
        for suffix, body in (("", text), (".star", fx.mapq255(text))):
            (folder / "{}{}.sam".format(name, suffix)).write_text(body)
            samtools("view", "-b", "-o", folder / "{}{}.bam".format(name, suffix), folder / "{}{}.sam".format(name, suffix))
            samtools("index", folder / "{}{}.bam".format(name, suffix))


@unittest.skipUnless(TOOLS, "native PSI-Sigma scripts, samtools, perl and pysam needed")
class NativeEquivalenceTests(unittest.TestCase):
    """The per-run caches give exactly what PSI-Sigma's own steps give on the
    same alignments marked the STAR way (MAPQ 255 for unique)."""

    @classmethod
    def setUpClass(cls):
        cls.temp = Path(tempfile.mkdtemp())
        write_fixture(cls.temp)
        ps.reference(cls.temp / "ann.gtf", cls.temp / "genome.fa", NATIVE, cls.temp / "ref")
        cls.universe = ps.read_universe(cls.temp / "ref" / "ir_universe.tsv.gz")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp)

    def test_universe_is_the_retained_intron(self):
        self.assertEqual(list(self.universe), [("chr1", "301", "1000")])

    def test_junctions_equal_psisigma_generate_sj(self):
        for sample in ("a1", "b1"):
            ours = self.temp / "ours.{}.sj".format(sample)
            with open(self.temp / "log", "w") as log:
                ps.junctions(self.temp / (sample + ".bam"), Path(NATIVE) / "sjFromSAMcollapseUandM_inclOverlaps.awk",
                             ours, self.temp, 2, log)
            native = subprocess.run("samtools sort -n -O sam {} | awk -f {}".format(
                self.temp / (sample + ".star.bam"), Path(NATIVE) / "sjFromSAMcollapseUandM_inclOverlaps.awk"),
                shell=True, check=True, capture_output=True, text=True).stdout
            self.assertEqual(sorted(ours.read_text().splitlines()), sorted(native.splitlines()), sample)

    def test_intron_medians_equal_psisigma_ir(self):
        work = self.temp / "native_ir"
        work.mkdir()
        shutil.copyfile(self.temp / "ref" / "annotation.psisigma.gtf", work / "a.gtf")
        for label in ("u1", "u2"):
            (work / (label + ".SJ.out.tab")).write_text("")
        (work / "x.groupa.txt").write_text("u1\n")
        (work / "x.groupb.txt").write_text("u2\n")
        with open(work / "log", "w") as log:
            database = ps.build_database(NATIVE, work / "a.gtf", ["chr1"], "x.groupa.txt", "x.groupb.txt", "x",
                                         work, 1, log)
        for sample in ("a1", "a2", "b1", "b2"):
            os.symlink(self.temp / (sample + ".star.bam"), work / (sample + ".bam"))
            os.symlink(self.temp / (sample + ".star.bam.bai"), work / (sample + ".bam.bai"))
            subprocess.run(["perl", str(Path(NATIVE) / "PSIsigma-ir-v.1.2.pl"), database.name, sample + ".bam", "1"],
                           cwd=work, check=True, capture_output=True)
            native = (work / (sample + ".IR.out.tab")).read_text().splitlines()
            qc = dict.fromkeys(("records", "primary", "proper_pair", "primary_soft_clipped", "missing_nh",
                                "mapq60_not_unique", "unique_not_mapq60", "eq_x_cigar", "intron_reads"), 0)
            counts = ps.coverage(self.temp / (sample + ".bam"), self.universe, True, qc)
            ours = ["{}\t{}\t{}\t{}".format(chrom, i2s, i2e, ps.median_value([counts.get((chrom, p), 0) for p in positions]))
                    for (chrom, i2s, i2e), positions in self.universe.items()]
            self.assertEqual(ours, native, sample)
            self.assertEqual((qc["missing_nh"], qc["mapq60_not_unique"], qc["unique_not_mapq60"]), (0, 0, 0))
            self.assertGreater(qc["primary_soft_clipped"], 0)

    def test_prep_drops_the_non_canonical_junction(self):
        out = self.temp / "cache_a1"
        ps.prep(self.temp / "a1.bam", self.temp / "ref" / "ir_universe.tsv.gz", NATIVE, out, "a1", "PE",
                self.temp / "genome.fa", self.temp / "ref" / "genome.fa.fai", threads=2, tmp_root=self.temp)
        with gzip.open(out / "run.SJ.tab.gz", "rt") as stream:
            junctions = [line.split("\t")[:3] for line in stream.read().splitlines()]
        self.assertNotIn(["chr1", "1500", "1599"], junctions)
        self.assertIn(["chr1", "301", "1000"], junctions)
        qc = json.loads((out / "qc.json").read_text())
        self.assertEqual(qc["junctions"]["junctions"], 1)
        self.assertEqual(qc["junctions"]["unique_fragments"], 12)
        record = json.loads((out / "cache.json").read_text())
        self.assertTrue(record["canonical_only"])


@unittest.skipUnless(TOOLS and os.environ.get("PSISIGMA_TEST_PERL5LIB"), "PSI-Sigma's Perl modules needed")
class NativeComparisonTests(unittest.TestCase):
    def test_two_vs_two_through_psisigma(self):
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            write_fixture(temp)
            ps.reference(temp / "ann.gtf", temp / "genome.fa", NATIVE, temp / "ref")
            caches = {}
            for sample in ("a1", "a2", "b1", "b2"):
                caches[sample] = str(temp / "caches" / sample)
                ps.prep(temp / (sample + ".bam"), temp / "ref" / "ir_universe.tsv.gz", NATIVE, caches[sample],
                        sample, "PE", temp / "genome.fa", temp / "ref" / "genome.fa.fai", tmp_root=temp)
            preflight = {"inference_supported": True, "reasons": [], "project_id": "p", "comparison_id": "c",
                         "replicates": {"a": ["a1", "a2"], "b": ["b1", "b2"]},
                         "collapse": {s: s for s in caches}, "runs": {"a": ["a1", "a2"], "b": ["b1", "b2"]}}
            ps.compare(preflight, caches, temp / "ref", NATIVE, temp / "out", threads=1,
                       perl5lib=os.environ["PSISIGMA_TEST_PERL5LIB"])
            status = json.loads((temp / "out" / "status.json").read_text())
            self.assertEqual(status["status"], "ok", status)
            with gzip.open(temp / "out" / "results.tsv.gz", "rt") as stream:
                import csv
                rows = list(csv.DictReader(stream, delimiter="\t"))
            exon = [r for r in rows if r["event_class"] in ("SES", "MES") and r["start"] == "1000"]
            self.assertTrue(exon, rows)
            self.assertGreater(float(exon[0]["effect"]), 0.3)          # A includes E2 more than B
            intron = [r for r in rows if r["event_class"].startswith("IR")]
            self.assertTrue(intron, rows)
            self.assertGreater(float(intron[0]["effect"]), 0)          # A retains the intron more


if __name__ == "__main__":
    unittest.main()
