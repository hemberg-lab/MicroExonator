import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src import run_ciriquant as rc
from tests import circrna_fixture as fixture

REPOSITORY = Path(__file__).resolve().parent.parent


def write_cache(folder, run_id, rows, mapped=1000, reference="m1"):
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / "counts.tsv.gz", "wt") as stream:
        stream.write("\t".join(rc.COUNT_COLUMNS) + "\n")
        for cid, bsj, fsj, measurement in rows:
            chrom = cid.split(":")[0]
            start, end = cid.split(":")[1].split("|")
            ratio = rc.junction_ratio(bsj, fsj) if measurement != "not_indexed" else None
            stream.write("\t".join(str(v) for v in (cid, chrom, start, end, "+", "" if bsj is None else bsj,
                                                       "" if fsj is None else fsj, "" if ratio is None else ratio,
                                                       measurement)) + "\n")
    (folder / "cache.json").write_text(json.dumps({"run_id": run_id, "total_reads": mapped * 2,
                                                   "mapped_reads": mapped, "module_reference_id": reference}))


def read_tsv(path):
    with gzip.open(path, "rt") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in stream]


class UnitTests(unittest.TestCase):
    def test_native_junction_ratio(self):
        self.assertEqual(rc.junction_ratio(3, 6), 0.5)
        self.assertIsNone(rc.junction_ratio(0, 0))
        self.assertEqual(rc.junction_ratio(0, 4), 0.0)

    def test_library_type_codes_and_auto(self):
        self.assertEqual(rc.library_type("unstranded"), (0, "unstranded"))
        self.assertEqual(rc.library_type("secondstrand"), (1, "secondstrand"))
        self.assertEqual(rc.library_type("firststrand"), (2, "firststrand"))
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / "s.json"
            summary.write_text(json.dumps({"inferred_strandedness": "firststrand"}))
            self.assertEqual(rc.library_type("auto", summary), (2, "firststrand"))
            summary.write_text(json.dumps({"inferred_strandedness": "undetermined"}))
            self.assertEqual(rc.library_type("auto", summary), (0, "unstranded"))
        with self.assertRaises(ValueError):
            rc.library_type("auto")

    def test_single_end_is_refused_before_running(self):
        with self.assertRaisesRegex(SystemExit, "paired-end"):
            rc.prep(["r1.fq.gz"], "x.bam", "c.bed", "g.fa", "out", "run", layout="SE")


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        circ, decoy, missing = "chr1:1001|1800", "chr1:201|1150", "chr9:1|50"
        write_cache(self.dir / "a1", "a1", [(circ, 6, 6, "measured"), (decoy, 0, 9, "measured"),
                                             (missing, None, None, "not_indexed")])
        write_cache(self.dir / "a2x", "a2x", [(circ, 2, 2, "measured"), (decoy, 0, 0, "no_support"),
                                               (missing, None, None, "not_indexed")])
        write_cache(self.dir / "a2y", "a2y", [(circ, 2, 4, "measured"), (decoy, 0, 3, "measured"),
                                               (missing, None, None, "not_indexed")])
        write_cache(self.dir / "b1", "b1", [(circ, 1, 8, "measured"), (decoy, 0, 5, "measured"),
                                             (missing, None, None, "not_indexed")])
        write_cache(self.dir / "b2", "b2", [(circ, 0, 0, "no_support"), (decoy, 0, 5, "measured"),
                                             (missing, None, None, "not_indexed")])
        self.preflight = {"comparison_id": "c1", "project_id": "p1",
                          "replicates": {"a": ["A1", "A2"], "b": ["B1", "B2"]},
                          "collapse": {"a1": "A1", "a2x": "A2", "a2y": "A2", "b1": "B1", "b2": "B2"}}
        self.caches = {unit: str(self.dir / unit) for unit in self.preflight["collapse"]}

    def tearDown(self):
        self.tmp.cleanup()

    def run_compare(self, preflight=None, out="out", **kwargs):
        rc.compare(preflight or self.preflight, self.caches, self.dir / out, ids={"analysis_id": "x"}, **kwargs)
        return self.dir / out

    def test_technical_units_are_summed_and_undefined_ratios_are_not_zero(self):
        out = self.run_compare()
        bsj = {row["circ_id"]: row for row in read_tsv(out / "bsj_counts.tsv.gz")}
        self.assertEqual(bsj["chr1:1001|1800"]["A2"], "4")         # 2 + 2 from two technical units
        self.assertEqual(bsj["chr9:1|50"]["A1"], "NA")              # not indexed stays missing
        rows = {row["feature_id"]: row for row in read_tsv(out / "results.tsv.gz")}
        circ = rows["chr1:1001|1800"]
        # A: 2*6/(12+6)=0.6667, 2*4/(8+6)=0.5714 ; B: 2/(2+8)=0.2 only (B2 has no support)
        self.assertAlmostEqual(float(circ["effect"]), (0.666667 + 0.571429) / 2 - 0.2, places=4)
        self.assertEqual((circ["n_a"], circ["n_b"]), ("2", "1"))
        self.assertEqual(circ["coverage"], "11")
        self.assertNotIn("chr1:201|1150", rows)                   # no BSJ anywhere: below min_bsj
        status = json.loads((out / "status.json").read_text())
        self.assertEqual(status["status"], "no_tests")
        self.assertIn("technical units summed within biological replicates", status["reasons"])

    def test_swapping_sides_reverses_the_effect(self):
        swapped = dict(self.preflight, replicates={"a": ["B1", "B2"], "b": ["A1", "A2"]})
        first = read_tsv(self.run_compare(out="o1") / "results.tsv.gz")[0]
        second = read_tsv(self.run_compare(preflight=swapped, out="o2") / "results.tsv.gz")[0]
        self.assertAlmostEqual(float(first["effect"]), -float(second["effect"]), places=6)

    def test_selection_metadata_reaches_status(self):
        out = self.run_compare(selection={"role": "polya_reference", "approved_by": "guillermo"})
        self.assertEqual(json.loads((out / "status.json").read_text())["selection"]["role"], "polya_reference")

    def test_caches_from_different_catalogues_are_refused(self):
        write_cache(self.dir / "b2", "b2", [("chr1:1001|1800", 1, 1, "measured")], reference="other")
        with self.assertRaisesRegex(SystemExit, "different catalogues"):
            self.run_compare()

    def test_empty_when_nothing_reaches_min_bsj(self):
        out = self.run_compare(settings={"min_bsj": 1000})
        self.assertEqual(json.loads((out / "status.json").read_text())["status"], "empty")


ENV = os.environ.get("CIRCRNA_CIRIQUANT_ENV")


@unittest.skipUnless(ENV, "set CIRCRNA_CIRIQUANT_ENV to a conda prefix made from envs/umbrella-ciriquant.yaml")
class NativeCiriquantTests(unittest.TestCase):
    """CIRIquant 1.1.3's own counting on a HISAT2 alignment of the synthetic fixture."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.tmp.name)
        cls.env = dict(os.environ, PATH=str(Path(ENV) / "bin") + os.pathsep + os.environ["PATH"])
        cls.python = str(Path(ENV) / "bin" / "python")
        cls.fx = fixture.build(cls.dir / "fx")
        catalogue_dir = cls.dir / "catalogue"
        subprocess.run([cls.python, "src/circrna_catalogue.py", "build", "--input", str(cls.fx["catalogue"]),
                        "--genome", str(cls.fx["genome"]), "--out-dir", str(catalogue_dir)],
                       cwd=REPOSITORY, env=cls.env, check=True, capture_output=True)
        cls.catalogue = catalogue_dir / "catalogue.bed"
        index = cls.dir / "hisat2" / "genome"
        index.parent.mkdir()
        subprocess.run(["hisat2-build", "-q", str(cls.fx["genome"]), str(index)], env=cls.env, check=True,
                       capture_output=True)
        cls.bam = cls.dir / "aligned.bam"
        subprocess.run("hisat2 -p 2 -x {} -1 {} -2 {} 2> {} | samtools sort -o {} - && samtools index {}".format(
            index, cls.fx["r1"], cls.fx["r2"], cls.dir / "hisat2.log", cls.bam, cls.bam), shell=True, env=cls.env,
            check=True)
        # as the module reference holds it: a link to the genome, with the .fai only beside the link
        link = cls.dir / "module_genome" / "genome.fa"
        link.parent.mkdir()
        link.symlink_to(cls.fx["genome"])
        Path(str(cls.fx["genome"]) + ".fai").rename(str(link) + ".fai")
        cls.genome_link = link
        cls.out = cls.dir / "cache"
        cls.result = subprocess.run(
            [cls.python, "src/run_ciriquant.py", "prep", "--reads", str(cls.fx["r1"]), str(cls.fx["r2"]),
             "--bam", str(cls.bam), "--catalogue", str(cls.catalogue), "--genome", str(link),
             "--out-dir", str(cls.out), "--run-id", "run1", "--threads", "2",
             "--identity", json.dumps({"module_reference_id": "test"}), "--tmp-root", str(cls.dir / "tmp")],
            cwd=REPOSITORY, env=cls.env, capture_output=True, text=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_prep_finishes(self):
        self.assertEqual(self.result.returncode, 0, self.result.stderr[-3000:])

    def test_bsj_counts_match_the_simulated_circles(self):
        rows = {row["circ_id"]: row for row in read_tsv(self.out / "counts.tsv.gz")}
        for cid, expected in self.fx["expected_bsj"].items():
            observed = int(rows[cid]["bsj"])
            # every fragment with a mate crossing the BSJ by >= 10 nt each side, give or take
            # the anchor and MAPQ filters CIRIquant applies
            self.assertGreaterEqual(observed, int(0.8 * expected), (cid, observed, expected))
            self.assertLessEqual(observed, expected + 10, (cid, observed, expected))
        self.assertEqual(rows["chr1:201|1150"]["bsj"], "0")        # decoy: linear reads only
        self.assertGreater(int(rows["chr1:1001|1800"]["fsj"]), 0)  # E1-E2 and E3-E4 linear reads
        self.assertEqual(rows["chr1:1001|1800"]["measurement"], "measured")

    def test_cache_records_provenance(self):
        record = json.loads((self.out / "cache.json").read_text())
        self.assertEqual(record["run_id"], "run1")
        self.assertEqual(record["patches"], list(rc.PATCHES))
        self.assertGreater(record["mapped_reads"], 0)
        self.assertEqual(record["catalogue_entries"], 3)

    def test_shared_index_gives_identical_counts(self):
        index_dir = self.dir / "shared_index"
        built = subprocess.run([self.python, "src/run_ciriquant.py", "index", "--catalogue", str(self.catalogue),
                                "--genome", str(self.genome_link), "--out-dir", str(index_dir), "--threads", "2"],
                               cwd=REPOSITORY, env=self.env, capture_output=True, text=True)
        self.assertEqual(built.returncode, 0, built.stderr[-2000:])
        out = self.dir / "cache_shared"
        result = subprocess.run(
            [self.python, "src/run_ciriquant.py", "prep", "--reads", str(self.fx["r1"]), str(self.fx["r2"]),
             "--bam", str(self.bam), "--catalogue", str(self.catalogue), "--genome", str(self.genome_link),
             "--out-dir", str(out), "--run-id", "run1", "--threads", "2", "--index-dir", str(index_dir),
             "--tmp-root", str(self.dir / "tmp")], cwd=REPOSITORY, env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr[-2000:])
        self.assertEqual(read_tsv(out / "counts.tsv.gz"), read_tsv(self.out / "counts.tsv.gz"))
        self.assertEqual(json.loads((out / "cache.json").read_text())["pseudo_index"], "shared")
        # an index built from another catalogue is refused
        other = self.dir / "other.bed"
        other.write_text(Path(self.catalogue).read_text().splitlines(True)[0])
        refused = subprocess.run(
            [self.python, "src/run_ciriquant.py", "prep", "--reads", str(self.fx["r1"]), str(self.fx["r2"]),
             "--bam", str(self.bam), "--catalogue", str(other), "--genome", str(self.genome_link),
             "--out-dir", str(self.dir / "refused"), "--run-id", "run1", "--index-dir", str(index_dir),
             "--tmp-root", str(self.dir / "tmp")], cwd=REPOSITORY, env=self.env, capture_output=True, text=True)
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("built from another catalogue", refused.stderr)

    def test_patched_helpers_are_still_broken_upstream(self):
        # if upstream fixes these, the patches can go: this test will then fail and say so
        probe = subprocess.run([self.python, "-c", (
            "from CIRIquant import circ; import itertools\n"
            "fai = circ.load_fai({!r})\n"
            "print(type(next(iter(fai.values()))[1]).__name__)\n"
            "try:\n    circ.grouper(range(3), 2.0); print('grouper-ok')\n"
            "except TypeError:\n    print('grouper-float-error')\n").format(str(self.genome_link) + ".fai")],
            env=self.env, capture_output=True, text=True, check=True)
        self.assertEqual(probe.stdout.split(), ["float", "grouper-float-error"])


if __name__ == "__main__":
    unittest.main()
