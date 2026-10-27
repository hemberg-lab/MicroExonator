import gzip
import json
import os
import platform
import subprocess
import tempfile
import unittest
from pathlib import Path

from src import run_circrna_discovery as rd
from tests import circrna_fixture as fixture

REPOSITORY = Path(__file__).resolve().parent.parent


def read_tsv(path):
    with gzip.open(path, "rt") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in stream]


def write_star_cache(folder, run_id, text):
    folder.mkdir(parents=True, exist_ok=True)
    with gzip.open(folder / "Chimeric.out.junction.gz", "wt") as stream:
        stream.write(text)
    with gzip.open(folder / "SJ.out.tab.gz", "wt") as stream:
        stream.write("")
    (folder / "Log.final.out").write_text("")
    (folder / "cache.json").write_text(json.dumps({"run_id": run_id}))


class UnitTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.fx = fixture.build(self.dir / "fx")

    def tearDown(self):
        self.tmp.cleanup()

    def test_refflat_from_gtf(self):
        out = self.dir / "ref.txt"
        self.assertEqual(rd.gtf_to_refflat(self.fx["gtf"], out), 2)
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0].split("\t"), ["G1", "G1.1", "chr1", "+", "200", "2800", "200", "2800", "4",
                                                "200,1000,1600,2600,", "400,1150,1800,2800,"])
        self.assertEqual(lines[1].split("\t")[3], "-")

    def test_star_header_and_comment_lines_are_dropped(self):
        path = self.dir / "j.gz"
        with gzip.open(path, "wt") as stream:
            stream.write("chr_donorA\tbrkpt_donorA\n" "chr1\t1801\t+\n" "# 2.7.11b STAR\n")
        self.assertEqual(list(rd.junction_lines(path)), ["chr1\t1801\t+\n"])

    def test_circexplorer2_rows_use_the_catalogue_convention(self):
        known = self.dir / "known.txt"
        known.write_text("chr1\t1000\t1800\tcircular_RNA/12\t0\t+\t1000\t1000\t0,0,0\t2\t150,200\t0,600\t12\tcircRNA"
                         "\tG1\tG1.1\t2,3\tchr1:400-1000|chr1:1800-2600\n")
        row = rd.read_circexplorer2(known)[0]
        self.assertEqual((row["circ_id"], row["reads"], row["gene"], row["strand"]), ("chr1:1001|1800", 12, "G1", "+"))

    def test_support_counts_replicates_and_catalogue(self):
        preflight = {"comparison_id": "c", "project_id": "p", "replicates": {"a": ["A1", "A2"], "b": ["B1"]},
                     "collapse": {"u1": "A1", "u2": "A2", "u3": "A2", "u4": "B1"}}
        caches = {}
        for unit, calls in (("u1", [("chr1:1001|1800", 3)]), ("u2", [("chr1:1001|1800", 1)]),
                            ("u3", [("chr1:1001|1800", 2), ("chr9:5|90", 4)]), ("u4", [])):
            folder = self.dir / unit
            folder.mkdir()
            (folder / "cache.json").write_text(json.dumps({"run_id": unit}))
            with gzip.open(folder / "calls.tsv.gz", "wt") as stream:
                stream.write("\t".join(rd.CALL_COLUMNS) + "\n")
                for cid, reads in calls:
                    chrom, rest = cid.split(":")
                    start, end = rest.split("|")
                    stream.write("\t".join([cid, chrom, start, end, "+", str(reads), "G", "circRNA", "2"]) + "\n")
            caches[unit] = str(folder)
        catalogue = self.dir / "catalogue.bed"
        catalogue.write_text("chr1\t1001\t1800\tchr1:1001|1800\t.\t+\n")
        rd.support(preflight, caches, self.dir / "out", catalogue=catalogue, selection={"role": "test"})
        support = {row["circ_id"]: row for row in read_tsv(self.dir / "out" / "support.tsv.gz")}
        self.assertEqual((support["chr1:1001|1800"]["A1"], support["chr1:1001|1800"]["A2"]), ("3", "3"))
        self.assertEqual(support["chr1:1001|1800"]["catalogue_supported"], "yes")
        self.assertEqual(support["chr9:5|90"]["catalogue_supported"], "no")
        results = {row["feature_id"]: row for row in read_tsv(self.dir / "out" / "results.tsv.gz")}
        self.assertEqual((results["chr1:1001|1800"]["n_a"], results["chr1:1001|1800"]["n_b"]), ("2", "0"))
        self.assertEqual(results["chr1:1001|1800"]["p_value"], "")
        status = json.loads((self.dir / "out" / "status.json").read_text())
        self.assertEqual((status["status"], status["selection"]["role"]), ("no_tests", "test"))


ENV = os.environ.get("CIRCRNA_DISCOVERY_ENV")


@unittest.skipUnless(ENV, "set CIRCRNA_DISCOVERY_ENV to a conda prefix made from envs/umbrella-circrna-discovery.yaml")
class NativeParserTests(unittest.TestCase):
    """CIRCexplorer2 and DCC themselves, on STAR-format junction caches (emulated lines)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.tmp.name)
        cls.env = dict(os.environ, PATH=str(Path(ENV) / "bin") + os.pathsep + os.environ["PATH"])
        cls.python = str(Path(ENV) / "bin" / "python")
        cls.fx = fixture.build(cls.dir / "fx")
        text, cls.truth = fixture.star_chimeric_lines()
        cls.caches = {}
        for unit in ("a1", "a2", "b1", "b2"):
            write_star_cache(cls.dir / "star" / unit, unit, text if unit != "b2" else text.splitlines(True)[0])
            cls.caches[unit] = cls.dir / "star" / unit
        subprocess.run([cls.python, "src/run_circrna_discovery.py", "refflat", "--gtf", str(cls.fx["gtf"]),
                        "--out", str(cls.dir / "ref.txt")], cwd=REPOSITORY, env=cls.env, check=True, capture_output=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_script(self, *args):
        result = subprocess.run([self.python, "src/run_circrna_discovery.py"] + [str(a) for a in args],
                                cwd=REPOSITORY, env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_circexplorer2_calls_both_circles_with_their_strands(self):
        out = self.dir / "cx2"
        self.run_script("circexplorer2", "--star-cache", self.caches["a1"], "--refflat", self.dir / "ref.txt",
                        "--genome", self.fx["genome"], "--out-dir", out, "--run-id", "a1", "--tmp-root", self.dir / "tmp")
        calls = {row["circ_id"]: row for row in read_tsv(out / "calls.tsv.gz")}
        self.assertEqual({cid: int(row["reads"]) for cid, row in calls.items()}, self.truth)
        self.assertEqual(calls["chr2:2201|2400"]["strand"], "-")
        self.assertEqual(calls["chr1:1001|1800"]["gene"], "G1")

    def test_dcc_detects_filters_by_replicate_and_matches_coordinates(self):
        preflight = {"comparison_id": "c", "project_id": "p", "replicates": {"a": ["A1", "A2"], "b": ["B1", "B2"]},
                     "collapse": {"a1": "A1", "a2": "A2", "b1": "B1", "b2": "B2"}}
        (self.dir / "preflight.json").write_text(json.dumps(preflight))
        out = self.dir / "dcc"
        catalogue = self.dir / "catalogue.bed"
        catalogue.write_text("chr1\t1001\t1800\tchr1:1001|1800\t.\t+\n")
        args = ["dcc", "--preflight", self.dir / "preflight.json", "--out-dir", out, "--gtf", self.fx["gtf"],
                "--catalogue", catalogue, "--tmp-root", self.dir / "tmp", "--settings", json.dumps({"min_reads": 2,
                                                                                                    "min_replicates": 2})]
        for unit, folder in self.caches.items():
            args += ["--cache", "{}={}".format(unit, folder)]
        self.run_script(*args)
        support = {row["circ_id"]: row for row in read_tsv(out / "support.tsv.gz")}
        self.assertEqual(set(support), set(self.truth))
        self.assertEqual(support["chr1:1001|1800"]["A1"], str(self.truth["chr1:1001|1800"]))
        self.assertEqual(support["chr2:2201|2400"]["strand"], "-")
        self.assertEqual(support["chr1:1001|1800"]["catalogue_supported"], "yes")
        self.assertEqual(support["chr2:2201|2400"]["catalogue_supported"], "no")
        self.assertEqual(json.loads((out / "status.json").read_text())["status"], "no_tests")


class StarBuildTests(unittest.TestCase):
    """bioconda's STAR dispatch script crashed on the cluster (Illegal instruction, exit 132)."""

    def builds(self, folder, crash=(), flags="fpu sse sse2 pni ssse3 sse4_1 sse4_2 avx avx2"):
        for name in [n for n, _ in rd.STAR_BUILDS] + ["STAR"]:
            path = folder / name
            # exit 132 is what the shell reports for SIGILL
            path.write_text("#!/bin/sh\n" + ("exit 132\n" if name in crash else "echo 2.7.11b\n"))
            path.chmod(0o755)
        cpuinfo = folder / "cpuinfo"
        cpuinfo.write_text("processor\t: 0\nflags\t\t: {}\n".format(flags))
        return cpuinfo

    def pick(self, **kwargs):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            cpuinfo = self.builds(folder, **kwargs)
            old = os.environ.get("PATH", "")
            os.environ["PATH"] = str(folder) + os.pathsep + old
            try:
                return Path(rd.star_binary(str(cpuinfo))).name
            finally:
                os.environ["PATH"] = old

    def test_fastest_supported_build_is_used(self):
        self.assertEqual(self.pick(), "STAR-avx2")

    def test_builds_the_cpu_lacks_are_never_tried(self):
        # rc06: sse4_2 and no avx/avx2
        self.assertEqual(self.pick(flags="fpu sse sse2 pni ssse3 sse4_1 sse4_2"), "STAR-sse4.1")

    def test_a_build_that_crashes_is_skipped(self):
        self.assertEqual(self.pick(crash=("STAR", "STAR-avx2", "STAR-avx")), "STAR-sse4.1")

    def test_no_runnable_build_stops_with_the_list_tried(self):
        crash = tuple(n for n, _ in rd.STAR_BUILDS) + ("STAR",)
        with self.assertRaises(SystemExit) as stop:
            self.pick(crash=crash, flags="sse sse2")
        self.assertIn("STAR-plain", str(stop.exception))

    def test_override_and_machines_without_cpuinfo(self):
        os.environ["CIRCRNA_STAR_BIN"] = "/opt/STAR"
        try:
            self.assertEqual(rd.star_binary("/nonexistent"), "/opt/STAR")
        finally:
            del os.environ["CIRCRNA_STAR_BIN"]
        self.assertEqual(rd.star_binary("/nonexistent/cpuinfo"), "STAR")

    def test_star_index_command(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            fake = folder / "STAR-fake"
            fake.write_text('#!/bin/sh\necho "$@" > "$(dirname "$0")/args"\n')
            fake.chmod(0o755)
            crash = folder / "STAR-crash"
            crash.write_text("#!/bin/sh\nexit 132\n")
            crash.chmod(0o755)
            os.environ["CIRCRNA_STAR_BIN"] = str(fake)
            try:
                rd.main(["star-index", "--genome", "g.fa", "--gtf", "a.gtf", "--out-dir", str(folder / "index"),
                         "--threads", "8", "--overhang", "100", "--sa-bases", "14"])
                args = (folder / "args").read_text().split()
                self.assertEqual(args[:2], ["--runMode", "genomeGenerate"])
                for flag, value in (("--runThreadN", "8"), ("--genomeDir", str(folder / "index")),
                                    ("--sjdbGTFfile", "a.gtf"), ("--genomeSAindexNbases", "14"),
                                    ("--outFileNamePrefix", str(folder / "index") + "/")):
                    self.assertEqual(args[args.index(flag) + 1], value)
                os.environ["CIRCRNA_STAR_BIN"] = str(crash)
                with self.assertRaises(SystemExit) as stop:
                    rd.main(["star-index", "--genome", "g.fa", "--gtf", "a.gtf", "--out-dir", str(folder / "i2")])
                self.assertIn("exit 132", str(stop.exception))
            finally:
                del os.environ["CIRCRNA_STAR_BIN"]


@unittest.skipUnless(ENV and platform.system() == "Linux",
                     "real STAR runs on Linux only (macOS STAR cannot spawn readFilesCommand); run on the cluster")
class NativeStarTests(unittest.TestCase):
    """Cluster smoke test: real STAR -> CIRCexplorer2 on the fixture reads."""

    def test_star_then_circexplorer2_finds_the_simulated_circles(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            env = dict(os.environ, PATH=str(Path(ENV) / "bin") + os.pathsep + os.environ["PATH"])
            python = str(Path(ENV) / "bin" / "python")
            fx = fixture.build(work / "fx")
            index = work / "star_index"
            index.mkdir()
            subprocess.run([rd.star_binary(), "--runMode", "genomeGenerate", "--genomeDir", str(index), "--genomeFastaFiles",
                            str(fx["genome"]), "--sjdbGTFfile", str(fx["gtf"]), "--sjdbOverhang", "74",
                            "--genomeSAindexNbases", "5", "--outFileNamePrefix", str(index) + "/"],
                           env=env, check=True, capture_output=True)
            subprocess.run([python, "src/run_circrna_discovery.py", "star", "--reads", str(fx["r1"]), str(fx["r2"]),
                            "--index", str(index), "--out-dir", str(work / "star"), "--run-id", "r1",
                            "--tmp-root", str(work / "tmp")], cwd=REPOSITORY, env=env, check=True)
            record = json.loads((work / "star" / "cache.json").read_text())
            self.assertGreater(int(record["input_reads"]), 0)
            subprocess.run([python, "src/run_circrna_discovery.py", "refflat", "--gtf", str(fx["gtf"]),
                            "--out", str(work / "ref.txt")], cwd=REPOSITORY, env=env, check=True)
            subprocess.run([python, "src/run_circrna_discovery.py", "circexplorer2", "--star-cache", str(work / "star"),
                            "--refflat", str(work / "ref.txt"), "--genome", str(fx["genome"]), "--out-dir",
                            str(work / "cx2"), "--run-id", "r1", "--tmp-root", str(work / "tmp")],
                           cwd=REPOSITORY, env=env, check=True)
            calls = {row["circ_id"]: int(row["reads"]) for row in read_tsv(work / "cx2" / "calls.tsv.gz")}
            for cid, expected in fx["expected_bsj"].items():
                # STAR needs >= 15 nt on each side (chimSegmentMin): fewer than the 10-nt truth
                self.assertGreaterEqual(calls.get(cid, 0), int(0.5 * expected), (cid, calls, expected))
                self.assertLessEqual(calls.get(cid, 0), expected + 10)


if __name__ == "__main__":
    unittest.main()
