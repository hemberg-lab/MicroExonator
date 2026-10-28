"""Native circRNA differential tests: CIRI_DE_replicate (CIRIquant) and CircTest on DCC counts.

The input-writing tests run anywhere (a fake Rscript). The native ones need:
  CIRCRNA_NATIVE_ENV    a conda prefix made from envs/umbrella-circrna-native-de.yaml
  CIRCRNA_NATIVE_ROOT   a folder made by `circrna_native_de.py software` (CIRI_DE.R, R/CircTest)
  CIRCRNA_DISCOVERY_ENV a prefix from envs/umbrella-circrna-discovery.yaml (DCC, pysam), for dcc-linear
"""
import csv
import gzip
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src import circrna_native_de as nd
from src.umbrella_modules import NORMALIZED

REPOSITORY = Path(__file__).resolve().parent.parent
NATIVE_ENV = os.environ.get("CIRCRNA_NATIVE_ENV")
NATIVE_ROOT = os.environ.get("CIRCRNA_NATIVE_ROOT")
DISCOVERY_ENV = os.environ.get("CIRCRNA_DISCOVERY_ENV")
PREFLIGHT = {"project_id": "p", "comparison_id": "c",
             "replicates": {"a": ["ra1", "ra2", "ra3"], "b": ["rb1", "rb2", "rb3"]},
             "collapse": {"ua1": "ra1", "ua2": "ra2", "ua3": "ra3", "ub1": "rb1", "ub2": "rb2", "ub3": "rb3"}}


def write_tsv(path, header, rows):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as stream:
        stream.write("\t".join(header) + "\n")
        for row in rows:
            stream.write("\t".join(str(value) for value in row) + "\n")


def ciriquant_source(folder, bsj):
    """A CIRIquant comparison folder as run_ciriquant.py compare leaves it."""
    folder = Path(folder)
    folder.mkdir(parents=True)
    replicates = PREFLIGHT["replicates"]["a"] + PREFLIGHT["replicates"]["b"]
    write_tsv(folder / "bsj_counts.tsv.gz", ["circ_id"] + replicates, [[cid] + values for cid, values in bsj])
    write_tsv(folder / "library.tsv", ["replicate", "side", "units", "total_reads", "mapped_reads"],
              [[r, "a" if r.startswith("ra") else "b", "u" + r[1:], 1000000, 900000] for r in replicates])
    rows = []
    for cid, _ in bsj:
        chrom, span = cid.split(":")
        start, end = span.split("|")
        rows.append({"feature_id": cid, "chrom": chrom, "start": start, "end": end, "strand": "+", "coverage": 10})
    write_tsv(folder / "results.tsv.gz", NORMALIZED, [[row.get(c, "") for c in NORMALIZED] for row in rows])
    return folder


def gene_shards(folder, genes):
    """Two featureCounts shards (one per group), counts per unit."""
    paths = []
    for group, units in (("case", ["ua1", "ua2", "ua3"]), ("control", ["ub1", "ub2", "ub3"])):
        path = Path(folder) / (group + ".featurecounts.tsv.gz")
        write_tsv(path, ["gene_id", "length"] + units, [[g, 1000] + [genes[g][u] for u in units] for g in genes])
        paths.append(str(path))
    return paths


FAKE_RSCRIPT = r'''#!{python}
import sys, json
args = dict(a[2:].split("=", 1) for a in sys.argv[2:])
open(args["out"].replace("circ_de.csv", "args.json"), "w").write(json.dumps(sys.argv[1:]))
lines = open(args["bsj"]).read().splitlines()[1:]
with open(args["out"], "w") as out:
    out.write(",logFC,logCPM,LR,PValue,DE,FDR\n")
    for i, line in enumerate(lines):
        out.write("{{}},1.5,3.2,12.0,0.0005,1,0.001\n".format(line.split(",")[0]))
open(args["out2"], "w").write(",logFC\n")
'''


class CiriquantInputTests(unittest.TestCase):
    def test_inputs_follow_prep_ciriquant_and_results_are_normalized(self):
        directory = Path(tempfile.mkdtemp())
        source = ciriquant_source(directory / "src", [
            ("chr1:100|200", [5, 6, 7, 1, 0, 2]),
            ("chr1:300|400", [0, 0, 0, 0, 0, 0]),              # no BSJ read: not in the matrix
            ("chr2:100|200", [3, "NA", 1, 1, 1, 1])])          # not indexed somewhere: not measured
        genes = {"G1": {u: 100 for u in PREFLIGHT["collapse"]}, "G2": {u: 7 for u in PREFLIGHT["collapse"]}}
        shards = gene_shards(directory, genes)
        rscript = directory / "Rscript"
        rscript.write_text(FAKE_RSCRIPT.format(python=sys.executable))
        rscript.chmod(rscript.stat().st_mode | stat.S_IEXEC)
        out = directory / "out"
        nd.ciriquant(PREFLIGHT, source, shards, "/sw/CIRI_DE.R", out, ids={"analysis_id": "a" * 16},
                     rscript=str(rscript))
        library = list(csv.DictReader(open(out / "library.csv")))
        self.assertEqual([(row["Sample"], row["Group"]) for row in library],
                         [("ra1", "T"), ("ra2", "T"), ("ra3", "T"), ("rb1", "C"), ("rb2", "C"), ("rb3", "C")])
        self.assertEqual(library[0]["Circular"], "10")                    # 2 x BSJ, as CIRIquant reports
        bsj = open(out / "circRNA_bsj.csv").read().splitlines()
        self.assertEqual(bsj, ["circ_id,ra1,ra2,ra3,rb1,rb2,rb3", "chr1:100|200,5,6,7,1,0,2"])
        gene = open(out / "gene_count_matrix.csv").read().splitlines()
        self.assertEqual(gene[0], "gene_id,ra1,ra2,ra3,rb1,rb2,rb3")
        self.assertEqual(gene[1], "G1,100,100,100,100,100,100")
        args = json.loads((out / "args.json").read_text())
        self.assertEqual(args[0], "/sw/CIRI_DE.R")
        with gzip.open(out / "results.tsv.gz", "rt") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["feature_id"], rows[0]["effect"], rows[0]["p_value"], rows[0]["q_value"],
                          rows[0]["chrom"], rows[0]["tool"]), ("chr1:100|200", "1.5", "0.0005", "0.001", "chr1",
                                                                "ciriquant"))
        self.assertEqual(json.loads((out / "status.json").read_text())["status"], "ok")

    def test_one_replicate_side_is_unsupported(self):
        directory = Path(tempfile.mkdtemp())
        preflight = dict(PREFLIGHT, replicates={"a": ["ra1"], "b": ["rb1", "rb2"]})
        nd.ciriquant(preflight, directory, [], "/sw/CIRI_DE.R", directory / "out")
        self.assertEqual(json.loads((directory / "out/status.json").read_text())["status"], "unsupported")

    def test_missing_runs_in_gene_shards_are_refused(self):
        directory = Path(tempfile.mkdtemp())
        shards = gene_shards(directory, {"G1": {u: 1 for u in PREFLIGHT["collapse"]}})
        with self.assertRaises(SystemExit):
            nd.gene_matrix(shards[:1], {"ra1": ["ua1"], "rb1": ["ub1"]})


def native_env():
    return dict(os.environ, PATH=str(Path(NATIVE_ENV) / "bin") + os.pathsep + os.environ["PATH"])


@unittest.skipUnless(NATIVE_ENV and NATIVE_ROOT, "set CIRCRNA_NATIVE_ENV and CIRCRNA_NATIVE_ROOT")
class NativeTests(unittest.TestCase):
    def run_native(self, *args):
        result = subprocess.run([str(Path(NATIVE_ENV) / "bin/python"), "src/circrna_native_de.py"] + list(args),
                                cwd=REPOSITORY, env=native_env(), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def rows(self, path):
        with gzip.open(path, "rt") as stream:
            return {row["feature_id"]: row for row in csv.DictReader(stream, delimiter="\t")}

    def test_ciri_de_replicate_finds_the_changed_circrna(self):
        import random
        rng = random.Random(1)
        bsj = [("chrX:1000|2000", [60, 70, 65, 10, 12, 8])]
        bsj += [("chr{}:{}|{}".format(1 + i % 5, 100 * i, 100 * i + 50),
                 [rng.randint(15, 25) for _ in range(6)]) for i in range(1, 60)]
        directory = Path(tempfile.mkdtemp())
        source = ciriquant_source(directory / "src", bsj)
        genes = {"G{}".format(g): {u: rng.randint(400, 600) for u in PREFLIGHT["collapse"]} for g in range(300)}
        shards = gene_shards(directory, genes)
        preflight = directory / "preflight.json"
        preflight.write_text(json.dumps(PREFLIGHT))
        out = directory / "out"
        self.run_native("ciriquant", "--preflight", str(preflight), "--source", str(source), "--gene-shard",
                        shards[0], "--gene-shard", shards[1], "--ciri-de", str(Path(NATIVE_ROOT) / "CIRI_DE.R"),
                        "--out-dir", str(out))
        rows = self.rows(out / "results.tsv.gz")
        self.assertEqual(len(rows), 60)
        changed = rows["chrX:1000|2000"]
        self.assertGreater(float(changed["effect"]), 2)                    # log2(A / B), A higher
        self.assertLess(float(changed["q_value"]), 0.01)
        others = [float(row["q_value"]) for cid, row in rows.items() if cid != "chrX:1000|2000"]
        self.assertGreater(min(others), 0.05)

    def circtest(self, n_a, n_b, circ_rows, linear_rows):
        directory = Path(tempfile.mkdtemp())
        dcc = directory / "dcc"
        dcc.mkdir()
        samples = ["s{:03d}.junction".format(i) for i in range(n_a + n_b)]
        write_tsv(dcc / "CircRNACount.tsv.gz", ["Chr", "Start", "End", "Strand"] + samples, circ_rows)
        write_tsv(directory / "LinearCount.tsv.gz", ["Chr", "Start", "End"] + samples, linear_rows)
        preflight = directory / "preflight.json"
        preflight.write_text(json.dumps({"replicates": {"a": ["a{}".format(i) for i in range(n_a)],
                                                        "b": ["b{}".format(i) for i in range(n_b)]}}))
        out = directory / "out"
        self.run_native("dcc-circtest", "--preflight", str(preflight), "--dcc-dir", str(dcc), "--linear",
                        str(directory / "LinearCount.tsv.gz"), "--circtest-lib", str(Path(NATIVE_ROOT) / "R"),
                        "--out-dir", str(out))
        return self.rows(out / "results.tsv.gz"), json.loads((out / "status.json").read_text())

    def test_circtest_separates_a_host_independent_change(self):
        circ = [["chr1", 100, 900, "+", 40, 42, 38, 8, 9, 7], ["chr2", 100, 900, "-", 20, 22, 18, 21, 19, 20],
                ["chr3", 100, 900, "+", 0, 1, 0, 0, 0, 1]]                  # too few reads: filtered
        linear = [["chr1", 100, 900, 60, 58, 62, 92, 91, 93], ["chr2", 100, 900, 80, 78, 82, 79, 81, 80],
                  ["chr3", 100, 900, 50, 50, 50, 50, 50, 50]]
        rows, status = self.circtest(3, 3, circ, linear)
        self.assertEqual(set(rows), {"chr1:100|900", "chr2:100|900"})
        self.assertLess(float(rows["chr1:100|900"]["q_value"]), 0.01)
        self.assertGreater(float(rows["chr1:100|900"]["effect"]), 0.25)    # ratio A minus B
        self.assertGreater(float(rows["chr2:100|900"]["p_value"]), 0.2)
        self.assertEqual(rows["chr2:100|900"]["strand"], "-")
        self.assertEqual((status["status"], status["tested"], status["filtered"]), ("ok", 2, 1))

    def test_circtest_with_unequal_groups(self):
        circ = [["chr1", 100, 900, "+", 40, 42, 38, 8, 9], ["chr2", 100, 900, "+", 20, 22, 18, 21, 19]]
        linear = [["chr1", 100, 900, 60, 58, 62, 92, 91], ["chr2", 100, 900, 80, 78, 82, 79, 81]]
        rows, status = self.circtest(3, 2, circ, linear)
        self.assertEqual(len(rows), 2)
        self.assertLess(float(rows["chr1:100|900"]["p_value"]), 0.01)
        self.assertIn("unequal groups", " ".join(status["reasons"]))


@unittest.skipUnless(DISCOVERY_ENV, "set CIRCRNA_DISCOVERY_ENV to run DCC's host-gene count")
class DccLinearTests(unittest.TestCase):
    def test_dcc_linear_counts_reads_at_circ_ends_as_dcc_does(self):
        directory = Path(tempfile.mkdtemp())
        genome = directory / "genome.fa"
        sequence = "ACGT" * 500
        genome.write_text(">chr1\n" + sequence + "\n")
        python = str(Path(DISCOVERY_ENV) / "bin/python")
        # unit ua1: 4 reads over the start (base 101), 2 over the end (base 900); ub1: 1 and 3
        layout = {"ua1": [(91, 4), (891, 2)], "ub1": [(91, 1), (891, 3)]}
        bams = {}
        for unit, blocks in layout.items():
            sam = directory / (unit + ".sam")
            lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:{}".format(len(sequence))]
            n = 0
            for position, copies in blocks:
                for _ in range(copies):
                    n += 1
                    read = sequence[position - 1:position + 19]
                    lines.append("\t".join(["r{}".format(n), "0", "chr1", str(position), "60", "20M", "*", "0", "0",
                                            read, "I" * 20]))
            sam.write_text("\n".join(lines) + "\n")
            bams[unit] = directory / (unit + ".bam")
            subprocess.run([python, "-c", "import pysam,sys; pysam.sort('-o', sys.argv[2], sys.argv[1]); "
                            "pysam.index(sys.argv[2]); pysam.faidx(sys.argv[3])", str(sam), str(bams[unit]),
                            str(genome)], check=True)
        dcc = directory / "dcc"
        dcc.mkdir()
        write_tsv(dcc / "CircRNACount.tsv.gz", ["Chr", "Start", "End", "Strand", "s000.junction", "s001.junction"],
                  [["chr1", 101, 900, "+", 5, 1]])
        preflight = directory / "preflight.json"
        preflight.write_text(json.dumps({"replicates": {"a": ["ra1"], "b": ["rb1"]},
                                         "collapse": {"ua1": "ra1", "ub1": "rb1"}}))
        out = directory / "LinearCount.tsv.gz"
        result = subprocess.run([python, "src/circrna_native_de.py", "dcc-linear", "--preflight", str(preflight),
                                 "--dcc-dir", str(dcc), "--bam", "ua1=" + str(bams["ua1"]), "--bam",
                                 "ub1=" + str(bams["ub1"]), "--genome", str(genome), "--out", str(out)],
                                cwd=REPOSITORY, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with gzip.open(out, "rt") as stream:
            lines = stream.read().splitlines()
        self.assertEqual(lines[0], "Chr\tStart\tEnd\ts000.junction\ts001.junction")
        self.assertEqual(lines[1], "chr1\t101\t900\t3\t2")                 # round((4 + 2) / 2), (1 + 3) / 2


if __name__ == "__main__":
    unittest.main()
