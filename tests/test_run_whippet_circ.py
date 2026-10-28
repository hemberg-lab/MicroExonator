"""Whippet --circ module: quant wrapper and the comparison (native diff and BS events)."""
import gzip
import json
import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from src import run_whippet_circ as wc

HEADER = ("Gene\tNode\tCoord\tStrand\tType\tPsi\tCI_Width\tCI_Lo,Hi\tTotal_Reads\tComplexity\tEntropy\t"
          "Inc_Paths\tExc_Paths\tEdges\n")


def linear(gene, node, psi, reads):
    return "{}\t{}\tchr1:{}-{}\t+\tCE\t{}\t0.1\t0.4,0.5\t{}\tK1\t0.0\t1-2-3:1\t1-3:0\t1-2:1\n".format(
        gene, node, 100 * node, 100 * node + 20, psi, reads)


def back_splice(gene, nodes, coord, psi, reads):
    return "{}\t{}\t{}\t+\tBS\t{}\t0.2\t0.1,0.3\t{}\tNA\tNA\tNA\tNA\tNA\n".format(gene, nodes, coord, psi, reads)


def write_psi(path, lines):
    with gzip.open(path, "wt") as stream:
        stream.write(HEADER + "".join(lines))


def write_cache(folder, run, lines, version="1.6.1"):
    os.makedirs(folder)
    write_psi(os.path.join(folder, "quant.psi.gz"), lines)
    Path(folder, "cache.json").write_text(json.dumps({"run_id": run, "module_reference_id": "m" * 16,
                                                      "whippet_version": version}))


def fake_native(julia, whippet_bin, files, prefix):
    with gzip.open(str(prefix) + ".diff.gz", "wt") as stream:
        stream.write(wc.DIFF_HEADER)
    Path(str(prefix) + ".args.json").write_text(json.dumps(files))


class CompareTests(unittest.TestCase):
    def run_compare(self, per_run, runs):
        directory = tempfile.mkdtemp()
        caches = {}
        for run, lines in per_run.items():
            caches[run] = os.path.join(directory, "cache", run)
            write_cache(caches[run], run, lines)
        preflight = {"project_id": "p", "comparison_id": "c", "runs": runs,
                     "collapse": {run: run for run in per_run}}
        out = os.path.join(directory, "out")
        wc.compare(preflight, caches, out, ids={"analysis_id": "a" * 16}, native=fake_native)
        return Path(out), caches

    def rows(self, out):
        with gzip.open(out / "results.tsv.gz", "rt") as stream:
            header = stream.readline().rstrip("\n").split("\t")
            return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in stream]

    def test_bs_events_are_matched_by_coordinates_not_line_order(self):
        bs1 = ("G1", "5-2", "chr1:500-200")
        bs2 = ("G2", "9-7", "chr2:900-700")
        a = [linear("G1", 1, 0.5, 20), back_splice(*bs1, 0.8, 50), back_splice(*bs2, 0.1, 40)]
        b = [linear("G1", 1, 0.5, 20), back_splice(*bs2, 0.1, 40), back_splice(*bs1, 0.2, 50)]
        out, caches = self.run_compare({"r1": a, "r2": a, "r3": b, "r4": b},
                                       {"a": ["r1", "r2"], "b": ["r3", "r4"]})
        rows = {row["gene_id"]: row for row in self.rows(out)}
        self.assertEqual(set(rows), {"G1", "G2"})
        self.assertGreater(float(rows["G1"]["effect"]), 0.5)
        self.assertGreater(float(rows["G1"]["native_statistic"]), 0.99)
        self.assertAlmostEqual(float(rows["G2"]["effect"]), 0.0, places=1)
        self.assertEqual(rows["G1"]["event_class"], "back_splice_junction")
        self.assertEqual((rows["G1"]["chrom"], rows["G1"]["start"], rows["G1"]["end"]), ("chr1", "200", "500"))
        self.assertEqual(rows["G1"]["tool_version"], "Whippet 1.6.1 --circ")
        status = json.loads((out / "status.json").read_text())
        self.assertEqual((status["status"], status["bs_events"], status["tested"]), ("ok", 2, 2))
        # the native diff gets the --circ psi files, A then B, in manifest order
        files = json.loads((out / "native.args.json").read_text())
        self.assertEqual([Path(path).parent.name for path in files["a"]], ["r1", "r2"])
        self.assertTrue(all(path.endswith("quant.psi.gz") for path in files["a"] + files["b"]))

    def test_a_sample_without_the_event_is_not_used_and_low_reads_are_skipped(self):
        bs = ("G1", "5-2", "chr1:500-200")
        out, _ = self.run_compare({"r1": [back_splice(*bs, 0.8, 50)], "r2": [linear("G1", 1, 0.5, 20)],
                                   "r3": [back_splice(*bs, 0.2, 3)]},
                                  {"a": ["r1", "r2"], "b": ["r3"]})
        self.assertEqual(self.rows(out), [])          # B has only 3 reads (min_reads 5)
        status = json.loads((out / "status.json").read_text())
        self.assertEqual(status["status"], "no_tests")
        with gzip.open(out / "bs_counts.tsv.gz", "rt") as stream:
            header, row = [line.rstrip("\n").split("\t") for line in stream]
        values = dict(zip(header, row))
        self.assertEqual((values["r1_bs_reads"], values["r2_psi"], values["r3_total_reads"]), ("40", "NA", "3"))

    def test_results_do_not_depend_on_the_other_events(self):
        bs = ("G1", "5-2", "chr1:500-200")
        extra = back_splice("G9", "3-1", "chr9:300-100", 0.5, 30)
        one, _ = self.run_compare({"r1": [back_splice(*bs, 0.7, 30)], "r2": [back_splice(*bs, 0.3, 30)]},
                                  {"a": ["r1"], "b": ["r2"]})
        two, _ = self.run_compare({"r1": [extra, back_splice(*bs, 0.7, 30)], "r2": [back_splice(*bs, 0.3, 30)]},
                                  {"a": ["r1"], "b": ["r2"]})
        pick = lambda out: [row for row in self.rows(out) if row["gene_id"] == "G1"][0]["native_statistic"]
        self.assertEqual(pick(one), pick(two))

    def test_caches_of_another_run_are_refused(self):
        directory = tempfile.mkdtemp()
        write_cache(os.path.join(directory, "x"), "other", [])
        write_cache(os.path.join(directory, "y"), "r2", [])
        with self.assertRaises(SystemExit):
            wc.compare({"runs": {"a": ["r1"], "b": ["r2"]}, "collapse": {}},
                       {"r1": os.path.join(directory, "x"), "r2": os.path.join(directory, "y")},
                       os.path.join(directory, "out"), native=fake_native)

    def test_a_side_without_runs_is_unsupported(self):
        directory = tempfile.mkdtemp()
        out = Path(directory, "out")
        wc.compare({"runs": {"a": ["r1"], "b": []}, "collapse": {}}, {}, out, native=fake_native)
        self.assertEqual(json.loads((out / "status.json").read_text())["status"], "unsupported")
        self.assertTrue((out / "native.diff.gz").is_file())


FAKE_JULIA = r'''#!{python}
import gzip, sys
args = sys.argv[2:]
assert "--circ" in args, args
prefix = args[args.index("-o") + 1]
open(prefix + ".args", "w").write(" ".join(args))
for suffix in ("gene.tpm.gz", "isoform.tpm.gz", "jnc.gz", "map.gz"):
    gzip.open(prefix + "." + suffix, "wt").write("x\n")
with gzip.open(prefix + ".psi.gz", "wt") as stream:
    stream.write({header!r})
    stream.write({linear!r})
    stream.write({bs!r})
'''


class NativeDeltaTests(unittest.TestCase):
    def test_one_run_sides_stay_file_lists(self):
        calls = []
        original = wc.subprocess.run
        wc.subprocess.run = lambda command, check: calls.append(command)
        try:
            wc.native_delta("julia", "/w/bin", {"a": ["/x/a.psi.gz"], "b": ["/x/b1.psi.gz", "/x/b2.psi.gz"]}, "/o/native")
        finally:
            wc.subprocess.run = original
        self.assertEqual(calls, [["julia", "/w/bin/whippet-delta.jl", "-a", "/x/a.psi.gz,",
                                  "-b", "/x/b1.psi.gz,/x/b2.psi.gz,", "-o", "/o/native"]])


class QuantTests(unittest.TestCase):
    def test_quant_adds_circ_keeps_flags_and_drops_the_read_map(self):
        directory = tempfile.mkdtemp()
        julia = os.path.join(directory, "julia")
        Path(julia).write_text(FAKE_JULIA.format(python=sys.executable, header=HEADER, linear=linear("G", 1, 0.5, 9),
                                                 bs=back_splice("G", "5-2", "chr1:500-200", 0.4, 10)))
        os.chmod(julia, os.stat(julia).st_mode | stat.S_IEXEC)
        bin_dir = os.path.join(directory, "Whippet.jl", "bin")
        os.makedirs(bin_dir)
        Path(directory, "Whippet.jl", "Project.toml").write_text('name = "Whippet"\nversion = "1.6.1"\n')
        out = os.path.join(directory, "cache")
        record = wc.quant(["r1.fastq.gz", "r2.fastq.gz"], "index.jls", out, "run_1", julia=julia, whippet_bin=bin_dir,
                          flags="--biascorrect", identity={"cache_id": "c" * 16})
        self.assertEqual(sorted(os.listdir(out)), ["cache.json", "quant.gene.tpm.gz", "quant.isoform.tpm.gz",
                                                   "quant.jnc.gz", "quant.log", "quant.psi.gz"])
        self.assertEqual((record["whippet_version"], record["events"], record["bs_events"]), ("1.6.1", 2, 1))
        self.assertIn("--biascorrect --circ", record["command"])
        self.assertIn("r1.fastq.gz r2.fastq.gz -x index.jls", record["command"])
        self.assertEqual(json.loads(Path(out, "cache.json").read_text())["cache_id"], "c" * 16)


if __name__ == "__main__":
    unittest.main()
