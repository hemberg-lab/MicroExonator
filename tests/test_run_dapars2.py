"""DaPars2 module helpers, and the storage admission test against real DaPars2.

The real test runs only with DAPARS2_DIR (the pinned DaPars2 src folder) and a
Python with numpy and scipy, e.g.
  DAPARS2_DIR=.../DaPars2/src <python with scipy> -m unittest tests.test_run_dapars2
It checks that coverage kept only inside the merged 3' UTR windows (gaps =
zero) gives DaPars2 exactly the same fitted sites and PDUI as whole-genome
coverage with explicit zero intervals.
"""

import csv
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src.run_dapars2 import (FORK_RUNNER, benjamini_hochberg, compare, overlapping_opposite_strands, paired_t,
                             sum_bedgraphs, welch)

CHROM_LENGTH = 20000
# two 3' UTRs: + strand 2000-5000 (proximal site ~3500), - strand 10000-13000 (~11500)
UTRS = [("chr1", 2000, 5000, "ENST1|GENEA|chr1|+", "+"), ("chr1", 10000, 13000, "ENST2|GENEB|chr1|-", "-")]


def coverage_profile(distal_level):
    """Per-base coverage: high on the short isoform part, `distal_level` beyond the site."""
    depth = [0] * CHROM_LENGTH
    for base in range(500, 1500):        # coverage outside any UTR
        depth[base] = 30
    for base in range(2000, 3500):
        depth[base] = 100
    for base in range(3500, 5000):
        depth[base] = distal_level
    for base in range(11500, 13000):
        depth[base] = 80
    for base in range(10000, 11500):
        depth[base] = distal_level
    return depth


def runs_of(depth, keep_zero):
    start = 0
    for base in range(1, CHROM_LENGTH + 1):
        if base == CHROM_LENGTH or depth[base] != depth[start]:
            if depth[start] or keep_zero:
                yield start, base, depth[start]
            start = base


class HelperTests(unittest.TestCase):
    def test_welch_and_bh(self):
        self.assertAlmostEqual(welch([0.1, 0.3, 0.25], [0.6, 0.7, 0.5]), 0.0100688814, places=8)
        self.assertIsNone(welch([0.1], [0.2, 0.3]))
        self.assertIsNone(welch([0.2, 0.2], [0.2, 0.2]))
        self.assertEqual(benjamini_hochberg([0.01, None, 0.04, 0.03]), [0.03, None, 0.04, 0.04])

    def test_technical_runs_sum_and_gaps_stay_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, rows in (("a", ["chr1\t0\t10\t2\n", "chr2\t5\t8\t1\n"]), ("b", ["chr1\t5\t20\t3\n"])):
                with gzip.open(root / (name + ".gz"), "wt") as stream:
                    stream.writelines(rows)
            sum_bedgraphs([root / "a.gz", root / "b.gz"], root / "sum.bedgraph", ["chr1", "chr2"])
            self.assertEqual((root / "sum.bedgraph").read_text().splitlines(),
                             ["chr1\t0\t5\t2", "chr1\t5\t10\t5", "chr1\t10\t20\t3", "chr2\t5\t8\t1"])

    def test_opposite_strand_overlaps_are_flagged(self):
        with tempfile.TemporaryDirectory() as directory:
            bed = Path(directory) / "utr.bed"
            bed.write_text("chr1\t100\t300\tA|a|chr1|+\t0\t+\nchr1\t250\t400\tB|b|chr1|-\t0\t-\n"
                           "chr1\t900\t950\tC|c|chr1|+\t0\t+\n")
            self.assertEqual(overlapping_opposite_strands(bed), {"chr1:100-300", "chr1:250-400"})


class PairedDesignTests(unittest.TestCase):
    """compare() with a fake DaPars2 run: paired t-test on block differences for designed comparisons."""

    def run_compare(self, design):
        replicates = {"a": ["A1", "A2", "A3"], "b": ["C1", "C2", "C3"]}
        pdui = {"A1": 0.5, "A2": 0.75, "A3": 0.68, "C1": 0.4, "C2": 0.5, "C3": 0.5}
        preflight = {"project_id": "p", "comparison_id": "c", "inference_supported": True, "reasons": [],
                     "replicates": replicates, "collapse": {r: r for side in replicates.values() for r in side}}
        if design:
            preflight["design"] = {"type": "paired", "block": "subject_id",
                                   "replicate_block": {"A1": "S1", "C1": "S1", "A2": "S2", "C2": "S2",
                                                       "A3": "S3", "C3": "S3"}}
            preflight["blocks"] = [{"block": "S{}".format(i), "a": ["A{}".format(i)], "b": ["C{}".format(i)]}
                                   for i in (1, 2, 3)]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            caches = {}
            for replicate in pdui:
                cache = root / "cache" / replicate
                cache.mkdir(parents=True)
                with gzip.open(cache / "utr_coverage.bedgraph.gz", "wt") as stream:
                    stream.write("chr1\t2000\t5000\t10\n")
                (cache / "depth.json").write_text(json.dumps({"mapped_reads": 1000}))
                caches[replicate] = str(cache)
            (root / "utr.bed").write_text("chr1\t2000\t5000\tENST1|GENEA|chr1|+\t0\t+\n")

            def fake_dapars2(command, cwd, **kwargs):
                folder = Path(cwd) / "dapars2_out_chr1"
                folder.mkdir()
                order = replicates["a"] + replicates["b"]
                (folder / "dapars2_result_temp.chr1.txt").write_text(
                    "Gene\tfit_value\tPredicted_Proximal_APA\tLoci\t" + "\t".join(r + "_PDUI" for r in order)
                    + "\nENST1|GENEA|chr1|+\t3.2\t3500\tchr1:2000-5000\t"
                    + "\t".join(str(pdui[r]) for r in order) + "\n")
            from unittest import mock
            with mock.patch("src.run_dapars2.subprocess.run", side_effect=fake_dapars2):
                (root / "pf.json").write_text(json.dumps(preflight))
                compare(root / "pf.json", caches, root / "utr.bed", root, root / "out")
            with gzip.open(root / "out" / "results.tsv.gz", "rt") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            return rows, json.loads((root / "out" / "status.json").read_text())

    def test_paired_t_matches_reference(self):
        # scipy.stats.ttest_1samp([0.1, 0.25, 0.18], 0).pvalue
        self.assertAlmostEqual(paired_t([0.1, 0.25, 0.18]), 0.055226946144194836, places=8)
        self.assertIsNone(paired_t([0.1]))
        self.assertIsNone(paired_t([0.1, 0.1]))

    def test_designed_comparison_uses_paired_t(self):
        rows, status = self.run_compare(design=True)
        self.assertAlmostEqual(float(rows[0]["p_value"]), 0.055226946144194836, places=5)
        self.assertAlmostEqual(float(rows[0]["effect"]), (0.1 + 0.25 + 0.18) / 3, places=5)
        self.assertIn("mean over blocks", rows[0]["effect_definition"])
        self.assertEqual(status["model"], "paired_test")
        unpaired_rows, unpaired_status = self.run_compare(design=False)
        self.assertAlmostEqual(float(unpaired_rows[0]["p_value"]),
                               welch([0.5, 0.75, 0.68], [0.4, 0.5, 0.5]), places=5)
        self.assertNotIn("model", unpaired_status)


@unittest.skipUnless(os.environ.get("DAPARS2_DIR"), "set DAPARS2_DIR to the pinned DaPars2 src folder")
class RealDaPars2Tests(unittest.TestCase):
    def test_utr_window_coverage_matches_whole_genome_coverage(self):
        dapars2 = Path(os.environ["DAPARS2_DIR"])
        samples = {"A1": 90, "A2": 85, "B1": 20, "B2": 25}      # A: long isoform, B: short
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed = root / "utr.bed"
            bed.write_text("".join("{}\t{}\t{}\t{}\t0\t{}\n".format(*utr) for utr in UTRS))
            windows = [(start, end) for _, start, end, _, _ in UTRS]
            # whole-genome coverage with explicit zeros, straight into DaPars2
            full = root / "full"
            full.mkdir()
            for sample, level in samples.items():
                with open(full / (sample + ".bedgraph"), "w") as stream:
                    for start, end, value in runs_of(coverage_profile(level), keep_zero=True):
                        stream.write("chr1\t{}\t{}\t{}\n".format(start, end, value))
            (full / "depth.txt").write_text("".join("{}.bedgraph\t1000000\n".format(s) for s in samples))
            (full / "chroms.txt").write_text("chr1\n")
            (full / "config.txt").write_text(
                "Annotated_3UTR={}\nAligned_Wig_files={}\nOutput_directory=out/\nOutput_result_file=r\n"
                "sequencing_depth_file=depth.txt\nNum_Threads=1\nCoverage_threshold=10\n".format(
                    bed, ",".join(s + ".bedgraph" for s in samples)))
            subprocess.run([sys.executable, "-c", FORK_RUNNER, str(dapars2 / "DaPars2_Multi_Sample_Multi_Chr.py"),
                            "config.txt", "chroms.txt"], cwd=full, check=True, capture_output=True)
            reference = (full / "out_chr1" / "r_result_temp.chr1.txt").read_text().splitlines()
            # the module: per-run caches clipped to the UTR windows, zeros dropped
            caches = {}
            for sample, level in samples.items():
                folder = root / "cache" / sample
                folder.mkdir(parents=True)
                with gzip.open(folder / "utr_coverage.bedgraph.gz", "wt") as stream:
                    for start, end, value in runs_of(coverage_profile(level), keep_zero=False):
                        for low, high in windows:
                            if start < high and end > low:
                                stream.write("chr1\t{}\t{}\t{}\n".format(max(start, low), min(end, high), value))
                (folder / "depth.json").write_text(json.dumps({"mapped_reads": 1000000}))
                caches[sample] = str(folder)
            preflight = root / "preflight.json"
            preflight.write_text(json.dumps({
                "project_id": "p", "comparison_id": "c", "inference_supported": True, "reasons": [],
                "replicates": {"a": ["A1", "A2"], "b": ["B1", "B2"]},
                "collapse": {s: s for s in samples}}))
            # as the workflow calls it: from the project root, every path relative
            import shutil
            shutil.copytree(dapars2, root / "software")
            previous = os.getcwd()
            os.chdir(root)
            self.addCleanup(os.chdir, previous)
            compare("preflight.json", {s: os.path.relpath(p, root) for s, p in caches.items()},
                    "utr.bed", "software", "module", threads=1, coverage_threshold=10)
            with gzip.open(root / "module" / "native.tsv.gz", "rt") as stream:
                module = stream.read().splitlines()
            self.assertEqual(sorted(module[1:]), sorted(reference[1:]))
            self.assertEqual(len(module), 3)
            with gzip.open(root / "module" / "results.tsv.gz", "rt") as stream:
                rows = [line.split("\t") for line in stream.read().splitlines()]
            header = rows[0]
            effects = {row[header.index("gene_name")]: float(row[header.index("effect")]) for row in rows[1:]}
            # A keeps more distal coverage: higher distal usage (PDUI) than B, on both strands
            self.assertGreater(effects["GENEA"], 0.3)
            self.assertGreater(effects["GENEB"], 0.3)


if __name__ == "__main__":
    unittest.main()
