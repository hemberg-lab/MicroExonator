import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from src import inventory_rmats_prep, reduce_featurecounts, reduce_qc


def featurecounts(path, bam, rows):
    path.write_text("# Program:featureCounts v2.0.6\n"
                    "Geneid\tChr\tStart\tEnd\tStrand\tLength\t{}\n".format(bam) +
                    "".join("{}\tchr1\t1\t10\t+\t{}\t{}\n".format(*row) for row in rows))
    return str(path)


class FeatureCountsReducerTests(unittest.TestCase):
    def test_merges_runs_in_sorted_order_and_rejects_mismatched_genes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts = {
                "run2": featurecounts(root / "b.txt", "b.bam", [("gA", 100, 3), ("gB", 50, 0)]),
                "run1": featurecounts(root / "a.txt", "a.bam", [("gA", 100, 7), ("gB", 50, 1)]),
            }
            prefix = str(root / "shard")
            reduce_featurecounts.reduce_group(counts, prefix, "ref", "hash")
            with gzip.open(prefix + ".featurecounts.tsv.gz", "rt") as stream:
                self.assertEqual(stream.read(), "gene_id\tlength\trun1\trun2\n"
                                 "gA\t100\t7\t3\ngB\t50\t1\t0\n")
            counts["run3"] = featurecounts(root / "c.txt", "c.bam", [("gB", 50, 1)])
            with self.assertRaisesRegex(ValueError, "differs"):
                reduce_featurecounts.shard_text(counts)


class QcReducerTests(unittest.TestCase):
    def test_one_row_of_raw_values_per_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "hisat2.txt").write_text("HISAT2 summary stats:\n\tOverall alignment rate: 92.50%\n")
            (root / "junctions.json").write_text(json.dumps({
                "inferred_strandedness": "firststrand", "sense_fraction": 0.03,
                "spliced_fragments": 40, "junction_fragments": 45, "ambiguous_strand": 2}))
            (root / "fc.summary").write_text("Status\tx.bam\nAssigned\t75\nUnassigned_NoFeatures\t25\n")
            run = {"run_id": "run1", "layout": "PE", "strandedness": "auto",
                   "hisat2_summary": str(root / "hisat2.txt"),
                   "junction_summary": str(root / "junctions.json"),
                   "featurecounts_summary": str(root / "fc.summary")}
            prefix = str(root / "shard")
            reduce_qc.reduce_group([run], prefix, "ref", "hash")
            with gzip.open(prefix + ".qc.tsv.gz", "rt") as stream:
                header, row = stream.read().splitlines()
            values = dict(zip(header.split("\t"), row.split("\t")))
            self.assertEqual(values["hisat2_overall_alignment"], "0.925")
            self.assertEqual(values["inferred_strandedness"], "firststrand")
            self.assertEqual(values["featurecounts_assigned_fraction"], "0.75")


class RmatsInventoryTests(unittest.TestCase):
    def test_keeps_one_prep_file_with_its_recorded_bam_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tmp = root / "tmp"
            tmp.mkdir()
            (tmp / "2026-09-28-10_00_00_1.rmats").write_text("bam\tumbrella/x.bam\n")
            destination = root / "keep" / "run1.rmats"
            inventory_rmats_prep.keep_run(tmp, destination, {
                "run_id": "run1", "bam": "umbrella/x.bam", "layout": "SE",
                "read_length": 100, "prep_file": str(destination)})
            self.assertEqual(destination.read_text(), "bam\tumbrella/x.bam\n")
            text = inventory_rmats_prep.inventory_text([str(destination) + ".json"])
            self.assertIn("run1\tumbrella/x.bam\tSE\t100", text)
            (tmp / "second.rmats").write_text("x")
            with self.assertRaisesRegex(ValueError, "expected one"):
                inventory_rmats_prep.single_prep_file(tmp)

    def test_entrypoints_load_without_pythonpath(self):
        for name in ("reduce_featurecounts", "reduce_qc", "inventory_rmats_prep"):
            script = Path(__file__).resolve().parents[1] / "src" / (name + ".py")
            result = subprocess.run([sys.executable, str(script), "--help"], text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(result.returncode, 0, result.stdout)


if __name__ == "__main__":
    unittest.main()
