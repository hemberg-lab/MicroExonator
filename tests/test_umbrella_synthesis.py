import gzip
import json
import tempfile
import unittest
from pathlib import Path

from src.synthesize_umbrella import evaluate, qc_outliers, report
from src.umbrella_tool_inputs import leafcutter_files, suppa_tables


def gz(path, text):
    path.write_bytes(gzip.compress(text.encode()))
    return str(path)


class EvaluateTests(unittest.TestCase):
    def test_tool_formats_and_statuses(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deseq = gz(root / "d.tsv.gz", "gene_id\tbaseMean\tlog2FoldChange\tpadj\n"
                       "g1\t10\t1.5\t0.01\ng2\t10\t-2\t0.001\ng3\t10\t0.1\tNA\n")
            self.assertEqual(evaluate({"format": "deseq2", "paths": [deseq]})["called"], 2)
            delta = root / "delta.tsv"
            delta.write_text("exon_ID\tDeltaPsi\tProbability\n"
                             "a\t0.3\t0.95\nb\t-0.05\t0.99\nc\t-0.4\t0.8\n")
            result = evaluate({"format": "delta", "paths": [str(delta)]})
            self.assertEqual((result["tested"], result["called"], result["up"]), (3, 1, 1))
            (root / "SE.MATS.JC.txt").write_text("ID\tFDR\tIncLevelDifference\n1\t0.01\t-0.2\n2\t0.5\t0.3\n")
            (root / "A3SS.MATS.JC.txt").write_text("ID\tFDR\tIncLevelDifference\n1\t0.001\t0.15\n")
            result = evaluate({"format": "rmats", "paths": [str(root / "SE.MATS.JC.txt"),
                                                            str(root / "A3SS.MATS.JC.txt")]})
            self.assertEqual((result["called"], result["up"], result["down"]), (2, 1, 1))
            suppa = root / "SE.dpsi"
            suppa.write_text("b-a_dPSI\tb-a_p-val\nev1\t0.3\t0.01\nev2\t0.05\t0.001\n")
            self.assertEqual(evaluate({"format": "suppa", "paths": [str(suppa)]})["called"], 1)
            leaf = root / "cluster_significance.txt"
            leaf.write_text("cluster\tstatus\tp\tp.adjust\nc1\tSuccess\t0.001\t0.01\nc2\tSuccess\t0.3\tNA\n")
            result = evaluate({"format": "leafcutter", "paths": [str(leaf)]})
            self.assertEqual((result["called"], result["up"]), (1, None))
            unsupported = gz(root / "u.tsv.gz", "# inference not supported\n")
            self.assertEqual(evaluate({"format": "deseq2", "paths": [unsupported]})["status"], "unsupported")
            self.assertEqual(evaluate({"format": "deseq2", "paths": [str(root / "none")]})["status"], "missing")
            skipped = root / "skip.txt"
            skipped.write_text("# not configured: leafcutter_dir is not set\n")
            self.assertEqual(evaluate({"format": "leafcutter", "paths": [str(skipped)]})["status"],
                             "skipped: not configured: leafcutter_dir is not set")


class QcOutlierTests(unittest.TestCase):
    def test_outliers_are_computed_within_groups_and_strandedness_mismatch_is_flagged(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            header = ("run_id\tdeclared_strandedness\tinferred_strandedness\thisat2_overall_alignment\t"
                      "featurecounts_assigned_fraction\tjunction_fragments\n")
            rows = "".join("r{}\tunstranded\tunstranded\t{}\t0.7\t1000\n".format(i, rate)
                           for i, rate in enumerate((0.95, 0.94, 0.96, 0.95, 0.5)))
            rows += "s1\tfirststrand\tsecondstrand\t0.9\t0.7\t1000\n"
            qc = gz(root / "qc.tsv.gz", header + rows)
            groups = {"r{}".format(i): "g" for i in range(5)}
            groups["s1"] = "h"
            flagged = qc_outliers([qc], groups)
            self.assertIn(("g", "r4", "hisat2_overall_alignment"), [row[:3] for row in flagged])
            self.assertTrue(any("inferred secondstrand" in row[2] for row in flagged))


class ToolInputTests(unittest.TestCase):
    PREFLIGHT = {"replicates": {"a": ["A1", "A2"], "b": ["B1"]},
                 "collapse": {"a1": "A1", "a1b": "A1", "a2": "A2", "b1": "B1"}}

    def test_leafcutter_files_per_replicate(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            joined = gz(root / "j.tsv.gz",
                        "chrom\tstart\tend\tstrand\tlabel\tmax_anchor\tmulti_total\tshort_anchor_total\tA1\tA2\tB1\n"
                        "chr1\t201\t300\t+\tindex_intron\t40\t0\t0\t5\t0\t2\n")
            files = leafcutter_files(joined, self.PREFLIGHT, root / "leaf")
            self.assertEqual(Path(files[0]).read_text(), "chr1\t200\t300\t.\t5\t+\n")
            self.assertEqual(Path(files[1]).read_text(), "")
            self.assertEqual((root / "leaf" / "groups.txt").read_text(), "A1.junc\ta\nA2.junc\ta\nB1.junc\tb\n")

    def test_suppa_tables_average_technical_runs(self):
        with tempfile.TemporaryDirectory() as directory:
            joined = gz(Path(directory) / "t.tsv.gz",
                        "transcript_id\ta1\ta1b\ta2\tb1\ntx1\t2\t4\t5\t1\n")
            tables = suppa_tables(joined, self.PREFLIGHT)
            self.assertEqual(tables["a"], "A1\tA2\ntx1\t3\t5\n")
            self.assertEqual(tables["b"], "B1\ntx1\t1\n")


class ReportTests(unittest.TestCase):
    def test_report_lists_tools_capture_and_reasons_without_a_vote(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rates = root / "rates.json"
            rates.write_text(json.dumps({"a1": 0.97, "b1": None, "other": 0.5}))
            text = report({"comparison_id": "c", "group_a": "case", "group_b": "ctrl",
                           "inference_supported": False, "reasons": ["too few replicates"],
                           "tools": [{"name": "deseq2_featurecounts", "readout": "expression",
                                      "family": "annotation_based", "upstream": "hisat2",
                                      "format": "deseq2", "paths": [str(root / "missing")]}],
                           "capture_rates": [str(rates)], "runs": ["a1", "b1"], "qc": [], "groups": {}})
            self.assertIn("# reason: too few replicates", text)
            self.assertIn("deseq2_featurecounts\texpression\tannotation_based\thisat2\tmissing", text)
            self.assertIn("a1\t0.97", text)
            self.assertNotIn("other", text)


if __name__ == "__main__":
    unittest.main()
