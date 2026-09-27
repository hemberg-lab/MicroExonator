import json
import tempfile
import unittest
from pathlib import Path

from src.comparison_preflight import load_comparisons, preflight
from src.umbrella_manifest import load_umbrella_manifest


HEADER = ("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\t"
          "source_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\n")


def manifest(root, rows):
    path = root / "manifest.tsv"
    lines = [HEADER]
    for run, replicate, batch, group, layout in rows:
        source_2 = "r2.fq.gz" if layout == "PE" else ""
        lines.append("{rep}\t{run}\t{rep}\tproj\t{batch}\t{group}\tfastq\tr1.fq.gz\t{s2}\t{layout}\t"
                     "unstranded\tref\ttrue\n".format(rep=replicate, run=run, batch=batch,
                                                       group=group, s2=source_2, layout=layout))
    path.write_text("".join(lines))
    return load_umbrella_manifest(path)


COMPARISON = {"comparison_id": "case_vs_ctrl", "project_id": "proj",
              "group_a": "case", "group_b": "ctrl", "exclude": []}


class PreflightTests(unittest.TestCase):
    def run_preflight(self, rows, comparison=COMPARISON):
        with tempfile.TemporaryDirectory() as directory:
            return preflight(manifest(Path(directory), rows), dict(comparison))

    def test_balanced_design_supports_inference_and_collapses_technical_runs(self):
        result = self.run_preflight([
            ("a1", "A1", "b1", "case", "SE"), ("a1b", "A1", "b1", "case", "SE"),
            ("a2", "A2", "b2", "case", "SE"),
            ("c1", "C1", "b1", "ctrl", "SE"), ("c2", "C2", "b2", "ctrl", "SE")])
        self.assertTrue(result["inference_supported"], result["reasons"])
        self.assertEqual(result["replicates"]["a"], ["A1", "A2"])
        self.assertEqual(result["collapse"]["a1b"], "A1")
        self.assertEqual(result["technical_run_replicates"], ["A1"])
        self.assertTrue(result["tools"]["rmats"]["supported"])
        self.assertEqual(result["shards"], ["ref/proj/case/b1", "ref/proj/case/b2",
                                            "ref/proj/ctrl/b1", "ref/proj/ctrl/b2"])

    def test_technical_runs_do_not_count_as_replicates(self):
        result = self.run_preflight([
            ("a1", "A1", "b1", "case", "SE"), ("a1b", "A1", "b1", "case", "SE"),
            ("c1", "C1", "b1", "ctrl", "SE"), ("c2", "C2", "b1", "ctrl", "SE")])
        self.assertFalse(result["inference_supported"])
        self.assertIn("1 biological replicate", result["reasons"][0])

    def test_confounding_with_batch_and_mixed_layout(self):
        result = self.run_preflight([
            ("a1", "A1", "b1", "case", "PE"), ("a2", "A2", "b1", "case", "PE"),
            ("c1", "C1", "b2", "ctrl", "SE"), ("c2", "C2", "b2", "ctrl", "SE")])
        self.assertFalse(result["inference_supported"])
        joined = " ".join(result["reasons"])
        self.assertIn("confounded with batch_id", joined)
        self.assertIn("confounded with layout", joined)
        self.assertFalse(result["tools"]["rmats"]["supported"])
        self.assertIn("mixed layouts", " ".join(result["tools"]["rmats"]["reasons"]))

    def test_exclusions_are_honoured(self):
        comparison = dict(COMPARISON, exclude=["A2"])
        result = self.run_preflight([
            ("a1", "A1", "b1", "case", "SE"), ("a2", "A2", "b1", "case", "SE"),
            ("a3", "A3", "b1", "case", "SE"),
            ("c1", "C1", "b1", "ctrl", "SE"), ("c2", "C2", "b1", "ctrl", "SE")], comparison)
        self.assertEqual(result["replicates"]["a"], ["A1", "A3"])
        self.assertEqual(result["excluded"], ["A2"])

    def test_comparisons_file_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "comparisons.yaml"
            path.write_text(json.dumps({"comparisons": [COMPARISON]}))
            self.assertEqual(load_comparisons(path)[0]["group_a"], "case")
            path.write_text(json.dumps({"comparisons": [dict(COMPARISON, group_b="case")]}))
            with self.assertRaisesRegex(ValueError, "the same"):
                load_comparisons(path)
            path.write_text(json.dumps({"comparisons": [COMPARISON, COMPARISON]}))
            with self.assertRaisesRegex(ValueError, "duplicate"):
                load_comparisons(path)


if __name__ == "__main__":
    unittest.main()
