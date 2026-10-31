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


def metadata_manifest(root):
    path = root / "manifest.tsv"
    lines = [HEADER.rstrip("\n") + "\tsex\tphenotype\n"]
    for run, group, sex, phenotype in (("r1", "ctrl", "F", "control"), ("r2", "ctrl", "M", "control"),
                                       ("r3", "ctrl", "F", "control"), ("r4", "case", "F", "ASD"),
                                       ("r5", "case", "F", "ASD"), ("r6", "case", "M", "ASD")):
        lines.append("{r}\t{r}\t{r}\tproj\tb1\t{g}\tfastq\tr1.fq.gz\t\tSE\tunstranded\tref\ttrue"
                     "\t{s}\t{p}\n".format(r=run, g=group, s=sex, p=phenotype))
    path.write_text("".join(lines))
    return load_umbrella_manifest(path)


def comparisons_file(root, text):
    path = root / "comparisons.yaml"
    path.write_text(text)
    return load_comparisons(path)


class SelectionTests(unittest.TestCase):
    def test_metadata_selection_keeps_runs_in_their_shards(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            comparison = comparisons_file(root, (
                "comparisons:\n  - comparison_id: female_asd\n    project_id: proj\n"
                "    a: {label: asd_f, where: {phenotype: ASD, sex: F}}\n"
                "    b: {label: ctrl_f, groups: ctrl, where: {sex: [F]}}\n"))[0]
            result = preflight(metadata_manifest(root), comparison)
        self.assertEqual(result["group_a"], "asd_f")
        self.assertEqual(result["runs"], {"a": ["r4", "r5"], "b": ["r1", "r3"]})
        self.assertEqual(result["shards"], ["ref/proj/case/b1", "ref/proj/ctrl/b1"])
        self.assertTrue(result["inference_supported"], result["reasons"])

    def test_group_shorthand_still_works_and_samples_select(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            shorthand, by_sample = comparisons_file(root, (
                "comparisons:\n  - {comparison_id: old, project_id: proj, group_a: case, group_b: ctrl}\n"
                "  - comparison_id: picked\n    project_id: proj\n"
                "    a: {label: two_cases, samples: [r4, r6]}\n    group_b: ctrl\n"))
            manifest = metadata_manifest(root)
            self.assertEqual(preflight(manifest, shorthand)["runs"]["a"], ["r4", "r5", "r6"])
            self.assertEqual(preflight(manifest, by_sample)["runs"]["a"], ["r4", "r6"])

    def test_bad_selections_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = metadata_manifest(root)
            overlap = comparisons_file(root, (
                "comparisons:\n  - comparison_id: x\n    project_id: proj\n"
                "    a: {label: f, where: {sex: F}}\n    b: {label: asd, where: {phenotype: ASD}}\n"))[0]
            with self.assertRaisesRegex(ValueError, "both sides: r4, r5"):
                preflight(manifest, overlap)
            unknown = comparisons_file(root, (
                "comparisons:\n  - comparison_id: x\n    project_id: proj\n"
                "    a: {label: f, where: {tissue: brain}}\n    group_b: ctrl\n"))[0]
            with self.assertRaisesRegex(ValueError, "no manifest column tissue"):
                preflight(manifest, unknown)
            for text, message in (
                    ("    a: {label: f}\n    group_b: ctrl\n", "selects nothing"),
                    ("    a: {label: f, colour: red}\n    group_b: ctrl\n", "unknown selection key"),
                    ("    group_a: case\n    a: {label: f, groups: case}\n    group_b: ctrl\n", "either group_a")):
                with self.assertRaisesRegex(ValueError, message):
                    comparisons_file(root, "comparisons:\n  - comparison_id: x\n    project_id: proj\n" + text)


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


DESIGN_HEADER = HEADER.rstrip("\n") + "\tsubject_id\tlibrary_id\n"


def design_manifest(root, rows):
    """rows: (run, replicate, group, subject, library)."""
    path = root / "manifest.tsv"
    lines = [DESIGN_HEADER] + [
        "{r}\t{r}\t{rep}\tproj\tb1\t{g}\tfastq\tr1.fq.gz\t\tSE\tunstranded\tref\ttrue\t{s}\t{l}\n".format(
            r=run, rep=replicate, g=group, s=subject, l=library)
        for run, replicate, group, subject, library in rows]
    path.write_text("".join(lines))
    return load_umbrella_manifest(path)


PAIRED_ROWS = [("a1", "A1", "case", "S1", ""), ("a2", "A2", "case", "S2", ""), ("a3", "A3", "case", "S3", ""),
               ("c1", "C1", "ctrl", "S1", ""), ("c2", "C2", "ctrl", "S2", ""), ("c3", "C3", "ctrl", "S3", "")]


class DesignTests(unittest.TestCase):
    def run_design(self, rows, design, **extra):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            comparison = dict(COMPARISON, **extra)
            if design is not None:
                comparison["design"] = design
            comparison = {key: value for key, value in comparison.items() if value is not None}
            (root / "comparisons.yaml").write_text(json.dumps({"comparisons": [comparison]}))
            loaded = load_comparisons(root / "comparisons.yaml")[0]
            return preflight(design_manifest(root, rows), loaded)

    def test_unpaired_preflight_is_byte_identical(self):
        # sha256 of the preflight JSON written by the code before designs existed
        # (feat/umbrella-circrna 4baed19), for a manifest that has a subject_id
        # column and a merged library: declaring no design changes nothing
        import hashlib
        rows = [("a1", "A1", "case", "S1", ""), ("a1b", "A1", "case", "S1", "LA1"),
                ("a1c", "A1", "case", "S1", "LA1"), ("a2", "A2", "case", "S2", ""),
                ("c1", "C1", "ctrl", "S1", ""), ("c2", "C2", "ctrl", "S2", "")]
        expected = {"x": "81071c8fd5ced9f271ad86ee1403e9ab3c328ff782c86e2afa0f78191c07c313",
                    "y": "5bbce8dec94680d2ef568095057f69e3be58d04076f4da739ff41cfccefd3492"}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest_ = design_manifest(root, rows)
            comparisons = comparisons_file(root, (
                "comparisons:\n  - {comparison_id: x, project_id: proj, group_a: case, group_b: ctrl}\n"
                "  - comparison_id: y\n    project_id: proj\n"
                "    a: {label: s, where: {subject_id: [S1, S2]}, groups: case}\n"
                "    group_b: ctrl\n    exclude: [c2]\n"))
            for comparison in comparisons:
                result = preflight(manifest_, comparison)
                self.assertNotIn("design", result)
                self.assertNotIn("blocks", result)
                text = json.dumps(result, sort_keys=True, indent=2) + "\n"
                self.assertEqual(hashlib.sha256(text.encode()).hexdigest(), expected[comparison["comparison_id"]])

    def test_paired_design_emits_ordered_blocks(self):
        result = self.run_design(PAIRED_ROWS, {"type": "paired", "block": "subject_id"})
        self.assertTrue(result["inference_supported"], result["reasons"])
        self.assertEqual(result["design"]["type"], "paired")
        self.assertEqual(result["design"]["block"], "subject_id")
        self.assertEqual(result["design"]["unbalanced"], "refuse")
        self.assertEqual(result["design"]["residual_df"], 2)
        self.assertEqual(result["design"]["replicate_block"]["C2"], "S2")
        self.assertEqual(result["blocks"], [{"block": "S1", "a": ["A1"], "b": ["C1"]},
                                            {"block": "S2", "a": ["A2"], "b": ["C2"]},
                                            {"block": "S3", "a": ["A3"], "b": ["C3"]}])
        # the comparison and its runs are those of the unpaired comparison
        unpaired = self.run_design(PAIRED_ROWS, None)
        self.assertEqual(result["runs"], unpaired["runs"])
        self.assertEqual(result["collapse"], unpaired["collapse"])

    def test_blocked_design_allows_several_replicates_per_block(self):
        rows = PAIRED_ROWS[:3] + [("a4", "A4", "case", "S1", "")] + PAIRED_ROWS[3:]
        result = self.run_design(rows, {"type": "blocked", "block": "subject_id"})
        self.assertEqual(result["blocks"][0], {"block": "S1", "a": ["A1", "A4"], "b": ["C1"]})
        self.assertEqual(result["design"]["residual_df"], 7 - 3 - 1)
        with self.assertRaisesRegex(ValueError, "exactly one biological replicate per side.*S1"):
            self.run_design(rows, {"type": "paired", "block": "subject_id"})

    def test_one_sided_blocks_are_refused_or_dropped(self):
        rows = PAIRED_ROWS + [("a9", "A9", "case", "S9", "")]
        with self.assertRaisesRegex(ValueError, "one side only \\(S9\\).*unbalanced: drop"):
            self.run_design(rows, {"type": "paired", "block": "subject_id"})
        result = self.run_design(rows, {"type": "paired", "block": "subject_id", "unbalanced": "drop"})
        self.assertEqual(result["design"]["dropped_blocks"], ["S9"])
        self.assertEqual(result["design"]["dropped_replicates"], ["A9"])
        self.assertEqual(result["runs"]["a"], ["a1", "a2", "a3"])
        self.assertNotIn("a9", result["collapse"])

    def test_too_few_blocks_and_no_shared_block_are_refused(self):
        two = [row for row in PAIRED_ROWS if row[3] != "S3"]
        with self.assertRaisesRegex(ValueError, "no subject_id block spans both sides"):
            self.run_design([("a1", "A1", "case", "S1", ""), ("a2", "A2", "case", "S1", ""),
                             ("c1", "C1", "ctrl", "S2", ""), ("c2", "C2", "ctrl", "S2", "")],
                            {"type": "blocked", "block": "subject_id", "unbalanced": "drop"})
        with self.assertRaisesRegex(ValueError, "1 complete subject_id block"):
            self.run_design([row for row in two if row[3] == "S1"] + [("a5", "A5", "case", "S1", "")],
                            {"type": "blocked", "block": "subject_id"})
        # two pairs: ~ block + group leaves 4 - 2 - 1 = 1 residual df, accepted
        self.assertEqual(self.run_design(two, {"type": "paired", "block": "subject_id"})["design"]["residual_df"], 1)

    def test_block_values_must_be_set_and_constant_within_a_replicate(self):
        empty = [("a1", "A1", "case", "", "")] + PAIRED_ROWS[1:]
        with self.assertRaisesRegex(ValueError, "run a1 has no subject_id"):
            self.run_design(empty, {"type": "paired", "block": "subject_id"})
        technical = PAIRED_ROWS + [("a1b", "A1", "case", "S2", "")]
        with self.assertRaisesRegex(ValueError, "replicate A1 has several subject_id values \\(S1, S2\\)"):
            self.run_design(technical, {"type": "paired", "block": "subject_id"})
        merged = PAIRED_ROWS + [("a1b", "A1", "case", "S1", "LA1"), ("a1c", "A1", "case", "S9", "LA1")]
        with self.assertRaisesRegex(ValueError, "library LA1 merges runs with different subject_id values \\(S1;S9\\)"):
            self.run_design(merged, {"type": "paired", "block": "subject_id"})
        with self.assertRaisesRegex(ValueError, "no manifest column donor"):
            self.run_design(PAIRED_ROWS, {"type": "paired", "block": "donor"})

    def test_a_replicate_on_both_sides_is_still_refused(self):
        rows = [("a1", "A1", "case", "S1", ""), ("a2", "A2", "case", "S2", ""),
                ("c1", "A1", "ctrl", "S1", ""), ("c2", "C2", "ctrl", "S2", "")]
        with self.assertRaisesRegex(ValueError, "split across sides: A1"):
            self.run_design(rows, {"type": "paired", "block": "subject_id"},
                            group_a=None, group_b=None,
                            a={"label": "case", "runs": ["a1", "a2"]}, b={"label": "ctrl", "runs": ["c1", "c2"]})

    def test_design_schema_is_validated(self):
        for design, message in (({"type": "matched", "block": "subject_id"}, "paired or blocked"),
                                ({"type": "paired"}, "needs block"),
                                ({"type": "paired", "block": "batch_id"}, "workflow column"),
                                ({"type": "paired", "block": "subject_id", "unbalanced": "keep"}, "refuse or drop"),
                                ({"type": "paired", "block": "subject_id", "covariates": ["sex"]}, "unknown design key"),
                                ("paired", "must be a mapping")):
            with self.assertRaisesRegex(ValueError, message):
                self.run_design(PAIRED_ROWS, design)


class DesignModelTests(unittest.TestCase):
    """src/comparison_design.py: which model each tool records, and the paired statistics."""

    def paired(self, rows=PAIRED_ROWS, design=None):
        with tempfile.TemporaryDirectory() as directory:
            return preflight(design_manifest(Path(directory), rows),
                             dict(COMPARISON, design=design or {"type": "paired", "block": "subject_id"}))

    def test_models_are_recorded_only_for_designed_comparisons(self):
        from src.comparison_design import tool_model
        with tempfile.TemporaryDirectory() as directory:
            unpaired = preflight(design_manifest(Path(directory), PAIRED_ROWS), dict(COMPARISON))
        self.assertEqual(tool_model(unpaired, "deseq2_featurecounts"), {})
        paired = self.paired()
        expected = {"deseq2_featurecounts": "blocked_glm", "deseq2_tximport": "blocked_glm",
                    "circrna_abundance": "blocked_glm", "qapa": "blocked_glm", "dapars2": "paired_test",
                    "suppa2": "paired_native", "rmats": "paired_native", "rmats_conservative": "paired_native",
                    "vasttools": "paired_native", "majiq": "unpaired_adaptation",
                    "whippet_delta": "unpaired_adaptation", "microexonator_delta": "unpaired_adaptation",
                    "leafcutter": "unpaired_adaptation", "psisigma": "unpaired_adaptation"}
        for tool, model in expected.items():
            self.assertEqual(tool_model(paired, tool)["model"], model, tool)
        self.assertEqual(tool_model(paired, "majiq")["model_reason"], "tool has no paired/blocked model")
        self.assertIn("~ subject_id + group", tool_model(paired, "deseq2_tximport")["model_reason"])
        blocked = self.paired(design={"type": "blocked", "block": "subject_id"})
        self.assertEqual(tool_model(blocked, "suppa2")["model"], "unpaired_adaptation")
        self.assertEqual(tool_model(blocked, "deseq2_featurecounts")["model"], "blocked_glm")

    def test_runs_in_block_order_and_technical_runs(self):
        from src.comparison_design import paired_replicates, paired_runs, tool_model
        rows = [("z1", "A1", "case", "S2", ""), ("y2", "A2", "case", "S1", ""),
                ("c1", "C1", "ctrl", "S1", ""), ("c2", "C2", "ctrl", "S2", "")]
        result = self.paired(rows)
        self.assertEqual(paired_replicates(result), (["A2", "A1"], ["C1", "C2"]))
        self.assertEqual(paired_runs(result), (["y2", "z1"], ["c1", "c2"]))
        technical = self.paired(rows + [("z1b", "A1", "case", "S2", "")])
        self.assertIsNone(paired_runs(technical))
        self.assertEqual(tool_model(technical, "rmats")["model"], "unpaired_adaptation")
        self.assertEqual(tool_model(technical, "suppa2")["model"], "paired_native")

    def test_signed_rank_and_block_differences(self):
        from src.comparison_design import block_differences, wilcoxon_signed_rank
        # exact values (scipy.stats.wilcoxon agrees on these untied samples)
        self.assertAlmostEqual(wilcoxon_signed_rank([1, 2, 3, 4, 5, 6]), 0.03125)
        self.assertAlmostEqual(wilcoxon_signed_rank([0.5, -1, 2, 3, -4, 5, 6, 7]), 0.1484375)
        self.assertEqual(wilcoxon_signed_rank([0, 0, 1]), None)
        self.assertEqual(wilcoxon_signed_rank([1, -1]), 1.0)
        # normal approximation beyond the exact limit stays close to the exact value
        values = [i * (-1) ** (i % 3 == 0) for i in range(1, 31)]
        self.assertAlmostEqual(wilcoxon_signed_rank(values), wilcoxon_signed_rank(values, exact_limit=10), 2)
        blocked = {"blocks": [{"block": "S1", "a": ["A1", "A4"], "b": ["C1"]},
                              {"block": "S2", "a": ["A2"], "b": ["C2"]},
                              {"block": "S3", "a": ["A3"], "b": ["C3"]}]}
        values = {"A1": 0.5, "A4": 0.7, "C1": 0.2, "A2": 0.4, "C2": None, "A3": 0.1, "C3": 0.3}
        differences = block_differences(blocked, values)
        self.assertEqual(len(differences), 2)
        self.assertAlmostEqual(differences[0], 0.4)
        self.assertAlmostEqual(differences[1], -0.2)


if __name__ == "__main__":
    unittest.main()
