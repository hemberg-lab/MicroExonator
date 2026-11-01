import csv
import json
from pathlib import Path
import tempfile
import unittest
import yaml

from src.umbrella_continuation import compose
from src.umbrella_manifest import REQUIRED
from src.umbrella_resume_check import audit
from src.whippet_backfill import archive
from src.umbrella_launch import build_command


class ContinuationTests(unittest.TestCase):
    def fixture(self, root):
        inputs = root / "inputs"
        inputs.mkdir(parents=True)
        columns = list(REQUIRED) + ["subject_id"]
        for name, project, offset in (("old", "P", 1000), ("new", "Q", 2000)):
            rows = []
            for i in range(4):
                run = "SRR" + str(offset + i)
                rows.append([run, run, run, project, "b", "A" if i < 2 else "B", "sra", run, "", "SE",
                             "unstranded", "auto", "true", "donor" + str(i % 2)])
            with (inputs / (name + ".tsv")).open("w") as stream:
                writer = csv.writer(stream, delimiter="\t"); writer.writerow(columns); writer.writerows(rows)
            comparison = dict(comparison_id=name, project_id=project, group_a="A", group_b="B")
            if name == "new":
                comparison["design"] = dict(type="paired", block="subject_id")
            (inputs / (name + ".yaml")).write_text(yaml.safe_dump(dict(comparisons=[comparison])))
        (inputs / "approved.tsv").write_bytes((inputs / "old.tsv").read_bytes())
        (inputs / "approved.yaml").write_bytes((inputs / "old.yaml").read_bytes())
        cfg = dict(working_directory="/previous/run/", umbrella_manifest="/previous/run/inputs/old.tsv",
                   umbrella_comparisons="/previous/run/inputs/old.yaml", max_read_len=150,
                   min_number_files_detected=3)
        (root / "config.yaml").write_text(yaml.safe_dump(cfg))
        return dict(additional_manifest=inputs / "new.tsv", additional_comparisons=inputs / "new.yaml",
                    approved_manifest=inputs / "approved.tsv", approved_comparisons=inputs / "approved.yaml",
                    previous_root="/previous/run", reference_id="ref")

    def test_combines_approved_inputs_keeps_old_identity_and_defers_pairing(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "run"
            args = self.fixture(root)
            before = (root / "inputs/old.tsv").read_bytes()
            result = compose(root, root / "continuation", **args)
            self.assertEqual((result["raw_runs"], result["processing_inputs"], result["comparisons"]), (8, 8, 2))
            self.assertEqual(result["old_processing_identities_unchanged"], 4)
            cfg = yaml.safe_load((root / "continuation/config.continue.yaml").read_text())
            self.assertEqual(cfg["max_read_len"], 150)
            self.assertEqual(cfg["min_number_files_detected"], 3)
            self.assertEqual(cfg["whippet_flags"], "--circ")
            comps = yaml.safe_load((root / "continuation/comparisons.yaml").read_text())["comparisons"]
            self.assertTrue(all("design" not in c for c in comps))
            self.assertEqual((root / "inputs/old.tsv").read_bytes(), before)
            with self.assertRaises(FileExistsError):
                compose(root, root / "continuation", **args)

    def test_refuses_changed_old_source_before_creating_plan(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "run"
            args = self.fixture(root)
            old = root / "inputs/old.tsv"
            text = old.read_text().replace("sra\tSRR1000", "sra\tSRR9999")
            old.write_text(text)
            with self.assertRaisesRegex(ValueError, "metadata differs"):
                compose(root, root / "continuation", **args)
            self.assertFalse((root / "continuation").exists())

    def test_resume_gate_catches_old_producer_and_ignores_new(self):
        log = "rule umbrella_salmon_quant:\n    wildcards: run_id=old, project_id=P\n\nrule correct_quant:\n    wildcards: sample=new\n"
        self.assertEqual(audit(log, {"old"}), [("umbrella_salmon_quant", "old")])

    def test_backfill_archive_preserves_originals_and_refuses_repeat(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); plan = root / "continuation"; plan.mkdir()
            old = root / "umbrella/work/ref/P/b/r/whippet"; old.mkdir(parents=True)
            (old / "quant.psi.gz").write_bytes(b"original")
            (plan / "backfill-archive.json").write_text(json.dumps([str(old.relative_to(root))]))
            archive(root, plan)
            self.assertTrue(old.exists())
            result = archive(root, plan, apply=True)
            self.assertEqual((Path(result["backup"]) / old.relative_to(root) / "quant.psi.gz").read_bytes(), b"original")
            with self.assertRaises(ValueError):
                archive(root, plan, apply=True)

    def test_launcher_defaults_to_dry_run_and_requires_completed_backup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "run"; args = self.fixture(root)
            plan = root / "continuation"; compose(root, plan, **args)
            with self.assertRaises(FileNotFoundError):
                build_command(root, plan, Path("profile"), "backfill")
            (plan / "backfill-receipt.json").write_text(json.dumps({"complete": True}))
            cmd = build_command(root, plan, Path("profile"), "backfill")
            self.assertIn("-n", cmd)
            self.assertEqual(cmd[1:3], ["-s", "MicroExonator.smk"])
            self.assertTrue(cmd[3].startswith("umbrella/work/"))
            self.assertNotIn("quant_umbrella", cmd)

    def test_continuation_launcher_refuses_changed_verified_selection(self):
        import hashlib
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "run"; args = self.fixture(root)
            plan = root / "continuation"; compose(root, plan, **args)
            record = {key: hashlib.sha256((plan / name).read_bytes()).hexdigest() for name, key in
                      (("manifest.tsv", "continuation_manifest_sha256"),
                       ("comparisons.yaml", "continuation_comparisons_sha256"))}
            (plan / "backfill-verified.json").write_text(json.dumps(record))
            cmd = build_command(root, plan, Path("profile"), "continue", submit=True)
            self.assertEqual(cmd[3], "quant_umbrella")
            self.assertNotIn("-n", cmd)
            with (plan / "comparisons.yaml").open("a") as stream:
                stream.write("# changed\n")
            with self.assertRaisesRegex(ValueError, "selection changed"):
                build_command(root, plan, Path("profile"), "continue")
