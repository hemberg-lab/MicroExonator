"""Lifecycle tests: losing guards must permit an observable unsafe transition."""
import importlib.util
import json
import tempfile
import unittest
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ExecutionLedgerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def api(self):
        spec = importlib.util.find_spec("src.umbrella_execution")
        self.assertIsNotNone(spec, "guarded packet execution controller is missing")
        from src.umbrella_execution import PacketLedger
        return PacketLedger

    def test_complement_keeps_discovery_hashes_and_history(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            discovery = ledger.begin("scout", "discovery", ["L1", "L2"], ["C1"],
                                     ["microexonator", "whippet"])
            me = self.root / "me.tsv"; me.write_text("event\tPSI\nE1\t0.4\n")
            wh = self.root / "wh.tsv"; wh.write_text("event\tPSI\nE1\t0.5\n")
            ledger.complete(discovery, {"microexonator": [me], "whippet": [wh]})
            before = me.read_bytes(), wh.read_bytes()
            complement = ledger.begin("confirm", "complement", ["L1", "L2"], ["C1"], ["salmon"])
            self.assertEqual(complement["completed_tools"], [])
            salmon = self.root / "salmon.tsv"; salmon.write_text("gene\tcount\nG1\t9\n")
            ledger.complete(complement, {"salmon": [salmon]})
            self.assertEqual((me.read_bytes(), wh.read_bytes()), before)
            self.assertEqual(ledger.completed_tools(), ["microexonator", "salmon", "whippet"])
            self.assertEqual(len(ledger.state["executions"]), 2)
            self.assertEqual(ledger.state["executions"]["scout"]["status"], "complete")

    def test_tampered_or_missing_discovery_rejected_before_complement(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            output = self.root / "psi.tsv"; output.write_text("original\n")
            ledger.complete(execution, {"whippet": [output]})
            output.chmod(0o644); output.write_text("changed\n")
            with self.assertRaisesRegex(ValueError, "retained"):
                ledger.begin("confirm", "complement", ["L1"], [], ["salmon"])
            output.unlink()
            with self.assertRaisesRegex(ValueError, "retained"):
                ledger.verify()

    def test_identity_change_and_unrequested_completion_rejected(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"genome": "v1"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            with self.assertRaisesRegex(ValueError, "requested"):
                ledger.complete(execution, {"salmon": []})
            ledger.fail(execution, "interrupted")
        with self.assertRaisesRegex(ValueError, "identity"):
            with Ledger(self.root / "packet", {"genome": "v2"}):
                pass

    def test_resume_preserves_attempts_and_refuses_selection_change(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            ledger.fail(execution, "interrupted")
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            self.assertEqual(len(execution["attempts"]), 2)
            self.assertEqual(execution["attempts"][0]["status"], "failed")
            with self.assertRaisesRegex(ValueError, "selection"):
                ledger.begin("scout", "discovery", ["L2"], [], ["whippet"])
            self.assertEqual(ledger.completed_tools(), [])

    def test_complement_requires_discovery_and_disjoint_tools(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            with self.assertRaisesRegex(ValueError, "discovery"):
                ledger.begin("confirm", "complement", ["L1"], [], ["salmon"])
            with self.assertRaisesRegex(ValueError, "profile"):
                ledger.begin("scout", "discovery", ["L1"], [], ["salmon"])

    def test_missing_partial_inventory_is_not_marked_complete(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["microexonator", "whippet"])
            with self.assertRaisesRegex(ValueError, "inventory"):
                ledger.complete(execution, {"whippet": []})
            self.assertEqual(ledger.completed_tools(), [])

    def test_packet_lock_excludes_concurrent_execution(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}):
            with self.assertRaisesRegex(ValueError, "locked"):
                with Ledger(self.root / "packet", {"manifest": "frozen"}):
                    pass

    def test_completed_discovery_cannot_be_repeated_under_new_id(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            output = self.root / "psi.tsv"; output.write_text("original\n")
            ledger.complete(execution, {"whippet": [output]})
            with self.assertRaisesRegex(ValueError, "completed discovery"):
                ledger.begin("scout2", "discovery", ["L1"], [], ["whippet"])

    def test_controller_selection_does_not_edit_master_and_keeps_lane_closure(self):
        self.api()
        from src import umbrella_execution as execution
        self.assertTrue(hasattr(execution, "execution_view"), "execution selection adapter is missing")
        manifest = self.root / "manifest.tsv"
        manifest.write_text("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\tlibrary_id\n"
                            "S1\tSRR1\tB1\tP\tbatch1\tA\tsra\tSRR1\t\tSE\tunstranded\tref\ttrue\tL1\n"
                            "S1\tSRR2\tB1\tP\tbatch1\tA\tsra\tSRR2\t\tSE\tunstranded\tref\ttrue\tL1\n"
                            "S2\tSRR3\tB2\tP\tbatch1\tB\tsra\tSRR3\t\tSE\tunstranded\tref\ttrue\tL2\n")
        comparisons = self.root / "comparisons.yaml"
        comparisons.write_text("comparisons:\n- comparison_id: C1\n  project_id: P\n  group_a: A\n  group_b: B\n")
        before = manifest.read_bytes()
        view = execution.execution_view(manifest, comparisons, ["L1", "SRR3"], ["C1"])
        self.assertEqual([row["run_id"] for row in view["rows"]], ["SRR1", "SRR2", "SRR3"])
        self.assertEqual(view["libraries"], ["L1", "L2"])
        self.assertEqual(manifest.read_bytes(), before)
        with self.assertRaisesRegex(ValueError, "unknown"):
            execution.execution_view(manifest, comparisons, ["typo"], ["C1"])
        with self.assertRaisesRegex(ValueError, "empty"):
            execution.execution_view(manifest, comparisons, [], ["C1"])

    def test_complement_dag_guard_refuses_discovery_producers(self):
        self.api()
        from src import umbrella_execution as execution
        self.assertTrue(hasattr(execution, "validate_dry_run"), "DAG guard is missing")
        with self.assertRaisesRegex(ValueError, "discovery producer"):
            execution.validate_dry_run("rule umbrella_whippet_quant:\n", "complement")
        execution.validate_dry_run("rule umbrella_salmon_quant:\nrule umbrella_hisat2:\n", "complement")

    def test_comparison_only_selection_limits_quantification(self):
        from src.umbrella_execution import execution_view
        manifest = self.root / "manifest.tsv"
        manifest.write_text("sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\n" +
                            "".join("S{0}\tSRR{0}\tB{0}\tP\tbatch1\t{1}\tsra\tSRR{0}\t\tSE\tunstranded\tref\ttrue\n".format(i, "A" if i % 2 else "B") for i in range(1, 5)))
        comparisons = self.root / "comparisons.yaml"
        comparisons.write_text("comparisons:\n- comparison_id: C1\n  project_id: P\n  a:\n    label: A\n    runs: [SRR1]\n  b:\n    label: B\n    runs: [SRR2]\n")
        view = execution_view(manifest, comparisons, comparisons=["C1"])
        self.assertEqual(view["libraries"], ["SRR1", "SRR2"])

    def test_invalid_gzip_cannot_be_certified_complete(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            product = self.root / "psi.tsv.gz"; product.write_bytes(b"not a gzip file")
            with self.assertRaisesRegex(ValueError, "invalid.*product"):
                ledger.complete(execution, {"whippet": [product]})
            self.assertEqual(ledger.completed_tools(), [])

    def test_execution_reports_are_hashed_and_collected_cumulatively(self):
        Ledger = self.api()
        with Ledger(self.root / "packet", {"manifest": "frozen"}) as ledger:
            execution = ledger.begin("scout", "discovery", ["L1"], [], ["whippet"])
            product = self.root / "psi.tsv"; product.write_text("E1\t0.4\n")
            report = self.root / "synthesis.tsv"; report.write_text("tool\tstatus\nwhippet\tcomplete\n")
            ledger.complete(execution, {"whippet": [product]}, artifacts={"synthesis": [report]})
            self.assertTrue((ledger.root / "cumulative_reports.json").exists())
            report.chmod(0o644); report.write_text("changed report\n")
            with self.assertRaisesRegex(ValueError, "retained"):
                ledger.verify()

    def test_index_companions_and_seed_assets_participate_in_identity(self):
        self.api()
        from src.umbrella_execution import packet_identity
        (self.root / "Snakefile").write_text("rule all:\n    input: []\n")
        (self.root / "whippet.jls").write_bytes(b"index")
        companion = self.root / "whippet.jls.exons.tab.gz"; companion.write_bytes(b"exon identities")
        (self.root / "data").mkdir()
        tag = self.root / "data/ME_canonical_SJ_tags.DB.fa"; tag.write_text(">ME1\nACGT\n")
        config = {"umbrella_reference": {"whippet_index": str(self.root / "whippet.jls")}}
        before = packet_identity(config, self.root, "Snakefile")
        companion.write_bytes(b"changed exon identities")
        self.assertNotEqual(packet_identity(config, self.root, "Snakefile"), before)
        companion.write_bytes(b"exon identities")
        tag.write_text(">ME1\nAAAA\n")
        self.assertNotEqual(packet_identity(config, self.root, "Snakefile"), before)

    def test_attached_short_options_cannot_redirect_checked_workspace(self):
        self.api()
        from src import umbrella_execution as execution
        self.assertTrue(hasattr(execution, "validate_arguments"), "bounded launcher options are missing")
        for arguments in (["-d/tmp/escape"], ["-sOtherfile"], ["--quiet"], ["other_target"], ["--forceall"], ["--dry-run"]):
            with self.subTest(arguments=arguments), self.assertRaises(ValueError):
                execution.validate_arguments(arguments)
        execution.validate_arguments(["--profile", "base", "--use-conda", "-j", "50", "--resources", "get_data=8", "-k"])

    def test_profile_cannot_override_phase_workspace_or_turn_run_into_inspection(self):
        self.api()
        from src import umbrella_execution as execution
        self.assertTrue(hasattr(execution, "validated_profile"), "profile validation is missing")
        profile = self.root / "base"; profile.mkdir()
        config = profile / "config.yaml"
        for text in ("directory: /tmp/escape\n", "dry-run: true\n", "config: [umbrella_optional.whippet=true]\n"):
            config.write_text(text)
            with self.assertRaisesRegex(ValueError, "profile"):
                execution.validated_profile("base", self.root)
        config.write_text("jobs: 50\ncluster: qsub\ncluster-config: cluster.PBS.json\n")
        self.assertEqual(execution.validated_profile("base", self.root), profile.resolve())


class TwoPassTransitionTests(unittest.TestCase):
    """Real Snakemake scheduling/cleanup; toy products stand in for costly tools."""
    def test_direct_launch_cannot_silently_ignore_enabled_execution_mode(self):
        from tests.test_workflow_dag import WorkflowSelectionTests
        result = WorkflowSelectionTests().run_quant_dry_run(
            umbrella=True, target="quant_umbrella", extra_config=(
                "umbrella_execution:", "  enabled: true", "  packet_id: P", "  execution_id: scout", "  mode: discovery"))
        self.assertNotEqual(result.returncode, 0, "direct launch silently ignored execution profile")
        self.assertIn("umbrella_execution.py", result.stdout)

    def test_native_workflow_adapter_discovery_dry_run(self):
        from src.umbrella_execution import run_execution
        from tests.test_workflow_dag import WorkflowSelectionTests
        import yaml
        results = []
        def prepare(folder):
            root = Path(folder)
            config_path = root / "config.yaml"
            config = yaml.safe_load(config_path.read_text())
            config["umbrella_execution"] = {"enabled": True, "packet_id": "native",
                                            "execution_id": "scout", "mode": "discovery"}
            config_path.write_text(yaml.safe_dump(config))
            command = [sys.executable, "-c", "import appdirs, runpy; appdirs.system='linux'; runpy.run_module('snakemake', run_name='__main__')"]
            results.append(run_execution(config_path, root, command, ["-j", "1"], dry_run=True))
        fixture = WorkflowSelectionTests()
        result = fixture.run_quant_dry_run(umbrella=True, target="quant_umbrella", prepare=prepare,
                                          extra_config=("skip_discovery: T", "only_db: T", "delta_method: whippet"),
                                          umbrella_comparisons="comparisons:\n- comparison_id: C1\n  project_id: project\n  group_a: control\n  group_b: case\n",
                                          umbrella_extra_rows=["sample_b\trun_b\trep_b\tproject\tbatch\tcase\tfastq\tr1.fastq.gz\tr2.fastq.gz\tPE\tunstranded\tref\ttrue"])
        self.assertNotEqual(result.returncode, 0, "direct launch must refuse the master execution config")
        self.assertIn("umbrella_execution.py", result.stdout)
        self.assertIn("umbrella_whippet_quant", results[0]["scheduled_rules"])
        self.assertIn("umbrella_microexonator_delta", results[0]["scheduled_rules"])
        self.assertIn("umbrella_whippet_delta", results[0]["scheduled_rules"])
        self.assertNotIn("umbrella_microexonator_whippet_delta", results[0]["scheduled_rules"])
        self.assertEqual(results[0]["delta_method"], "microexonator")
        self.assertNotIn("umbrella_hisat2", results[0]["scheduled_rules"])

    def test_completed_discovery_cleanup_then_selected_complement(self):
        from src import umbrella_execution as execution
        self.assertTrue(hasattr(execution, "run_execution"), "guarded launcher is missing")
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "genome.fa").write_text(">chr1\nACGT\n")
            (root / "index.bin").write_bytes(b"immutable index")
            (root / "index.bin.exons.tab.gz").write_bytes(b"synthetic exon sidecar")
            header = "sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\n"
            (root / "manifest.tsv").write_text(header + "".join(
                "S{0}\tSRR{0}\tB{0}\tP\tbatch1\t{1}\tsra\tSRR{0}\t\tSE\tunstranded\tref\ttrue\n".format(i, group)
                for i, group in ((1, "A"), (2, "B"), (3, "A"), (4, "B"))))
            (root / "comparisons.json").write_text(json.dumps({"comparisons": [
                {"comparison_id": "C1", "project_id": "P", "group_a": "A", "group_b": "B"}]}))
            (root / "Snakefile").write_text('''import csv
configfile: "config.yaml"
with open(config["umbrella_manifest"]) as stream:
    libraries = [row["run_id"] for row in csv.DictReader(stream, delimiter="\\t")]
rule quant_umbrella:
    input:
        (["comparisons/microexonator_delta.tsv", "comparisons/whippet_delta.diff.gz"]
         if config["umbrella_optional"]["microexonator"] else ["genes/salmon_counts.tsv.gz"])
rule umbrella_stage_reads:
    output: temp("reads.fastq")
    shell: "printf 'tiny reads' > {output}"
rule umbrella_microexonator_delta:
    input: "reads.fastq"
    output: "comparisons/microexonator_delta.tsv"
    run:
        with open(output[0], "w") as stream:
            stream.write("event\\tPSI\\nE1\\t0.4\\n")
rule umbrella_whippet_delta:
    input: "reads.fastq"
    output: "comparisons/whippet_delta.diff.gz"
    run:
        import gzip
        with gzip.open(output[0], "wt") as stream:
            stream.write("event\\tPSI\\nE1\\t0.5\\n")
rule umbrella_salmon_shard:
    input: "reads.fastq"
    output: "genes/salmon_counts.tsv.gz"
    run:
        import gzip
        with gzip.open(output[0], "wt") as stream:
            stream.write("gene\\t" + "\\t".join(libraries) + "\\nG1\\t" + "\\t".join("9" for _ in libraries) + "\\n")
''')
            config = {"umbrella_manifest": "manifest.tsv", "umbrella_comparisons": "comparisons.json",
                      "umbrella_reference": {"genome_fasta": "genome.fa", "whippet_index": "index.bin"},
                      "umbrella_execution": {"enabled": True, "packet_id": "tiny", "execution_id": "scout",
                                             "mode": "discovery"}}
            config_path = root / "master.json"
            config_path.write_text(json.dumps(config))
            command = [sys.executable, "-c", "import appdirs, runpy; appdirs.system='linux'; runpy.run_module('snakemake', run_name='__main__')"]
            (root / "base").mkdir()
            (root / "base/config.yaml").write_text("jobs: 1\n")
            arguments = ["--profile", "base", "-j", "1", "--force-use-threads"]  # avoid macOS child-process source cache outside sandbox
            discovery = execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile")
            self.assertEqual(discovery["status"], "complete")
            work = root / "umbrella/executions/tiny/scout/work"
            self.assertFalse((work / "reads.fastq").exists(), "normal temp cleanup did not happen")
            retained = {path: Path(path).read_bytes() for files in discovery["outputs"].values() for path in files}
            master_before = (root / "manifest.tsv").read_bytes()
            config["umbrella_execution"].update(execution_id="confirm", mode="complement", tools=["salmon"],
                                                libraries=["SRR1", "SRR2"])
            config_path.write_text(json.dumps(config))
            plan = execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile", dry_run=True)
            self.assertNotIn("umbrella_microexonator_delta", plan["scheduled_rules"])
            self.assertNotIn("umbrella_whippet_delta", plan["scheduled_rules"])
            complement = execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile")
            self.assertEqual(complement["completed_tools"], ["salmon"])
            import gzip
            product = next(iter(complement["outputs"]["salmon"]))
            with gzip.open(product, "rt") as stream:
                self.assertEqual(stream.readline().strip(), "gene\tSRR1\tSRR2")
            self.assertEqual((root / "manifest.tsv").read_bytes(), master_before)
            for path, before in retained.items():
                self.assertEqual(Path(path).read_bytes(), before)
            packet = root / "umbrella/executions/tiny"
            cumulative = json.loads((packet / "cumulative_inventory.json").read_text())
            self.assertEqual({row["tool"] for row in cumulative}, {"microexonator", "whippet", "salmon"})
            config["umbrella_execution"]["execution_id"] = "bad"
            (root / "genome.fa").write_text(">chr1\nAAAA\n")
            config_path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "identity"):
                execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile", dry_run=True)
            (root / "genome.fa").write_text(">chr1\nACGT\n")
            (root / "index.bin").write_bytes(b"changed index")
            with self.assertRaisesRegex(ValueError, "identity"):
                execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile", dry_run=True)
            (root / "index.bin").write_bytes(b"immutable index")
            (root / "manifest.tsv").write_bytes(master_before.replace(b"\tA\t", b"\tNEW_GROUP\t"))
            with self.assertRaises(ValueError):
                execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile", dry_run=True)
            (root / "manifest.tsv").write_bytes(master_before)
            config["whippet_flags"] = "--different-producer-option"
            config_path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "identity"):
                execution.run_execution(config_path, root, command, arguments, snakefile="Snakefile", dry_run=True)
            config.pop("whippet_flags")
            config_path.write_text(json.dumps(config))
            with self.assertRaisesRegex(ValueError, "force"):
                execution.run_execution(config_path, root, command, arguments + ["--forcerun", "umbrella_whippet_delta"],
                                        snakefile="Snakefile", dry_run=True)


if __name__ == "__main__":
    unittest.main()
