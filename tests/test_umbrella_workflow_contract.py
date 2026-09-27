"""Static DAG interface checks when Snakemake is not installed locally."""

import pathlib
import re
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]


class UmbrellaWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.main = (ROOT / "MicroExonator.smk").read_text()
        cls.data = (ROOT / "rules/Get_data.smk").read_text()
        cls.intake = (ROOT / "rules/umbrella_inputs.smk").read_text()

    def test_shared_reference_rules_and_manifest_runs_are_available(self):
        self.assertRegex(self.main, r'DATA\.update\([^\n]*UMBRELLA_MANIFEST\.included_runs\(')
        self.assertRegex(self.main, r'(?m)^include\s*:\s*"rules/Get_data.smk"')
        self.assertIn('if "umbrella_manifest" in config:', self.data)
        self.assertIn('rule generate_fasta_from_bed12:', self.data)

    def test_legacy_fastq_has_a_manifest_backed_bridge(self):
        self.assertIn('rule umbrella_legacy_bridge:', self.intake)
        self.assertIn('UMBRELLA_MANIFEST.legacy_fastq', self.intake)
        self.assertIn('else temp("FASTQ/{sample}.fastq.gz")', self.intake)
        self.assertIn('rule umbrella_legacy_fastq:', self.intake)

    def test_native_reads_are_temporary_and_manifest_invalidates_outputs(self):
        # Staged reads are removed once all their consumers have run; the
        # legacy FASTQ is relinked (source symlink or hard link) so it
        # outlives them.
        self.assertRegex(self.intake, r'temp\("umbrella/work/\{reference_id\}/\{project_id\}/\{batch_id\}/\{run_id\}/R1\.fastq\.gz"\)')
        self.assertIn('link_legacy_fastq(input.fastq, output[0])', self.intake)
        self.assertIn('manifest=str(UMBRELLA_MANIFEST.path)', self.intake)
        self.assertIn('reads=umbrella_validation_inputs', self.intake)

if __name__ == "__main__":
    unittest.main()
