import os
import re
from types import SimpleNamespace
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_psi_selector(config, paired_samples=()):
    """Load the workflow selector without evaluating the Snakemake workflow."""
    workflow_path = os.path.join(REPOSITORY, "rules", "Whippet_quant.smk")
    with open(workflow_path) as handle:
        source = handle.read()
    match = re.search(
        r"^def get_downstream_PSI\(wildcards\):.*?(?=^rule ME_psi_to_quant:)",
        source,
        flags=re.MULTILINE | re.DOTALL,
    )
    if not match:
        raise AssertionError("get_downstream_PSI is missing from Whippet_quant.smk")
    namespace = {
        "config": config,
        "pe_samples": set(paired_samples),
        "str2bool": lambda value: str(value).lower() in ("true", "t", "1", "yes"),
    }
    exec(compile(match.group(), workflow_path, "exec"), namespace)
    return namespace["get_downstream_PSI"]


class WhippetPSISelectionTests(unittest.TestCase):
    def test_corrected_single_end_psi_is_the_default(self):
        selector = load_psi_selector({})

        self.assertEqual(
            selector(SimpleNamespace(sample="sample_a")),
            "Report/quant/corrected/PSI_sparse/bulk/se/{sample}.corrected.PSI.gz",
        )

    def test_uncorrected_psi_requires_explicit_opt_out(self):
        selector = load_psi_selector({"use_uncorrected_PSI": True})

        self.assertEqual(
            selector(SimpleNamespace(sample="sample_a")),
            "Report/quant/{sample}.out_filtered_ME.PSI.uncorrected.gz",
        )

    def test_corrected_paired_end_psi_remains_paired(self):
        selector = load_psi_selector({}, paired_samples=("sample_a",))

        self.assertEqual(
            selector(SimpleNamespace(sample="sample_a")),
            "Report/quant/corrected/PSI_sparse/bulk/pe/{sample}.corrected.PSI.gz",
        )


if __name__ == "__main__":
    unittest.main()
