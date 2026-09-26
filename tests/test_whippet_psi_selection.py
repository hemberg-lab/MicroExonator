import os
import re
from types import SimpleNamespace
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_psi_selector(config, paired_samples=()):
    """Load the PSI table selector (shared by Whippet and the Whippet-free delta) without evaluating the workflow."""
    workflow_path = os.path.join(REPOSITORY, "MicroExonator.smk")
    with open(workflow_path) as handle:
        source = handle.read()
    match = re.search(
        r"^def downstream_PSI\(sample\):.*?(?=^def partition)",
        source,
        flags=re.MULTILINE | re.DOTALL,
    )
    if not match:
        raise AssertionError("downstream_PSI is missing from MicroExonator.smk")
    namespace = {
        "config": config,
        "pe_samples": set(paired_samples),
        "str2bool": lambda value: str(value).lower() in ("true", "t", "1", "yes"),
    }
    exec(compile(match.group(), workflow_path, "exec"), namespace)
    return lambda wildcards: namespace["downstream_PSI"](wildcards.sample)


class WhippetPSISelectionTests(unittest.TestCase):
    def test_corrected_single_end_psi_is_the_default(self):
        selector = load_psi_selector({})

        self.assertEqual(
            selector(SimpleNamespace(sample="sample_a")),
            "Report/quant/corrected/PSI_sparse/bulk/se/sample_a.corrected.PSI.gz",
        )

    def test_uncorrected_psi_requires_explicit_opt_out(self):
        selector = load_psi_selector({"use_uncorrected_PSI": True})

        self.assertEqual(
            selector(SimpleNamespace(sample="sample_a")),
            "Report/quant/sample_a.out_filtered_ME.PSI.uncorrected.gz",
        )

    def test_corrected_paired_end_psi_remains_paired(self):
        selector = load_psi_selector({}, paired_samples=("sample_a",))

        self.assertEqual(
            selector(SimpleNamespace(sample="sample_a")),
            "Report/quant/corrected/PSI_sparse/bulk/pe/sample_a.corrected.PSI.gz",
        )


if __name__ == "__main__":
    unittest.main()
