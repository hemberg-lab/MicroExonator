"""Outputs that can only be remade from the reads are write-protected."""

import os
import re
import unittest


REPOSITORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def rule_outputs(path, rule):
    text = open(os.path.join(REPOSITORY, "rules", path)).read()
    body = re.search(r"^\s*rule {}:(.*?)(?=^\s*rule \w+:|\Z)".format(rule), text, re.S | re.M).group(1)
    return re.search(r"output:(.*?)(?=^\s+(?:params|wildcard_constraints|log|threads|conda|run|shell|resources):)",
                     body, re.S | re.M).group(1)


class ProtectionTests(unittest.TestCase):
    def test_outputs_needing_the_reads_or_hours_to_rebuild_are_protected(self):
        for path, rule, count in (
                ("umbrella_quantifiers.smk", "umbrella_whippet_quant", 5),
                ("umbrella_inputs.smk", "umbrella_validate_reads", 1),
                ("umbrella_reference.smk", "umbrella_reference_identity", 1),
                ("umbrella_reference.smk", "umbrella_hisat2_index", 1),
                ("umbrella_reference.smk", "umbrella_whippet_index", 1),
                ("umbrella_reference.smk", "umbrella_salmon_index", 1),
                ("umbrella_qc.smk", "umbrella_fastqc", 1),
                ("umbrella_qc.smk", "umbrella_qc_keep", 1)):
            outputs = rule_outputs(path, rule)
            self.assertEqual(outputs.count("protected("), count, rule)

    def test_staging_command_text_is_unchanged(self):
        # staged reads left on disk keep Snakemake's record of this command;
        # editing its text would restage (download) every such run
        text = open(os.path.join(REPOSITORY, "rules", "umbrella_inputs.smk")).read()
        self.assertIn('"--threads {threads} {params.tmpdir}"', text)
        self.assertIn('"--run-id {wildcards.run_id:q} --r1 {output:q} --threads {threads} {params.tmpdir}"', text)

    def test_comparison_outputs_stay_recomputable(self):
        text = open(os.path.join(REPOSITORY, "rules", "umbrella_comparisons.smk")).read()
        self.assertNotIn("protected(", text)


if __name__ == "__main__":
    unittest.main()
