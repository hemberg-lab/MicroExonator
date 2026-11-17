"""Paired statistics use the dedicated native comparison rules, never new alignments."""
import re
import unittest
from tests import test_workflow_dag as dag


class PairedDagTests(unittest.TestCase):
    def check_route(self, target, expected, mixed=False):
        rows = []
        for name, group, donor in [('b','control','d2'), ('e','control','d1'),
                                    ('c','case','d1'), ('d','case','d2')]:
            rows.append('sample_{0}\trun_{0}\trep_{0}\tproject\tbatch\t{1}\tfastq\tr1.fastq.gz\tr2.fastq.gz\tPE\tunstranded\tref\ttrue\t{2}'.format(name,group,donor))
        comparisons = ('comparisons:\n  - comparison_id: paired\n    project_id: project\n'
                       '    a: {label: case, samples: [sample_c, sample_d]}\n'
                       '    b: {label: control, samples: [sample_e, sample_b]}\n'
                       '    design: {type: paired, block: donor}\n')
        if mixed:
            comparisons += ('  - comparison_id: paired_reference\n    project_id: project\n'
                            '    a: {label: case, samples: [sample_c, sample_d]}\n'
                            '    b: {label: control, samples: [sample_e, sample_a]}\n')
        result = dag.WorkflowSelectionTests.run_quant_dry_run(self, umbrella=True, target=target,
            print_shell=True, umbrella_extra_column='donor', umbrella_extra_rows=rows,
            umbrella_comparisons=comparisons,
            extra_config=['umbrella_modules: [rmats_conservative]', 'umbrella_modules_mode: ingest'])
        self.assertEqual(result.returncode, 0, result.stdout[-5000:])
        self.assertIn('rule ' + expected + ':', result.stdout)
        return result.stdout

    def test_conservative_paired_native_rule_is_selected(self):
        self.check_route('quant_rmats_conservative', 'umbrella_rmats_conservative_compare_paired')

    def test_core_paired_native_rule_is_selected(self):
        self.check_route('quant_umbrella', 'umbrella_rmats_post_paired')

    def test_one_umbrella_run_routes_each_comparison_independently(self):
        output = self.check_route('quant_umbrella', 'umbrella_rmats_post_paired', mixed=True)
        jobs = re.findall(r'^rule (\w+):\n(.*?)(?=\n\n|\Z)', output, re.M | re.S)
        for rule in ('umbrella_rmats_post', 'umbrella_rmats_conservative_compare', 'umbrella_suppa'):
            for suffix, comparison in (('_paired', 'paired'), ('', 'paired_reference')):
                matching = [body for name, body in jobs if name == rule + suffix]
                self.assertEqual(len(matching), 1, rule + suffix)
                self.assertRegex(matching[0], rf'\bcomparison_id={comparison}(?:,|\n|$)')
        self.assertEqual(sum(name == 'umbrella_hisat2' for name, _ in jobs), 5)

    def test_module_without_paired_model_receives_explicit_label(self):
        output = self.check_route('quant_psisigma', 'umbrella_psisigma_compare')
        self.assertIn('"model": "unpaired_adaptation"', output)
