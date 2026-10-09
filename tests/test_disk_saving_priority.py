"""Scheduling-only priorities must drain readers without invalidating results."""
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]


class DiskSavingPriorityTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue((ROOT / 'src/disk_saving_priority.py').exists(),
                        'disk-saving scheduling policy has not been implemented')

    def command(self, root, *options):
        return subprocess.run([sys.executable, '-c',
            "import appdirs,runpy; appdirs.system='linux'; runpy.run_module('snakemake',run_name='__main__')",
            '-s', 'Snakefile', '-j', '1', '--scheduler', 'greedy', '--nocolor',
            '--force-use-threads', *options], cwd=root, text=True,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            env=dict(os.environ, XDG_CACHE_HOME=str(root / '.cache')))

    def fixture(self, root, enabled=True, fail=False):
        (root / 'src').symlink_to(ROOT / 'src', target_is_directory=True)
        (root / 'a.fq').write_text('ACGT\n')
        (root / 'enabled').write_text('true' if enabled else 'false')
        if fail:
            (root / 'fail').touch()
        (root / 'Snakefile').write_text('''
from pathlib import Path
from src.disk_saving_priority import apply_disk_saving_priorities
config = {'umbrella_manifest': 'unused', 'umbrella_disk_saving_priority': Path('enabled').read_text()}
rule all:
    input: 'a.hisat', 'a.whippet', 'b.done'
rule umbrella_stage_reads:
    output: temp('b.fq')
    shell: "echo stage >> order; echo ACGT > {output}"
rule umbrella_hisat2:
    input: 'a.fq'
    output: 'a.hisat'
    threads: 2
    shell: "echo hisat >> order; cp {input} {output}"
rule umbrella_whippet_quant:
    input: 'a.fq'
    output: 'a.whippet'
    shell: "echo whippet >> order; test ! -e fail; cp {input} {output}"
rule b_consumer:
    input: 'b.fq'
    output: 'b.done'
    shell: "cp {input} {output}"
apply_disk_saving_priorities(workflow, config)
''')
        # a.fq represents already-acquired temporary input.
        text = (root / 'Snakefile').read_text().replace(
            "rule umbrella_stage_reads:", "rule existing_a:\n    output: temp('a.fq')\n    shell: 'echo ACGT > {output}'\nrule umbrella_stage_reads:")
        (root / 'Snakefile').write_text(text)

    def test_ready_fastq_readers_run_before_new_acquisition_and_cleanup(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root)
            result = self.command(root)
            self.assertEqual(result.returncode, 0, result.stdout)
            self.assertEqual((root / 'order').read_text().splitlines()[:3],
                             ['hisat', 'whippet', 'stage'])
            self.assertFalse((root / 'a.fq').exists())
            self.assertFalse((root / 'b.fq').exists())

    def test_failed_consumer_keeps_reads_and_toggle_does_not_rerun_completed_job(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            self.fixture(root, fail=True)
            failed = self.command(root)
            self.assertNotEqual(failed.returncode, 0)
            self.assertTrue((root / 'a.fq').exists())
            before = (root / 'a.hisat').stat().st_mtime_ns
            (root / 'fail').unlink()
            (root / 'enabled').write_text('false')
            resumed = self.command(root, '--rerun-incomplete')
            self.assertEqual(resumed.returncode, 0, resumed.stdout)
            self.assertEqual((root / 'a.hisat').stat().st_mtime_ns, before)
            self.assertFalse((root / 'a.fq').exists())
            (root / 'enabled').write_text('true')
            complete = self.command(root, '-n')
            self.assertEqual(complete.returncode, 0, complete.stdout)
            self.assertIn('Nothing to be done', complete.stdout)

    def test_disabled_preserves_legacy_priorities_and_invalid_setting_fails(self):
        from snakemake.rules import Rule
        from src.disk_saving_priority import apply_disk_saving_priorities
        rule = Rule('Round2_alingment_pre_processing', None)
        workflow = SimpleNamespace(rules=[rule])
        rule.priority = 500
        apply_disk_saving_priorities(workflow, {'umbrella_manifest': 'x'})
        self.assertEqual(rule.priority, 500)
        apply_disk_saving_priorities(workflow, {'umbrella_disk_saving_priority': True})
        self.assertEqual(rule.priority, 500)  # legacy non-umbrella route unaffected
        with self.assertRaises(ValueError):
            apply_disk_saving_priorities(workflow, {'umbrella_manifest': 'x',
                                                   'umbrella_disk_saving_priority': 'maybe'})

    def test_native_core_and_full_dags_keep_the_same_jobs(self):
        from tests.test_workflow_dag import WorkflowSelectionTests
        helper = WorkflowSelectionTests()
        core = ['umbrella_optional: {salmon: false, hisat2: false, rmats: false, coverage: false, leafcutter: false, suppa2: false, multiqc: false}']
        for options in (core, ['umbrella_optional: {leafcutter: false}']):
            off = helper.run_quant_dry_run(umbrella=True, target='quant_umbrella', extra_config=options)
            on = helper.run_quant_dry_run(umbrella=True, target='quant_umbrella',
                extra_config=options + ['umbrella_disk_saving_priority: true'])
            self.assertEqual(off.returncode, 0, off.stdout)
            self.assertEqual(on.returncode, 0, on.stdout)
            self.assertEqual(sorted(re.findall(r'^rule (\w+):', off.stdout, re.M)),
                             sorted(re.findall(r'^rule (\w+):', on.stdout, re.M)))
            block = re.search(r'rule umbrella_whippet_quant:.*?(?=\n\[|\nrule |\Z)', on.stdout, re.S)
            self.assertIsNotNone(block)
            self.assertRegex(block.group(), r'priority: 90')


if __name__ == '__main__':
    unittest.main()
