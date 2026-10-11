"""Cached PCA must never acquire reads; plots retain biological group labels."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
import shutil
RSCRIPT = os.environ.get('MICROEXONATOR_RSCRIPT', shutil.which('Rscript') or '')


class ExpressionPCATests(unittest.TestCase):
    def test_standalone_dag_has_only_pca_jobs_and_missing_cache_fails(self):
        try:
            import snakemake
        except ImportError:
            self.skipTest('Snakemake 7 required for DAG test')
        import sys
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'src').symlink_to(ROOT / 'src', target_is_directory=True)
            (root / 'ExpressionPCA.smk').write_text((ROOT / 'ExpressionPCA.smk').read_text())
            comparison = root / 'comparisons/hg38/project/case-control'
            (comparison / 'joined').mkdir(parents=True)
            (comparison / 'preflight.json').write_text('{}')
            counts = comparison / 'joined/featurecounts.tsv.gz'
            counts.write_bytes(b'cached')
            command = [sys.executable, '-c',
                "import appdirs,runpy; appdirs.system='linux'; runpy.run_module('snakemake',run_name='__main__')",
                '-s', 'ExpressionPCA.smk', '-n', '--cores', '1', '--nocolor']
            result = subprocess.run(command, cwd=root, capture_output=True, text=True,
                                    env=dict(os.environ, XDG_CACHE_HOME=str(root / '.cache')))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('cached_expression_pca', result.stdout + result.stderr)
            self.assertNotIn('umbrella_stage_reads', result.stdout + result.stderr)
            self.assertNotIn('umbrella_hisat2', result.stdout + result.stderr)
            counts.unlink()
            failed = subprocess.run(command, cwd=root, capture_output=True, text=True,
                                    env=dict(os.environ, XDG_CACHE_HOME=str(root / '.cache')))
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn('Cached-only PCA', failed.stdout + failed.stderr)

    def test_cached_plan_requires_selected_matrices(self):
        self.assertTrue((ROOT / 'src/expression_pca_cache.py').exists(), 'cached PCA planner missing')
        from src.expression_pca_cache import cached_comparisons
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            comparison = root / 'hg38/project/case-control'
            (comparison / 'joined').mkdir(parents=True)
            (comparison / 'preflight.json').write_text('{}')
            with self.assertRaisesRegex(ValueError, 'featurecounts'):
                cached_comparisons(root, ['featurecounts'])
            (comparison / 'joined/featurecounts.tsv.gz').write_bytes(b'cached')
            jobs = cached_comparisons(root, ['featurecounts'])
            self.assertEqual(list(jobs), ['hg38/project/case-control/featurecounts'])
            with self.assertRaisesRegex(ValueError, 'salmon'):
                cached_comparisons(root, ['tximport'])

    def test_r_plot_small_groups_and_sample_mismatch(self):
        self.assertTrue((ROOT / 'src/expression_pca.R').exists(), 'R PCA implementation missing')
        if not RSCRIPT or not Path(RSCRIPT).is_file():
            self.skipTest('set MICROEXONATOR_RSCRIPT for R plot integration')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            samples = ['a1', 'a2', 'b1', 'b2']
            preflight = root / 'preflight.json'
            preflight.write_text(json.dumps({'group_a': 'Treatment', 'group_b': 'Control',
                'replicates': {'a': samples[:2], 'b': samples[2:]}}))
            counts = root / 'counts.tsv.gz'
            with gzip.open(counts, 'wt') as stream:
                stream.write('gene_id\tlength\t' + '\t'.join(samples) + '\n')
                for i in range(100):
                    stream.write('g%d\t1000\t%d\t%d\t%d\t%d\n' %
                                 (i, 10+i*2, 12+i*3, 18+i, 20+i*4))
            command = [RSCRIPT, str(ROOT / 'src/expression_pca.R'), '--route', 'featurecounts',
                       '--preflight', str(preflight), '--counts', str(counts), '--out', str(root / 'pca')]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((root / 'pca.pdf').read_bytes().startswith(b'%PDF'))
            with gzip.open(root / 'pca.scores.tsv.gz', 'rt') as stream:
                import csv
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual([r['sample_id'] for r in rows], samples)
            self.assertEqual([r['group'] for r in rows], ['Treatment']*2 + ['Control']*2)
            self.assertEqual({r['qc_status'] for r in rows}, {'Not assessed'})
            # Salmon counts use replicate columns; TPM/length use processing-run
            # columns and must be mapped through the preflight collapse table.
            salmon = root / 'salmon.tsv.gz'
            tpm = root / 'tpm.tsv.gz'
            lengths = root / 'length.tsv.gz'
            for path, kind in ((salmon, 'counts'), (tpm, 'tpm'), (lengths, 'length')):
                with gzip.open(path, 'wt') as stream:
                    columns = samples if kind == 'counts' else ['run_'+s for s in samples]
                    stream.write('gene_id\t' + '\t'.join(columns) + '\n')
                    for i in range(100):
                        values = [10+i*2, 12+i*3, 18+i, 20+i*4] if kind != 'length' else [1000]*4
                        stream.write('g%d\t%s\n' % (i, '\t'.join(map(str, values))))
            preflight.write_text(json.dumps({'group_a': 'Treatment', 'group_b': 'Control',
                'replicates': {'a': samples[:2], 'b': samples[2:]},
                'collapse': {'run_'+s: s for s in samples}}))
            salmon_command = [RSCRIPT, str(ROOT / 'src/expression_pca.R'), '--route', 'tximport',
                              '--preflight', str(preflight), '--counts', str(salmon),
                              '--tpm', str(tpm), '--length', str(lengths), '--out', str(root / 'salmon_pca')]
            result = subprocess.run(salmon_command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((root / 'salmon_pca.pdf').read_bytes().startswith(b'%PDF'))
            # Both count sources must produce the complete seven-page report,
            # without changing the retained matrices or classifying samples.
            import re
            for report_command, prefix in ((command, 'pca'), (salmon_command, 'salmon_pca')):
                before = counts.read_bytes()
                report = subprocess.run(report_command + ['--qc-report', 'true', '--cooks', 'true'],
                                        capture_output=True, text=True, cwd=root)
                self.assertEqual(report.returncode, 0, report.stderr)
                self.assertEqual(len(re.findall(rb'/Type /Page\b', (root / (prefix+'.pdf')).read_bytes())), 7)
                with gzip.open(root / (prefix+'.summary.tsv.gz'), 'rt') as stream:
                    summaries = list(csv.DictReader(stream, delimiter='\t'))
                self.assertEqual([r['sample_id'] for r in summaries], samples)
                self.assertEqual({r['processing_qc_status'] for r in summaries}, {'Unavailable'})
                self.assertEqual({r['cooks_status'] for r in summaries}, {'Available'})
                with gzip.open(root / (prefix+'.scores.tsv.gz'), 'rt') as stream:
                    self.assertEqual({r['qc_status'] for r in csv.DictReader(stream, delimiter='\t')}, {'Not classified'})
                self.assertEqual(counts.read_bytes(), before)
                self.assertFalse((root / 'Rplots.pdf').exists())
            preflight.write_text(json.dumps({'group_a': 'Treatment', 'group_b': 'Control',
                'replicates': {'a': ['missing', 'a2'], 'b': samples[2:]}}))
            failed = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn('missing', failed.stderr)

    def test_r_distance_flags_are_group_specific_and_degenerate_groups_unassessed(self):
        if not RSCRIPT or not Path(RSCRIPT).is_file():
            self.skipTest('set MICROEXONATOR_RSCRIPT for R integration')
        code = '''
source(commandArgs(TRUE)[1])
x <- rbind(c(0, 1, 2, 3, 4, 100, 1000, 1001, 1002, 1003, 1004, 1100))
q <- group_distance_qc(x, rep(c("A", "B"), each=6))
stopifnot(identical(q$qc_status, rep(c(rep("Unflagged", 5), "Review"), 2)))
q <- group_distance_qc(matrix(1, nrow=2, ncol=6), rep("A", 6))
stopifnot(all(q$qc_status == "Not assessed"), all(is.na(q$distance_threshold)))
'''
        result = subprocess.run([RSCRIPT, '-e', code, str(ROOT / 'src/expression_pca.R')],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
