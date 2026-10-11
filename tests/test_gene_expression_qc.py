"""Source selection and cached-only gene expression QC integration."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class GeneExpressionQCTests(unittest.TestCase):
    def test_source_selection_supports_none_one_both_and_aliases(self):
        from src import expression_pca_cache as cache
        self.assertTrue(hasattr(cache, 'qc_routes'), 'QC source selection not implemented')
        self.assertEqual(cache.qc_routes([]), [])
        self.assertEqual(cache.qc_routes('none'), [])
        self.assertEqual(cache.qc_routes('featurecounts'), ['featurecounts'])
        self.assertEqual(cache.qc_routes('salmon'), ['tximport'])
        self.assertEqual(cache.qc_routes('featurecounts,salmon'), ['featurecounts', 'tximport'])
        self.assertEqual(cache.qc_routes(['featurecounts', 'salmon']), ['featurecounts', 'tximport'])
        with self.assertRaises(ValueError):
            cache.qc_routes('typo')

    def test_disabled_qc_needs_no_caches_and_enabled_qc_cannot_acquire_reads(self):
        self.assertTrue((ROOT / 'GeneExpressionQC.smk').exists(), 'QC workflow missing')
        try:
            import snakemake
        except ImportError:
            self.skipTest('Snakemake 7 required')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'src').symlink_to(ROOT / 'src', target_is_directory=True)
            (root / 'GeneExpressionQC.smk').write_text((ROOT / 'GeneExpressionQC.smk').read_text())
            cmd = [sys.executable, '-c',
                "import appdirs,runpy; appdirs.system='linux'; runpy.run_module('snakemake',run_name='__main__')",
                '-s', 'GeneExpressionQC.smk', '-n', '--cores', '1', '--nocolor']
            env = dict(os.environ, XDG_CACHE_HOME=str(root / '.cache'))
            disabled = subprocess.run(cmd, cwd=root, env=env, capture_output=True, text=True)
            self.assertEqual(disabled.returncode, 0, disabled.stdout + disabled.stderr)
            self.assertNotIn('cached_gene_expression_qc', disabled.stdout + disabled.stderr)
            enabled = subprocess.run(cmd + ['--config', 'umbrella_expression_qc_sources=featurecounts'],
                                     cwd=root, env=env, capture_output=True, text=True)
            self.assertNotEqual(enabled.returncode, 0)
            self.assertIn('cached', (enabled.stdout + enabled.stderr).lower())
            self.assertNotIn('umbrella_stage_reads', enabled.stdout + enabled.stderr)
            import json
            comparison = root / 'comparisons/ref/project/contrast'
            (comparison / 'joined').mkdir(parents=True)
            (comparison / 'preflight.json').write_text(json.dumps({'shards': ['ref/project/group/batch']}))
            for kind in ('featurecounts', 'salmon_counts', 'salmon_tpm', 'salmon_length'):
                (comparison / 'joined' / (kind+'.tsv.gz')).write_bytes(b'cached')
            both = subprocess.run(cmd + ['--config', 'umbrella_expression_qc_sources=featurecounts,salmon'],
                                  cwd=root, env=env, capture_output=True, text=True)
            self.assertEqual(both.returncode, 0, both.stdout + both.stderr)
            import re
            self.assertRegex(both.stdout + both.stderr, r'cached_gene_expression_qc\s+2')
            self.assertNotIn('umbrella_stage_reads', both.stdout + both.stderr)

    def test_processing_qc_discovery_uses_only_existing_shards(self):
        from src.expression_pca_cache import retained_qc_inputs
        import json
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            preflight = root / 'preflight.json'
            preflight.write_text(json.dumps({'shards': ['ref/project/A/batch', 'ref/project/B/batch']}))
            qc = root / 'qc/ref/project/A/batch.qc.tsv.gz'
            qc.parent.mkdir(parents=True)
            qc.write_bytes(b'cached')
            self.assertEqual(retained_qc_inputs(root / 'comparisons', preflight), [str(qc)])


if __name__ == '__main__':
    unittest.main()
