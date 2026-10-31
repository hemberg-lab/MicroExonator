"""DESeq2 per comparison: ~ group, or ~ block + group for a paired/blocked design.

Needs an Rscript with DESeq2 and jsonlite (envs/umbrella-comparisons.yaml):
  MICROEXONATOR_RSCRIPT=/path/to/env/bin/Rscript python -m unittest tests.test_expression_deseq2
"""
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get('MICROEXONATOR_RSCRIPT', shutil.which('Rscript') or '')


def has_deseq2():
    if not RSCRIPT:
        return False
    check = subprocess.run([RSCRIPT, '-e', 'quit(status = !(requireNamespace("DESeq2", quietly=TRUE) && '
                                          'requireNamespace("jsonlite", quietly=TRUE)))'],
                           capture_output=True)
    return check.returncode == 0


@unittest.skipUnless(has_deseq2(), 'set MICROEXONATOR_RSCRIPT to an Rscript with DESeq2 and jsonlite')
class DESeq2DesignTests(unittest.TestCase):
    def run_route(self, design):
        replicates = {'a': ['A1', 'A2', 'A3'], 'b': ['C1', 'C2', 'C3']}
        donor = {'A1': 'd1', 'C1': 'd1', 'A2': 'd2', 'C2': 'd2', 'A3': 'd3', 'C3': 'd3'}
        preflight = {'replicates': replicates, 'inference_supported': True, 'reasons': [],
                     'collapse': {r: r for side in replicates.values() for r in side}}
        if design:
            preflight['design'] = {'type': 'paired', 'block': 'donor', 'replicate_block': donor}
        random.seed(5)
        level = {g: {d: random.lognormvariate(0, 1.2) for d in ('d1', 'd2', 'd3')} for g in range(300)}
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'preflight.json').write_text(json.dumps(preflight))
            with gzip.open(root / 'counts.tsv.gz', 'wt') as stream:
                writer = csv.writer(stream, delimiter='\t')
                order = replicates['a'] + replicates['b']
                writer.writerow(['gene_id', 'length'] + order)
                for gene in range(300):
                    fold = 2.0 if gene < 20 else 1.0
                    writer.writerow(['g{}'.format(gene), 1000] + [
                        int(random.gauss(1, 0.03) * 200 * level[gene][donor[r]] * (fold if r[0] == 'A' else 1))
                        for r in order])
            result = subprocess.run([RSCRIPT, str(ROOT / 'src/expression_deseq2.R'), '--route', 'featurecounts',
                                     '--preflight', str(root / 'preflight.json'),
                                     '--counts', str(root / 'counts.tsv.gz'), '--out', str(root / 'out')],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(root / 'out.results.tsv.gz', 'rt') as stream:
                rows = {row['gene_id']: row for row in csv.DictReader(stream, delimiter='\t')}
            return rows, (root / 'out.status.txt').read_text()

    def test_block_term_recovers_paired_changes(self):
        paired, status = self.run_route(design=True)
        self.assertEqual(status.splitlines(), ['supported', 'model: blocked_glm (~ donor + group)'])
        unpaired, plain = self.run_route(design=False)
        self.assertEqual(plain.splitlines(), ['supported'])
        called = lambda rows: sum(rows['g{}'.format(g)]['padj'] not in ('NA', '') and
                                  float(rows['g{}'.format(g)]['padj']) < 0.05 for g in range(20))
        self.assertTrue(all(float(paired['g{}'.format(g)]['log2FoldChange']) > 0.5 for g in range(20)))
        self.assertGreater(called(paired), called(unpaired))


if __name__ == '__main__':
    unittest.main()
