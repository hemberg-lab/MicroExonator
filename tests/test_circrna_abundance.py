"""Tiny native edgeR smoke; set CIRCRNA_DIFF_RSCRIPT to the pinned Rscript."""
import csv
import gzip
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np

RSCRIPT = os.environ.get('CIRCRNA_DIFF_RSCRIPT')
ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(RSCRIPT, 'set CIRCRNA_DIFF_RSCRIPT for native edgeR tests')
class AbundanceTests(unittest.TestCase):
    def run_fixture(self, zero=False, rows=200, identical=False, blocks=None, donor_effect=False):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            samples = ['a1', 'a2', 'a3', 'b1', 'b2', 'b3']
            rng = np.random.default_rng(157)
            counts = rng.poisson(100, size=(rows, 6))
            counts[:min(10, rows), 3:] += 250
            if donor_effect:
                # donors differ strongly; the group effect is consistent within each donor
                # (a per-circle donor level, so library-size normalization cannot remove it)
                level = rng.lognormal(0, 1, size=(rows, 3))
                for column in range(6):
                    fold = np.where(np.arange(rows) < 10, 1.6, 1.0) if column >= 3 else 1.0
                    counts[:, column] = rng.poisson(100 * level[:, column % 3] * fold)
            if zero:
                counts[:] = 0
            if identical:
                counts[:] = np.repeat(rng.poisson(50, size=(rows, 1)), 6, axis=1)
            with gzip.open(path / 'counts.tsv.gz', 'wt') as stream:
                writer = csv.writer(stream, delimiter='\t')
                writer.writerow(['circ_id'] + samples)
                for i, row in enumerate(counts):
                    writer.writerow(['circle{}'.format(i)] + list(row))
            with open(path / 'samples.tsv', 'w') as stream:
                stream.write('sample\tgroup' + ('\tblock' if blocks else '') + '\n')
                for index, sample in enumerate(samples):
                    stream.write(sample + '\t' + sample[0] + ('\t' + blocks[index] if blocks else '') + '\n')
            result = subprocess.run([RSCRIPT, str(ROOT / 'src/circrna_abundance.R'), str(path)],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(path / 'abundance.tsv.gz', 'rt') as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            return rows, (path / 'status.tsv').read_text()

    def test_known_abundance_changes(self):
        rows, status = self.run_fixture()
        self.assertIn('status\ttested', status)
        selected = {row['circ_id']: row for row in rows}
        self.assertLess(float(selected['circle0']['logFC']), -1)
        self.assertLess(float(selected['circle0']['FDR']), .05)
        self.assertEqual(len(rows), 200)

    def test_few_events_are_a_status_with_a_reason(self):
        rows, status = self.run_fixture(rows=2)
        self.assertEqual(rows, [])
        self.assertIn('status\tinsufficient_data', status)
        self.assertIn('fewer than 10', status)

    def test_degenerate_replicates_never_crash(self):
        # identical replicates broke glmQLFit with two events on the pilot fixture
        rows, status = self.run_fixture(rows=50, identical=True)
        first = status.splitlines()[0]
        self.assertIn(first, ('status\ttested', 'status\tfit_failed'))
        if first == 'status\ttested':
            self.assertTrue(all(float(row['FDR']) > .05 for row in rows))

    def test_zero_support_is_not_significant(self):
        rows, status = self.run_fixture(zero=True)
        self.assertEqual(rows, [])
        self.assertIn('status\tinsufficient_data', status)

    def test_blocked_design_adds_the_block_term(self):
        donors = ['d1', 'd2', 'd3', 'd1', 'd2', 'd3']
        rows, status = self.run_fixture(blocks=donors, donor_effect=True)
        self.assertIn('status\ttested', status)
        self.assertIn('model\tblocked_glm (~ block + group)', status)
        selected = {row['circ_id']: row for row in rows}
        changed = [selected['circle{}'.format(i)] for i in range(10)]
        self.assertTrue(all(float(row['logFC']) < 0 for row in changed))   # higher in B: log2 A/B < 0
        blocked = sum(float(row['FDR']) < .05 for row in changed)
        unblocked, plain = self.run_fixture(donor_effect=True)
        self.assertNotIn('model', plain)
        # donor variation swamps the unblocked test; the block term recovers the changes
        plain_called = sum(float(r['FDR']) < .05 for r in unblocked if int(r['circ_id'][6:]) < 10)
        self.assertGreater(blocked, plain_called)
        self.assertGreaterEqual(blocked, 8)
