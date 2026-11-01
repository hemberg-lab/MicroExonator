"""Tiny native PAIRADISE smoke, using the script shipped with pinned rMATS."""
import csv
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ENV = os.environ.get('RMATS_PAIRED_ENV')


@unittest.skipUnless(ENV, 'set RMATS_PAIRED_ENV to envs/umbrella-rmats-paired.yaml prefix')
class NativePairedTests(unittest.TestCase):
    def test_native_paired_model_emits_probabilities(self):
        prefix = Path(ENV)
        script = prefix/'rMATS/rMATS_R/paired_model.R'
        self.assertTrue(script.is_file())
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            source = work/'counts.tsv'
            source.write_text('ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n'
                'changed1\t80,72,90,85,77,88\t20,28,10,15,23,12\t20,12,30,25,17,28\t80,88,70,75,83,72\t200\t100\n'
                'changed2\t20,12,30,25,17,28\t80,88,70,75,83,72\t80,72,90,85,77,88\t20,28,10,15,23,12\t200\t100\n'
                'null\t45,55,40,60,47,53\t55,45,60,40,53,47\t47,53,42,58,45,55\t53,47,58,42,55,45\t200\t100\n')
            result = subprocess.run([str(prefix/'bin/Rscript'), str(script), str(source), '1', str(work/'out.tsv')],
                                    cwd=work, capture_output=True, text=True, timeout=120)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with open(work/'out.tsv') as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            self.assertEqual({r['ID'] for r in rows}, {'changed1','changed2','null'})
            for row in rows:
                for column in ('PValue','FDR'):
                    p = float(row[column])
                    self.assertTrue(math.isfinite(p) and 0 <= p <= 1)
