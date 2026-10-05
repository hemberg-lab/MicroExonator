"""Real staging/validation rules must reuse completed libraries across manifest paths."""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class ManifestRelocationTests(unittest.TestCase):
    def command(self, directory, *args):
        return subprocess.run([sys.executable, '-c',
            "import appdirs,runpy; appdirs.system='linux'; runpy.run_module('snakemake',run_name='__main__')",
            '-s', 'Snakefile', '-j', '1', '--nocolor', '--force-use-threads', *args], cwd=directory,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            env=dict(os.environ, XDG_CACHE_HOME=str(Path(directory) / '.cache')))

    def test_relocated_manifest_reuses_protected_validation_and_new_library_only(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for name in ['src', 'envs', 'rules']:
                (root / name).symlink_to(ROOT / name, target_is_directory=True)
            for mate in ['one', 'two']:
                (root / (mate + '.fq')).write_text('@read/1\nACGT\n+\nIIII\n')
            header = ['sample_id','run_id','biological_replicate_id','project_id','batch_id','group',
                      'source_type','source_1','source_2','layout','strandedness','reference_id','include']
            row = ['old','old','old','P','batch1','A','fastq','one.fq','','SE','auto','ref','true']
            def write(path, rows):
                with path.open('w') as handle:
                    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
                    writer.writerow(header)
                    writer.writerows(rows)
            write(root / 'manifest.tsv', [row])
            (root / 'selection.txt').write_text('manifest.tsv')
            (root / 'Snakefile').write_text('''
from pathlib import Path
from src.umbrella_manifest import load_umbrella_manifest
UMBRELLA_MANIFEST = load_umbrella_manifest(Path('selection.txt').read_text().strip())
config = {}
def str2bool(value):
    return str(value).lower() in ('true','t','1')
rule all:
    input: [r.work_dir + '/reads.valid' for r in UMBRELLA_MANIFEST.included_runs()]
include: 'rules/umbrella_inputs.smk'
''')
            first = self.command(root)
            self.assertEqual(first.returncode, 0, first.stdout)
            valid = root / 'umbrella/work/ref/P/batch1/old/reads.valid'
            before = valid.stat().st_mtime_ns
            write(root / 'relocated.tsv', [row])
            (root / 'selection.txt').write_text('relocated.tsv')
            moved = self.command(root, '-n')
            self.assertEqual(moved.returncode, 0, moved.stdout)
            self.assertIn('Nothing to be done', moved.stdout)
            new = ['new','new','new','P','batch1','B','fastq','two.fq','','SE','auto','ref','true']
            write(root / 'relocated.tsv', [row, new])
            added = self.command(root, '-n')
            self.assertEqual(added.returncode, 0, added.stdout)
            self.assertIn('run_id=new', added.stdout)
            self.assertNotIn('run_id=old', added.stdout)
            self.assertEqual(valid.stat().st_mtime_ns, before)
            changed = list(row)
            changed[7] = 'two.fq'
            write(root / 'relocated.tsv', [changed])
            altered = self.command(root, '-n')
            self.assertNotEqual(altered.returncode, 0)
            self.assertIn('ProtectedOutputException', altered.stdout)

    def test_migration_of_real_legacy_paired_metadata_preserves_outputs(self):
        from src.staging_provenance import migrate
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for name in ['src', 'envs']:
                (root / name).symlink_to(ROOT / name, target_is_directory=True)
            (root / 'rules').mkdir()
            current = (ROOT / 'rules/umbrella_inputs.smk').read_text()
            legacy = current.replace('{UMBRELLA_MANIFEST.path:q}', '{params.manifest:q}').replace(
                'tmpdir=" ".join(filter(None, (UMBRELLA_STAGE_TMPDIR, UMBRELLA_ALLOW_RESTAGE))),\n',
                'tmpdir=" ".join(filter(None, (UMBRELLA_STAGE_TMPDIR, UMBRELLA_ALLOW_RESTAGE))),\n        manifest=str(UMBRELLA_MANIFEST.path),\n')
            intake = root / 'rules/umbrella_inputs.smk'
            intake.write_text(legacy)
            for mate in [1, 2]:
                (root / ('r%d.fq' % mate)).write_text('@read/%d\nACGT\n+\nIIII\n' % mate)
            text = 'sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\nold\told\told\tP\tb1\tA\tfastq\tr1.fq\tr2.fq\tPE\tauto\tref\ttrue\n'
            (root / 'manifest.tsv').write_text(text)
            (root / 'new.tsv').write_text(text)
            snake = '''from pathlib import Path
from src.umbrella_manifest import load_umbrella_manifest
UMBRELLA_MANIFEST = load_umbrella_manifest(Path('selection.txt').read_text().strip())
config = {}
def str2bool(value):
    return str(value).lower() in ('true','t','1')
rule all:
    input: [r.work_dir + '/reads.valid' for r in UMBRELLA_MANIFEST.included_runs()]
include: 'rules/umbrella_inputs.smk'
'''
            (root / 'Snakefile').write_text(snake)
            (root / 'selection.txt').write_text('manifest.tsv')
            first = self.command(root)
            self.assertEqual(first.returncode, 0, first.stdout)
            valid = root / 'umbrella/work/ref/P/b1/old/reads.valid'
            before = valid.stat().st_mtime_ns
            intake.write_text(current)
            receipt = migrate(root, root / 'manifest.tsv', root / 'new.tsv', apply=True)
            self.assertEqual(receipt['records'], 2)
            (root / 'selection.txt').write_text('new.tsv')
            after = self.command(root, '-n')
            self.assertEqual(after.returncode, 0, after.stdout)
            self.assertIn('Nothing to be done', after.stdout)
            self.assertEqual(valid.stat().st_mtime_ns, before)
            intake.write_text(current.replace('--threads {threads}', '--threads 2'))
            future_change = self.command(root, '-n')
            self.assertNotEqual(future_change.returncode, 0, future_change.stdout)
            self.assertIn('ProtectedOutputException', future_change.stdout)


if __name__ == '__main__':
    unittest.main()
