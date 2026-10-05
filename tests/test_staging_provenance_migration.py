"""Migration must back up exact staging records and refuse semantic changes."""
import base64
import json
import tempfile
import unittest
from pathlib import Path

try:
    from src import staging_provenance
except ImportError:
    staging_provenance = None
from src.umbrella_manifest import load_umbrella_manifest


class MigrationTests(unittest.TestCase):
    def fixture(self, root):
        self.assertIsNotNone(staging_provenance, 'guarded staging migration is missing')
        rules = root / 'rules'
        rules.mkdir()
        (rules / 'umbrella_inputs.smk').write_text((Path(__file__).resolve().parents[1] / 'rules/umbrella_inputs.smk').read_text())
        codes = staging_provenance.staging_fingerprints(rules / 'umbrella_inputs.smk')
        header = 'sample_id\trun_id\tbiological_replicate_id\tproject_id\tbatch_id\tgroup\tsource_type\tsource_1\tsource_2\tlayout\tstrandedness\treference_id\tinclude\n'
        text = header + 'old\tSRR1\told\tP\tb1\tA\tsra\tSRR1\t\tSE\tauto\tref\ttrue\n'
        prior = root / 'prior.tsv'
        proposed = root / 'new.tsv'
        prior.write_text(text)
        proposed.write_text(text)
        manifest = load_umbrella_manifest(proposed)
        output = 'umbrella/work/ref/P/b1/SRR1/R1.fastq.gz'
        metadata = root / '.snakemake/metadata'
        metadata.mkdir(parents=True)
        name = base64.urlsafe_b64encode(output.encode()).decode()
        record = dict(rule='umbrella_stage_single', params=sorted([
            repr(''), repr(str(prior)), repr(manifest.run_sha256('SRR1'))]),
            code=codes[('legacy', 'umbrella_stage_single')], starttime=12, input=[], conda_env='env', incomplete=False,
            shellcmd='python3 src/umbrella_stage_reads.py --manifest ' + str(prior)
                     + ' --run-id SRR1 --r1 ' + output + ' --threads 4 ')
        path = metadata / name
        path.write_text(json.dumps(record))
        valid = root / 'umbrella/work/ref/P/b1/SRR1/reads.valid'
        valid.parent.mkdir(parents=True)
        valid.write_text('validated')
        return prior, proposed, path, record

    def test_preview_no_write_apply_backup_and_idempotence(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prior, proposed, path, record = self.fixture(root)
            preview = staging_provenance.migrate(root, prior, proposed)
            self.assertEqual(preview['records'], 1)
            self.assertEqual(json.loads(path.read_text()), record)
            applied = staging_provenance.migrate(root, prior, proposed, apply=True)
            self.assertEqual(json.loads((Path(applied['backup']) / path.name).read_text()), record)
            changed = json.loads(path.read_text())
            self.assertEqual(changed['params'], sorted(p for p in record['params'] if p != repr(str(prior))))
            self.assertEqual(changed['code'], staging_provenance.staging_fingerprints(
                root / 'rules/umbrella_inputs.smk')[('current', 'umbrella_stage_single')])
            self.assertEqual(changed['conda_env'], 'env')
            self.assertEqual(changed['starttime'], 12)
            self.assertEqual(staging_provenance.migrate(root, prior, proposed, apply=True)['records'], 0)

    def test_refuse_hash_mismatch_without_modifying_record(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prior, proposed, path, record = self.fixture(root)
            record['params'][0] = repr('unknown-hash')
            path.write_text(json.dumps(record))
            with self.assertRaises(ValueError):
                staging_provenance.migrate(root, prior, proposed, apply=True)
            self.assertEqual(json.loads(path.read_text()), record)

    def test_older_manifest_location_certified_from_recorded_command(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prior, proposed, path, record = self.fixture(root)
            older = root / 'first-batch.tsv'
            older.write_text(prior.read_text())
            record['params'] = [repr(str(older)) if p == repr(str(prior)) else p for p in record['params']]
            record['shellcmd'] = record['shellcmd'].replace(str(prior), str(older))
            path.write_text(json.dumps(record))
            self.assertEqual(staging_provenance.migrate(root, prior, proposed)['records'], 1)

    def test_refuse_changed_processing_or_locked_workspace(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prior, proposed, path, record = self.fixture(root)
            proposed.write_text(proposed.read_text().replace('sra\tSRR1', 'sra\tSRR2'))
            with self.assertRaises(ValueError):
                staging_provenance.migrate(root, prior, proposed, apply=True)
            proposed.write_text(prior.read_text())
            locks = root / '.snakemake/locks'
            locks.mkdir()
            (locks / '0.input.lock').write_text('locked')
            with self.assertRaises(ValueError):
                staging_provenance.migrate(root, prior, proposed, apply=True)
            self.assertEqual(json.loads(path.read_text()), record)
