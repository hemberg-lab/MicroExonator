"""Explicit, backed-up migration of path-only staging provenance (Snakemake 7)."""
import argparse
import ast
import base64
import io
import json
import os
import shlex
import tempfile
from pathlib import Path
from types import SimpleNamespace

if __package__:
    from src.umbrella_manifest import load_umbrella_manifest
    from src.umbrella_identity import check_existing_shards
else:
    from umbrella_manifest import load_umbrella_manifest
    from umbrella_identity import check_existing_shards


def staging_fingerprints(path):
    """Compile only staging function bodies; never execute workflow or stored pickle."""
    from snakemake.parser import parse
    from snakemake.persistence import pickle_code
    from snakemake.sourcecache import LocalSourceFile
    text = Path(path).read_text()
    fingerprints = {}
    for version, source in [('current', text), ('legacy', text.replace(
            '{UMBRELLA_MANIFEST.path:q}', '{params.manifest:q}'))]:
        reader = SimpleNamespace(sourcecache=SimpleNamespace(open=lambda unused: io.StringIO(source)))
        parsed = parse(LocalSourceFile(str(path)), reader)[0]
        for node in ast.parse(parsed).body:
            if isinstance(node, ast.FunctionDef) and node.name in {
                    '__rule_umbrella_stage_reads', '__rule_umbrella_stage_single'}:
                node.decorator_list = []
                namespace = {}
                exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
                code = namespace[node.name].__code__
                fingerprints[(version, node.name.removeprefix('__rule_'))] = base64.b64encode(pickle_code(code)).decode()
    if len(fingerprints) != 4:
        raise ValueError('Cannot certify staging fingerprints for this Snakemake parser.')
    return fingerprints


def migrate(root, prior, proposed, apply=False, scratch_options=''):
    root, prior, proposed = Path(root).resolve(), Path(prior), Path(proposed)
    locks = root / '.snakemake/locks'
    if locks.exists() and any(locks.iterdir()):
        raise ValueError('Workspace has Snakemake locks; wait for its driver/jobs to finish.')
    old, new = load_umbrella_manifest(prior), load_umbrella_manifest(proposed)
    for run in old.included_runs():
        if run.run_id not in new.by_run or not new.by_run[run.run_id].include:
            raise ValueError('Existing run removed/excluded: ' + run.run_id)
        if old.run_sha256(run.run_id) != new.run_sha256(run.run_id):
            raise ValueError('Processing identity changed: ' + run.run_id)
        if old.shard_sha256(run.reference_id, run.project_id, run.group, run.batch_id) != new.shard_sha256(
                run.reference_id, run.project_id, run.group, run.batch_id):
            raise ValueError('Stored group/shard changed: ' + run.run_id)
    problems = check_existing_shards(new, root)
    if problems:
        raise ValueError('\n'.join(problems))
    metadata = root / '.snakemake/metadata'
    changes = []
    fingerprints = None
    for path in sorted(metadata.glob('*')):
        if not path.is_file():
            continue
        try:
            record = json.loads(path.read_text())
        except (ValueError, UnicodeError):
            continue
        if record.get('rule') not in {'umbrella_stage_reads', 'umbrella_stage_single'}:
            continue
        params = record.get('params', [])
        try:
            output = base64.urlsafe_b64decode(path.name).decode()
        except (ValueError, UnicodeError):
            raise ValueError('Unrecognised metadata filename: ' + str(path))
        parts = output.split('/')
        if len(parts) != 7 or parts[:2] != ['umbrella', 'work'] or parts[-1] not in {'R1.fastq.gz','R2.fastq.gz'}:
            raise ValueError('Unrecognised staging output: ' + output)
        run_id = parts[-2]
        if run_id not in old.by_run:
            continue  # not an existing library of this reviewed run
        run = new.by_run[run_id]
        tokens = shlex.split(record.get('shellcmd', ''))
        if len(tokens) < 4 or tokens[:3] != ['python3','src/umbrella_stage_reads.py','--manifest']:
            raise ValueError('Recorded staging command not recognised: ' + run_id)
        recorded_path = tokens[3]
        if repr(recorded_path) not in params:
            continue  # path parameter already removed; retain all other provenance
        historical_path = Path(recorded_path)
        if not historical_path.is_absolute():
            historical_path = root / historical_path
        historical = load_umbrella_manifest(historical_path)
        if run_id not in historical.by_run or historical.run_sha256(run_id) != new.run_sha256(run_id):
            raise ValueError('Historical manifest does not certify this processing identity: ' + run_id)
        expected = sorted([repr(scratch_options), repr(recorded_path), repr(new.run_sha256(run_id))])
        if sorted(params) != expected or record.get('incomplete'):
            raise ValueError('Staging hash/settings differ or incomplete: ' + run_id)
        if fingerprints is None:
            fingerprints = staging_fingerprints(root / 'rules/umbrella_inputs.smk')
        if record.get('code') != fingerprints[('legacy', record['rule'])]:
            raise ValueError('Recorded staging code is not the reviewed legacy template: ' + run_id)
        if output != run.work_dir + '/' + parts[-1] or (root / run.work_dir / 'reads.valid').is_file() is False:
            raise ValueError('Completed validation/path missing: ' + run_id)
        paired = record['rule'] == 'umbrella_stage_reads'
        if paired != (run.layout == 'PE'):
            raise ValueError('Recorded rule/layout disagree: ' + run_id)
        expected_tokens = ['python3','src/umbrella_stage_reads.py','--manifest',recorded_path,
                           '--run-id',run_id,'--r1',run.work_dir + '/R1.fastq.gz']
        if paired:
            expected_tokens += ['--r2',run.work_dir + '/R2.fastq.gz']
        if tokens[:len(expected_tokens)] != expected_tokens:
            raise ValueError('Recorded staging command not recognised: ' + run_id)
        tail = tokens[len(expected_tokens):]
        if len(tail) < 2 or tail[0] != '--threads' or not tail[1].isdigit() or int(tail[1]) < 1 or tail[2:] != shlex.split(scratch_options):
            raise ValueError('Recorded command options differ: ' + run_id)
        updated = dict(record)
        updated['params'] = sorted(p for p in params if p != repr(recorded_path))
        # Only these exact commands changed their shell placeholder, not processing.
        # Preserve input/environment/timestamps; never clean all provenance or outputs.
        updated['code'] = fingerprints[('current', record['rule'])]
        changes.append((path, record, updated))
    report = dict(records=len(changes), apply=apply, prior_manifest=str(prior),
                  proposed_manifest=str(proposed), backup=None)
    if apply and changes:
        parent = root / 'umbrella/provenance-migrations'
        parent.mkdir(parents=True, exist_ok=True)
        backup = Path(tempfile.mkdtemp(prefix='staging-manifest-path-v1-', dir=parent))
        for path, record, updated in changes:
            (backup / path.name).write_bytes(path.read_bytes())
        report['backup'] = str(backup)
        (backup / 'receipt.json').write_text(json.dumps(report, indent=2) + '\n')
        for path, record, updated in changes:
            # Refuse a race rather than overwrite unexpectedly changed metadata.
            if json.loads(path.read_text()) != record:
                raise ValueError('Metadata changed during migration; backups kept at ' + str(backup))
        for path, record, updated in changes:
            temporary = path.with_name(path.name + '.manifest-path.tmp')
            temporary.write_text(json.dumps(updated))
            os.replace(temporary, path)
    return report


def main():
    import yaml
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior-config', required=True)
    parser.add_argument('--configfile', required=True)
    parser.add_argument('--apply', action='store_true', help='Default is a read-only preview.')
    args = parser.parse_args()
    old = yaml.safe_load(Path(args.prior_config).read_text())
    new = yaml.safe_load(Path(args.configfile).read_text())
    keys = set(old) | set(new)
    changed = {key for key in keys if old.get(key) != new.get(key)}
    if changed - {'umbrella_manifest', 'umbrella_comparisons'}:
        raise ValueError('Unexpected processing config changes: ' + str(changed))
    scratch = '--tmpdir ' + shlex.quote(str(new['umbrella_stage_tmpdir'])) if new.get('umbrella_stage_tmpdir') else ''
    if str(new.get('umbrella_allow_restage', False)).lower() in {'true','t','1'}:
        raise ValueError('Do not enable restaging for this migration.')
    print(json.dumps(migrate(Path.cwd(), old['umbrella_manifest'], new['umbrella_manifest'],
                             apply=args.apply, scratch_options=scratch), indent=2))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        raise SystemExit('STOP: ' + str(error))
