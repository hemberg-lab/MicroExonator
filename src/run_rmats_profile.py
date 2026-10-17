"""Isolated no-novel-splice-site rMATS profile, preserving the original route.

Library strandedness intentionally remains fr-unstranded, as in the original
umbrella rMATS workflow. This profile changes novelty handling, not the library
policy or statistical thresholds. Both profiles belong to the rmats tool family.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path

if __package__:
    from src.umbrella_modules import NORMALIZED, write_status
else:
    from umbrella_modules import NORMALIZED, write_status

KINDS = ('SE', 'A5SS', 'A3SS', 'RI', 'MXE')
EVENT_COORDINATES = {
    'SE': ('exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'),
    'A3SS': ('longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE'),
    'A5SS': ('longExonStart_0base', 'longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE'),
    'RI': ('riExonStart_0base', 'riExonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'),
    'MXE': ('1stExonStart_0base', '1stExonEnd', '2ndExonStart_0base', '2ndExonEnd',
            'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'),
}
POLICY = {'tool': 'rmats_conservative', 'tool_family': 'rmats',
          'novel_splice_sites': False, 'library_type': 'fr-unstranded',
          'variable_read_length': True}


def prep(bam, gtf, read_length, out_dir, run_id, layout='PE', threads=4, identity=None):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    bam = str(Path(bam).resolve())
    if layout not in ('PE', 'SE') or int(read_length) < 1:
        raise ValueError('rMATS prep needs PE/SE and positive read length')
    with tempfile.TemporaryDirectory(prefix='prep-', dir=out) as directory:
        work = Path(directory).resolve()
        (work/'tmp').mkdir(); (work/'od').mkdir()
        (work/'b1.txt').write_text(bam + '\n')
        command = ['rmats.py', '--b1', str(work/'b1.txt'), '--gtf', str(Path(gtf).resolve()),
                   '--od', str(work/'od'), '--tmp', str(work/'tmp'), '-t',
                   'paired' if layout == 'PE' else 'single', '--readLength', str(read_length),
                   '--variable-read-length', '--libType', 'fr-unstranded',
                   '--nthread', str(threads), '--task', 'prep']
        with (out/'prep.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        files = list((work/'tmp').glob('*.rmats'))
        if len(files) != 1:
            raise RuntimeError('rMATS prep must emit exactly one .rmats file')
        shutil.copyfile(files[0], out/'run.rmats')
    record = dict(identity or {}, **POLICY, run_id=run_id, bam=bam, layout=layout,
                  read_length=int(read_length),
                  sha256=hashlib.sha256((out/'run.rmats').read_bytes()).hexdigest())
    (out/'cache.json').write_text(json.dumps(record, sort_keys=True, indent=2)+'\n')


def number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def check_provenance(record, cache, run, preflight, ids):
    """A cache made for another reference, module profile or layout is refused,
    even if it was copied under a valid-looking path (its folder name is the
    cache_id the workflow expects)."""
    expected = {'reference_id': ids.get('reference_id'),
                'module_reference_id': ids.get('module_reference_id'),
                'cache_id': Path(cache).name if ids else None}
    for key, value in expected.items():
        if value is not None and record.get(key) != value:
            raise ValueError('rMATS cache {} of {} has {} {}, expected {}'.format(
                cache, run, key, record.get(key), value))
    layouts = preflight.get('layouts')
    if layouts and record.get('layout') not in layouts:
        raise ValueError('rMATS cache {} of {} is {}, the comparison is {}'.format(
            cache, run, record.get('layout'), '/'.join(layouts)))


def normalized_rows(path, kind, preflight, ids):
    rows = []
    with open(path) as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if kind == 'SE':
                start, end = int(row['exonStart_0base'])+1, int(row['exonEnd'])
            elif kind in ('A3SS', 'A5SS'):
                start, end = int(row['longExonStart_0base'])+1, int(row['longExonEnd'])
            elif kind == 'RI':
                start, end = int(row['upstreamEE'])+1, int(row['downstreamES'])
            else:
                start, end = int(row['1stExonStart_0base'])+1, int(row['2ndExonEnd'])
            p, q, effect = (number(row[key]) for key in ('PValue', 'FDR', 'IncLevelDifference'))
            # Preserve every native coordinate and count, not an interval-only event identity.
            context = {key: value for key, value in row.items()
                       if key not in ('PValue', 'FDR', 'IncLevelDifference', 'IncLevel1', 'IncLevel2')}
            identity = ':'.join([kind, row['GeneID'].strip('"'), row['chr'], row['strand'],
                                 hashlib.sha256(json.dumps({key: row[key] for key in EVENT_COORDINATES[kind]},
                                    sort_keys=True).encode()).hexdigest()[:16]])
            rows.append(dict(ids, tool='rmats_conservative', feature_id=identity,
                gene_id=row['GeneID'].strip('"'), gene_name=row['geneSymbol'].strip('"'),
                chrom=row['chr'], strand=row['strand'], start=start, end=end,
                coordinates=json.dumps(context, sort_keys=True), event_class=kind,
                effect='' if effect is None else effect,
                effect_definition='native rMATS IncLevelDifference, PSI(A) - PSI(B)',
                n_a=len(preflight['replicates']['a']), n_b=len(preflight['replicates']['b']),
                coverage=json.dumps({k: row[k] for k in ('IJC_SAMPLE_1','SJC_SAMPLE_1',
                                                       'IJC_SAMPLE_2','SJC_SAMPLE_2')}, sort_keys=True),
                statistic_type='rMATS likelihood-ratio p-value; native event-class FDR',
                p_value='' if p is None else p, q_value='' if q is None else q,
                status='tested' if p is not None and effect is not None else 'untested'))
    return rows


def compare(preflight, caches, gtf, out_dir, threads=4, ids=None):
    if not isinstance(preflight, dict):
        preflight = json.loads(Path(preflight).read_text())
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    base = dict(ids or {}, **POLICY, project_id=preflight.get('project_id',''),
                comparison_id=preflight.get('comparison_id',''))
    reasons = list(preflight['tools']['rmats'].get('reasons', []))
    if not preflight.get('inference_supported', False) or not preflight['tools']['rmats']['supported']:
        reasons += preflight.get('reasons', []) or ['unsupported comparison design']
    replicates = list(preflight['collapse'].values())
    if len(replicates) != len(set(replicates)):
        reasons.append('multiple processing units per biological replicate; merge technical libraries before rMATS')
    if any(len(preflight['replicates'][side]) < 2 for side in ('a','b')):
        reasons.append('at least two biological replicates per side are required')
    rows = []
    with tarfile.open(out/'native.tar.gz', 'w:gz') as archive:
        if not reasons:
            records = {}
            for run in sorted(preflight['collapse']):
                cache = Path(caches[run])
                record = json.loads((cache/'cache.json').read_text())
                if any(record.get(k) != v for k, v in POLICY.items()):
                    raise ValueError('incompatible conservative rMATS cache: '+run)
                if record.get('run_id') != run or hashlib.sha256((cache/'run.rmats').read_bytes()).hexdigest() != record['sha256']:
                    raise ValueError('rMATS cache identity/checksum mismatch: '+run)
                check_provenance(record, cache, run, preflight, ids or {})
                records[run] = record
            lengths = {r['read_length'] for r in records.values()}
            base['prep_read_lengths'] = sorted(lengths)
            base['warnings'] = (['mixed prep read lengths; post uses maximum, matching original rMATS route']
                                if len(lengths) > 1 else [])
            with tempfile.TemporaryDirectory(prefix='post-', dir=out) as directory:
                work=Path(directory).resolve(); (work/'tmp').mkdir(); (work/'od').mkdir()
                for index, run in enumerate(sorted(records)):
                    shutil.copyfile(Path(caches[run])/'run.rmats', work/'tmp'/('{}-run.rmats'.format(index)))
                for side in ('a','b'):
                    (work/(side+'.txt')).write_text(','.join(records[run]['bam'] for run in preflight['runs'][side])+'\n')
                command=['rmats.py','--b1',str(work/'a.txt'),'--b2',str(work/'b.txt'),
                         '--gtf',str(Path(gtf).resolve()),'--od',str(work/'od'),'--tmp',str(work/'tmp'),
                         '-t','paired' if preflight['layouts']==['PE'] else 'single',
                         '--readLength',str(max(lengths)), '--variable-read-length',
                         '--nthread',str(threads),'--individual-counts']
                with (out/'compare.log').open('w') as log:
                    for task in ('inte','post'):
                        subprocess.run(command+['--task',task],check=True,stdout=log,stderr=subprocess.STDOUT)
                for kind in KINDS:
                    rows += normalized_rows(work/'od'/('{}.MATS.JC.txt'.format(kind)),kind,preflight,base)
                for path in sorted((work/'od').iterdir()):
                    if path.is_file() and path.suffix=='.txt':
                        archive.add(path,arcname=path.name)
    with gzip.open(out/'results.tsv.gz','wt') as stream:
        writer=csv.DictWriter(stream,fieldnames=NORMALIZED,delimiter='\t',extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    tested = sum(row['status'] == 'tested' for row in rows)
    if reasons:
        status = 'unsupported'
    elif not rows:
        status = 'empty'
    elif not tested:
        status = 'no_tests'
        reasons = ['rMATS reported events, but none has a p-value and an effect']
    else:
        status = 'ok'
    write_status(out/'status.json',base.pop('tool'),status,
                 reasons=reasons,events=len(rows),tested=tested,**base)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('prep')
    p.add_argument('--bam',required=True);p.add_argument('--gtf',required=True)
    p.add_argument('--read-length',required=True,type=int);p.add_argument('--run-id',required=True)
    p.add_argument('--layout',choices=('PE','SE'),required=True)
    c=sub.add_parser('compare');c.add_argument('--preflight',required=True)
    c.add_argument('--cache',action='append',default=[]);c.add_argument('--gtf',required=True)
    for command in (p,c):
        command.add_argument('--out-dir',required=True);command.add_argument('--threads',type=int,default=4)
        command.add_argument('--ids',default='{}')
    args=parser.parse_args()
    if args.command=='prep':
        prep(args.bam,args.gtf,args.read_length,args.out_dir,args.run_id,args.layout,args.threads,json.loads(args.ids))
    else:
        compare(args.preflight,dict(item.split('=',1) for item in args.cache),args.gtf,args.out_dir,args.threads,json.loads(args.ids))


if __name__=='__main__':
    main()
