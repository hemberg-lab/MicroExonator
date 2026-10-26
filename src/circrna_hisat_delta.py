"""Experimental boundary-wise circRNA usage from retained HISAT2 junction shards.

Each boundary is tested separately: never sum donor and acceptor evidence.
HISAT2 introns are 1-based closed; circRNA start/end are 1-based exon bases.
This is a cross-aligner junction-support ratio, not native CIRIquant PSI.
"""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    from src.me_delta import event_rng, posterior_draws, group_posterior
except ImportError:
    from me_delta import event_rng, posterior_draws, group_posterior


def read_rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def boundary_counts(paths, units, sites=None):
    """Index unique-fragment counts at genomic left/right exon boundaries.

Duplicate run columns across shards are refused, including identical copies:
they could otherwise silently double-count a sample after regrouping.
"""
    seen = set()
    counts = defaultdict(int)
    for path in paths:
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as stream:
            reader = csv.DictReader(stream, delimiter='\t')
            required = {'chrom', 'start', 'end', 'strand', 'max_anchor', 'multi_total', 'short_anchor_total', 'label'}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError('not a retained junction shard: ' + str(path))
            present = set(reader.fieldnames) & set(units)
            if seen & present:
                raise ValueError('duplicate junction runs: ' + ','.join(sorted(seen & present)))
            seen.update(present)
            keys = set()
            for row in reader:
                key = (row['chrom'], int(row['start']), int(row['end']), row['strand'])
                if key in keys:
                    raise ValueError('duplicate intron in shard: ' + str(key))
                keys.add(key)
                if key[1] > key[2] or key[3] not in ('+', '-'):
                    raise ValueError('invalid junction coordinates/strand')
                left = (key[0], key[2] + 1, key[3], 'left')
                right = (key[0], key[1] - 1, key[3], 'right')
                wanted = [site for site in (left, right) if sites is None or site in sites]
                if not wanted:
                    continue
                for unit in present:
                    n = int(row[unit])
                    if n < 0:
                        raise ValueError('negative junction count')
                    # Upstream linear intron ending before the circle's left exon;
                    # downstream linear intron starting after its right exon.
                    if n:
                        for site in wanted:
                            counts[(unit,) + site] += n
    if seen != set(units):
        raise ValueError('missing HISAT2 runs: ' + ','.join(sorted(set(units) - seen)))
    return counts


def delta(a, b, seed=123456, event='event', size=10000, min_reads=5, min_samples=2, threshold=0.1):
    if size < 100 or min_samples < 2 or min_reads < 1 or not 0 <= threshold <= 1:
        raise ValueError('invalid Bayesian settings')
    used = []
    for group in (a, b):
        values = []
        for circular, linear in group:
            if circular is None or linear is None:
                continue
            if circular < 0 or linear < 0:
                raise ValueError('negative evidence')
            n = circular + linear
            if n >= min_reads:
                values.append((circular / n, n))
        # Stable draw assignment when manifest/sample order changes.
        used.append(sorted(values))
    if any(len(group) < min_samples for group in used):
        return {'status': 'insufficient_coverage', 'samples_a': len(used[0]), 'samples_b': len(used[1])}
    rng = event_rng(seed, event)
    sample_draws = [[posterior_draws(q, n, rng, size) for q, n in group] for group in used]
    means, draws = [], []
    for group in sample_draws:
        mean, draw = group_posterior(group, rng, size)
        means.append(mean)
        draws.append(draw)
    difference = draws[0] - draws[1]
    return dict(status='experimental', ratio_a=means[0], ratio_b=means[1], delta=means[0]-means[1],
                probability=float(max(np.mean(difference > 0), np.mean(difference < 0))),
                probability_increase=float(np.mean(difference > threshold)),
                probability_decrease=float(np.mean(difference < -threshold)),
                lower=float(np.quantile(difference, .025)), upper=float(np.quantile(difference, .975)),
                samples_a=len(used[0]), samples_b=len(used[1]))


def prepare(tool, preflight, source, out, shards=(), settings=None):
    settings = settings or {}
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    reps = preflight['replicates']['a'] + preflight['replicates']['b']
    if len(set(reps)) != len(reps):
        raise ValueError('a replicate occurs on both sides')
    if tool == 'ciriquant':
        rows = list(read_rows(Path(source) / 'bsj_counts.tsv.gz'))
        info = {r['feature_id']: r for r in read_rows(Path(source) / 'results.tsv.gz')}
    else:
        rows = list(read_rows(Path(source) / 'support.tsv.gz'))
        info = {r['circ_id']: r for r in rows}
    matrices = []
    seen = set()
    for row in rows:
        cid = row['circ_id']
        if cid in seen:
            raise ValueError('duplicate circRNA identity (check strand collisions): ' + cid)
        seen.add(cid)
        if cid not in info:
            continue
        values = []
        for rep in reps:
            value = row[rep]
            values.append(None if value in ('NA', '', 'nan') else int(value))
        if any(value is not None and value < 0 for value in values):
            raise ValueError('negative BSJ count')
        matrices.append((cid, values))
    with gzip.open(out / 'counts.tsv.gz', 'wt') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['circ_id'] + reps)
        for cid, values in matrices:
            writer.writerow([cid] + ['NA' if v is None else v for v in values])
    with open(out / 'samples.tsv', 'w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['sample', 'group'])
        for side in ('a', 'b'):
            writer.writerows((r, side) for r in preflight['replicates'][side])
    result = []
    evidence_rows = []
    if shards:
        sites = set()
        for cid, _ in matrices:
            row = info[cid]
            sites.add((row['chrom'], int(row['start']), row['strand'], 'left'))
            sites.add((row['chrom'], int(row['end']), row['strand'], 'right'))
        linear = boundary_counts(shards, preflight['collapse'], sites)
        members = defaultdict(list)
        for unit, rep in preflight['collapse'].items():
            members[rep].append(unit)
        for cid, values in matrices:
            row = info[cid]
            chrom, start, end, strand = row['chrom'], int(row['start']), int(row['end']), row['strand']
            if strand not in ('+', '-'):
                raise ValueError('unknown circle strand: ' + cid)
            for boundary, position in (('left', start), ('right', end)):
                name = ('acceptor' if strand == '+' else 'donor') if boundary == 'left' else ('donor' if strand == '+' else 'acceptor')
                evidence = {r: (v, sum(linear.get((u, chrom, position, strand, boundary), 0) for u in members[r]))
                            for r, v in zip(reps, values)}
                for rep, (circular, linear_count) in evidence.items():
                    evidence_rows.append((cid, strand, name, rep, 'NA' if circular is None else circular, linear_count))
                stats = delta([evidence[r] for r in preflight['replicates']['a']],
                              [evidence[r] for r in preflight['replicates']['b']],
                              event='{}:{}:{}'.format(cid, strand, name), **settings)
                result.append(dict(circ_id=cid, strand=strand, boundary=name, **stats))
    columns = ['circ_id', 'strand', 'boundary', 'status', 'ratio_a', 'ratio_b', 'delta', 'probability',
               'probability_increase', 'probability_decrease', 'lower', 'upper', 'samples_a', 'samples_b']
    with gzip.open(out / 'hisat_delta.tsv.gz', 'wt') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(result)
    with gzip.open(out / 'hisat_evidence.tsv.gz', 'wt') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['circ_id', 'strand', 'boundary', 'replicate', 'bsj', 'linear_unique'])
        writer.writerows(evidence_rows)
    (out / 'preparation.json').write_text(json.dumps(dict(tool=tool, settings=settings,
        endpoint='BSJ/(BSJ+unique HISAT2 linear junction support), separate splice boundaries',
        bayesian_status='experimental_cross_aligner' if shards else 'disabled',
        abundance_universe='native retained calls; detection/annotation filtering applies',
        caveat='Junction aggregates cannot deduplicate across distinct junctions or aligners; probabilities are not calibrated FDR.'), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tool', choices=['ciriquant', 'circexplorer2', 'dcc'], required=True)
    parser.add_argument('--preflight', required=True)
    parser.add_argument('--source', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--shard', action='append', default=[])
    parser.add_argument('--settings', default='{}')
    args = parser.parse_args()
    prepare(args.tool, json.loads(Path(args.preflight).read_text()), args.source, args.out,
            args.shard, json.loads(args.settings))
