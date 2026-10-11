"""Discover comparison PCA inputs without importing any read-processing rules."""
from pathlib import Path
import json

ROUTES = {'featurecounts': ('featurecounts',),
          'tximport': ('salmon_counts', 'salmon_tpm', 'salmon_length')}


def qc_routes(sources):
    """Public source names, including CLI comma syntax; empty means disabled."""
    if isinstance(sources, str):
        sources = [s.strip().lower() for s in sources.split(',') if s.strip()]
    if sources is None or sources == ['none']:
        return []
    if not isinstance(sources, (list, tuple)):
        raise ValueError('QC sources must be a list or comma-separated source names')
    aliases = {'featurecounts': 'featurecounts', 'salmon': 'tximport', 'tximport': 'tximport'}
    routes = []
    for source in sources:
        if source not in aliases:
            raise ValueError('Unknown QC source: {} (use featurecounts or salmon)'.format(source))
        route = aliases[source]
        if route not in routes:
            routes.append(route)
    return routes


def retained_qc_inputs(root, preflight):
    """Use only available shard QC; never regenerate missing processing metrics."""
    with open(preflight) as stream:
        shards = json.load(stream).get('shards', [])
    inputs = []
    for shard in shards:
        # A preflight's shard identifiers are relative, not arbitrary paths.
        parts = shard.split('/')
        if len(parts) != 4 or any(p in ('', '.', '..') for p in parts):
            raise ValueError('Invalid retained QC shard identifier: ' + shard)
        path = Path(root).parent / 'qc' / (shard + '.qc.tsv.gz')
        if path.is_file():
            inputs.append(str(path))
    return inputs


def cached_comparisons(root, routes, comparisons=None):
    root = Path(root)
    if not routes or len(routes) != len(set(routes)) or set(routes) - set(ROUTES):
        raise ValueError('PCA routes must be a nonempty unique list of featurecounts and/or tximport')
    selected = set(comparisons) if comparisons is not None else None
    jobs = {}
    found = set()
    for preflight in sorted(root.glob('*/*/*/preflight.json')):
        key = preflight.parent.relative_to(root).as_posix()
        if selected is not None and key not in selected:
            continue
        found.add(key)
        for route in routes:
            matrices = [preflight.parent / 'joined' / (kind + '.tsv.gz') for kind in ROUTES[route]]
            missing = [str(path) for path in matrices if not path.is_file()]
            if missing:
                raise ValueError('Cached-only PCA: missing retained matrices: ' + ', '.join(missing) +
                                 '. No acquisition or alignment will be scheduled.')
            jobs[key + '/' + route] = {'preflight': str(preflight),
                                      'matrices': [str(path) for path in matrices]}
    if selected is not None and selected - found:
        raise ValueError('Missing cached comparison preflights: ' + ', '.join(sorted(selected - found)))
    if not jobs:
        raise ValueError('No cached comparison preflights found under ' + str(root))
    return jobs
