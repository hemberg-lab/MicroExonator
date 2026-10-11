"""Discover comparison PCA inputs without importing any read-processing rules."""
from pathlib import Path

ROUTES = {'featurecounts': ('featurecounts',),
          'tximport': ('salmon_counts', 'salmon_tpm', 'salmon_length')}


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
