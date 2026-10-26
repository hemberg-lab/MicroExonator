"""Small count-level diagnostic, not a benchmark of native circRNA detectors.

Print held-out null/signal performance and agreement for two simulated counters.
The shared-bias case illustrates why agreement cannot replace known truth.
"""
import json
import numpy as np
from src.circrna_hisat_delta import delta


def simulate(seed=63721, events=120):
    rng = np.random.default_rng(seed)
    report = {}
    for scenario in ('null', 'signal', 'moderate_signal', 'heterogeneous_null', 'novel_only', 'shared_bias'):
        calls = [set(), set()]
        eligible = [set(), set()]
        for event in range(events):
            means = (.2, .7) if scenario in ('signal', 'shared_bias', 'novel_only') else (.35, .35)
            if scenario == 'moderate_signal':
                means = (.3, .5)
            concentration = 8 if scenario == 'heterogeneous_null' else 80
            base = [[rng.beta(p*concentration, (1-p)*concentration) for _ in range(3)] for p in means]
            for detector in (0, 1):
                if scenario == 'novel_only' and detector == 1 and event % 2:
                    continue
                eligible[detector].add(event)
                evidence = []
                for group in base:
                    group_counts = []
                    for p in group:
                        n = int(rng.integers(30, 201))
                        b = int(rng.binomial(n, p))
                        group_counts.append((b, n-b))
                    evidence.append(group_counts)
                stats = delta(*evidence, event=str(event), size=2000)
                if max(stats.get('probability_increase', 0), stats.get('probability_decrease', 0)) >= .95:
                    calls[detector].add(event)
        union = calls[0] | calls[1]
        common = eligible[0] & eligible[1]
        common_union = union & common
        report[scenario] = dict(events=events, calls=[len(x) for x in calls],
                                consensus=len(calls[0] & calls[1]),
                                jaccard=len(calls[0] & calls[1])/len(union) if union else None,
                                jointly_measurable=len(common),
                                common_jaccard=len(calls[0] & calls[1] & common)/len(common_union) if common_union else None,
                                truth='no biological change' if scenario in ('null', 'heterogeneous_null', 'shared_bias') else 'usage change')
    return report


if __name__ == '__main__':
    print(json.dumps(simulate(), indent=2))
