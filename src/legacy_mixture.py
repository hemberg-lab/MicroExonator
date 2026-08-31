"""Validation helpers for the historical Gaussian-mixture filter."""

import csv


LEGACY_MIXTURE_FAILURE = """\
The legacy Gaussian-mixture filter could not fit two U2-score populations.
This commonly occurs because too few microexons were quantified or their U2
scores lack sufficient variation. Use `filter_method: robustness`, which
applies reproducibility and spanning-read support filters without fitting a
mixture model.
""".strip()


def validate_legacy_mixture_scores(scores):
    """Reject score sets that cannot define two mixture components."""

    scores = list(scores)
    if len(scores) < 2 or len(set(scores)) < 2:
        raise ValueError(LEGACY_MIXTURE_FAILURE)
    return scores


def load_legacy_mixture_scores(
    coverage_file, matches_file, min_number_files_detected
):
    """Load the U2 scores that the historical R mixture receives."""

    covered_microexons = set()
    with open(coverage_file) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["N_samples"]) >= min_number_files_detected:
                covered_microexons.add(row["ME"])

    scores = []
    seen_max_scores = set()
    with open(matches_file) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            max_score = row["ME_max_U2"]
            if row["ME"] in covered_microexons and max_score not in seen_max_scores:
                seen_max_scores.add(max_score)
                scores.append(float(row["U2_score"]))
    return scores
