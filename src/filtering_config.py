"""Configuration helpers for selecting and grouping confidence filters."""

import csv
import re
import warnings
from collections import defaultdict


FILTER_METHODS = frozenset(("robustness", "legacy_mixture"))
FILTER_OUTPUTS = {
    "robustness": "Report/out.robustly_detected.txt",
    "legacy_mixture": "Report/out.high_quality.txt",
}
GROUP_ID_PATTERN = re.compile(r"^[A-Za-z0-9_.-]+$")


def as_bool(value):
    """Interpret the boolean spellings accepted by historical configs."""

    if isinstance(value, bool):
        return value
    return str(value).lower() in ("yes", "true", "t", "1")


def resolve_filter_method(config):
    """Return the configured confidence-filtering method."""

    if "filter_method" in config:
        method = config["filter_method"]
    elif as_bool(config.get("skip_mixture_model_filter", False)):
        warnings.warn(
            "skip_mixture_model_filter is deprecated; use "
            "filter_method: robustness",
            DeprecationWarning,
            stacklevel=2,
        )
        method = "robustness"
    elif "filter_mode" in config:
        legacy_mode = config["filter_mode"]
        legacy_methods = {
            "unbiased": "robustness",
            "original": "legacy_mixture",
        }
        if legacy_mode not in legacy_methods:
            raise ValueError(
                "filter_mode must be 'unbiased' or 'original' "
                "(received {!r})".format(legacy_mode)
            )
        warnings.warn(
            "filter_mode is deprecated; use filter_method: {}".format(
                legacy_methods[legacy_mode]
            ),
            DeprecationWarning,
            stacklevel=2,
        )
        method = legacy_methods[legacy_mode]
    else:
        method = "robustness"

    if method not in FILTER_METHODS:
        raise ValueError(
            "filter_method must be one of: {} (received {!r})".format(
                ", ".join(sorted(FILTER_METHODS)), method
            )
        )
    return method


def selected_microexon_output(method):
    """Return the final microexon list produced by a filtering method."""

    try:
        return FILTER_OUTPUTS[method]
    except KeyError:
        raise ValueError("Unknown filter method {!r}".format(method))


def paired_read_samples(first_reads, paired_dict):
    """Expand paired-end sample identifiers to both mates for read evidence."""

    samples = []
    for first_read in first_reads:
        samples.extend((first_read, paired_dict[first_read]))
    return samples


def normalize_group_id(value, source, normalized_labels):
    """Create a path-safe group ID and reject ambiguous normalization."""

    stripped = value.strip()
    normalized = re.sub(r"\s+", "_", stripped)
    if not normalized or not GROUP_ID_PATTERN.match(normalized):
        raise ValueError(
            "{} group names must be non-empty and path-safe after replacing "
            "spaces with underscores (received {!r})".format(source, value)
        )
    previous = normalized_labels.get(normalized)
    if previous is not None and previous != stripped:
        raise ValueError(
            "{} group-name collision: {!r} and {!r} both normalize to {!r}".format(
                source, previous, stripped, normalized
            )
        )
    normalized_labels[normalized] = stripped
    return normalized


def load_filter_groups(config, samples, paired_dict=None):
    """Load the biological sample groups used by robustness filtering."""

    samples = set(samples)
    paired_dict = paired_dict or {}
    single_cell = as_bool(config.get("Single_Cell", False))
    if not single_cell and "bulk_samples" not in config:
        raise ValueError(
            "Bulk RNA-seq robustness filtering requires bulk_samples.tsv "
            "with sample and condition columns. Set bulk_samples in config.yaml."
        )
    if single_cell and "bulk_samples" in config and "cluster_metadata" not in config:
        raise ValueError(
            "A mixed bulk and single-cell run requires cluster_metadata so cell "
            "samples can be distinguished from bulk samples."
        )
    bulk_se = defaultdict(list)
    bulk_pe = defaultdict(list)
    manifest_samples = set()
    bulk_group_labels = {}
    if "bulk_samples" in config:
        with open(config["bulk_samples"]) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required_columns = {"sample", "condition"}
            if not required_columns.issubset(set(reader.fieldnames or [])):
                raise ValueError(
                    "bulk_samples.tsv must contain sample and condition columns"
                )
            for row in reader:
                sample = row["sample"]
                condition = row["condition"]
                if not sample.strip() or not condition.strip():
                    raise ValueError(
                        "bulk_samples.tsv contains an empty sample or condition"
                    )
                condition = normalize_group_id(
                    condition, "bulk_samples.tsv", bulk_group_labels
                )
                if sample in manifest_samples:
                    raise ValueError(
                        "bulk_samples.tsv contains duplicate sample {!r}".format(sample)
                    )
                manifest_samples.add(sample)
                if sample not in samples:
                    raise ValueError(
                        "bulk_samples.tsv contains unknown sample {!r}".format(sample)
                    )
                if sample in paired_dict:
                    bulk_pe[condition].append(sample)
                else:
                    bulk_se[condition].append(sample)

    single_cell_groups = {}
    single_cell_samples = set()
    if single_cell:
        if "cluster_metadata" in config:
            cell_column = config.get("file_basename", "sample")
            cluster_column = config.get("cluster_name", "cluster")
            grouped_cells = defaultdict(list)
            cluster_group_labels = {}
            with open(config["cluster_metadata"]) as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required_columns = {cell_column, cluster_column}
                if not required_columns.issubset(set(reader.fieldnames or [])):
                    raise ValueError(
                        "cluster_metadata must contain configured cell {!r} and "
                        "cluster {!r} columns".format(cell_column, cluster_column)
                    )
                for row in reader:
                    cell = row[cell_column]
                    cluster = row[cluster_column]
                    if not cell.strip() or not cluster.strip():
                        raise ValueError(
                            "cluster_metadata contains an empty cell or cluster"
                        )
                    cluster = normalize_group_id(
                        cluster, "cluster_metadata", cluster_group_labels
                    )
                    if cell not in samples:
                        raise ValueError(
                            "cluster_metadata contains unknown cell {!r}".format(cell)
                        )
                    if cell in single_cell_samples:
                        raise ValueError(
                            "cluster_metadata contains duplicate cell {!r}".format(cell)
                        )
                    grouped_cells[cluster].append(cell)
                    single_cell_samples.add(cell)
            single_cell_groups = {
                group: sorted(group_samples)
                for group, group_samples in grouped_cells.items()
            }
        else:
            single_cell_groups["all_cells"] = sorted(samples)
            single_cell_samples.update(samples)

    classified_samples = set(manifest_samples)
    classified_samples.update(single_cell_samples)
    for first_read in manifest_samples.intersection(set(paired_dict)):
        classified_samples.add(paired_dict[first_read])
    ungrouped_samples = samples.difference(classified_samples)
    if ungrouped_samples:
        raise ValueError(
            "Input samples missing from filtering metadata: {}".format(
                ", ".join(sorted(ungrouped_samples))
            )
        )

    return {
        "bulk_se": {
            group: sorted(group_samples)
            for group, group_samples in bulk_se.items()
        },
        "bulk_pe": {
            group: sorted(group_samples)
            for group, group_samples in bulk_pe.items()
        },
        "single_cell": single_cell_groups,
    }
