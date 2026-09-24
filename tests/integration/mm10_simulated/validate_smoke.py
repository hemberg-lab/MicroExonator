#!/usr/bin/env python3
"""Validate the committed mm10 fixtures and completed smoke-test outputs."""

import argparse
import csv
import gzip
import hashlib
import math
import pathlib
import sys
from collections import Counter
from collections import defaultdict


READ_LENGTH = 100
MINIMUM_VALID_SUPPORT = 12
# ME053 is a 3-nt CAG in a long intron that occurs again in the same intron.
# The current reporting stage collapses it into the other CAG coordinate, so it
# is deliberately the sole allowed absence from the final coordinate list.
EXPECTED_COLLAPSED_TRUTH_EVENTS = frozenset({"chr17_-_30598522_30598525"})
# Largest accepted gap between the mean corrected PSI of a group and the
# simulated target PSI. In the Python 3 baseline run, events with complete
# skipping evidence have a median gap of 0.027 and a 95th percentile of 0.10.
PSI_TOLERANCE = 0.15


class ValidationError(RuntimeError):
    """A fixture or pipeline output violates its documented contract."""


def _open_text(path):
    path = pathlib.Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def _fastq_records(path):
    with gzip.open(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip("\r\n")
            separator = handle.readline().rstrip("\r\n")
            qualities = handle.readline().rstrip("\r\n")
            if not sequence or not separator or not qualities:
                raise ValidationError("{} is truncated".format(path))
            if not header.startswith("@") or not separator.startswith("+"):
                raise ValidationError("{} is not valid four-line FASTQ".format(path))
            yield header[1:].strip().split()[0], sequence, qualities


def inspect_fastq(path, expected_length=READ_LENGTH):
    """Validate one gzip FASTQ and return basic read/quality statistics."""
    seen = set()
    read_count = 0
    minimum_quality = None
    maximum_quality = None
    problems = []
    for read_id, sequence, qualities in _fastq_records(path):
        read_count += 1
        if read_id in seen:
            problems.append("duplicate read ID {}".format(read_id))
        seen.add(read_id)
        if len(sequence) != expected_length:
            problems.append(
                "read {} is {} nt, not {} nt".format(
                    read_id, len(sequence), expected_length
                )
            )
        if len(sequence) != len(qualities):
            problems.append(
                "read {} sequence and quality lengths differ".format(read_id)
            )
        if set(sequence.upper()) - set("ACGT"):
            problems.append("read {} contains a non-ACGT base".format(read_id))
        numeric_qualities = [ord(character) - 33 for character in qualities]
        if numeric_qualities:
            minimum_quality = min(
                numeric_qualities
                if minimum_quality is None
                else numeric_qualities + [minimum_quality]
            )
            maximum_quality = max(
                numeric_qualities
                if maximum_quality is None
                else numeric_qualities + [maximum_quality]
            )
        if any(value < 0 or value > 41 for value in numeric_qualities):
            problems.append("read {} has a quality outside Phred+33 0-41".format(read_id))
    if read_count == 0:
        problems.append("FASTQ contains no reads")
    if problems:
        raise ValidationError("{}: {}".format(path, "; ".join(problems)))
    return {
        "reads": read_count,
        "unique_ids": len(seen),
        "minimum_quality": minimum_quality,
        "maximum_quality": maximum_quality,
    }


def _pair_stem(read_id):
    if read_id.endswith("/1") or read_id.endswith("/2"):
        return read_id[:-2]
    return read_id


def inspect_paired_fastqs(read1_path, read2_path, expected_length=READ_LENGTH):
    """Validate two mates, including order and identifier synchronization."""
    read1_stats = inspect_fastq(read1_path, expected_length)
    read2_stats = inspect_fastq(read2_path, expected_length)
    if read1_stats["reads"] != read2_stats["reads"]:
        raise ValidationError(
            "paired FASTQs are not synchronized: {} versus {} reads".format(
                read1_stats["reads"], read2_stats["reads"]
            )
        )
    for index, (read1, read2) in enumerate(
        zip(_fastq_records(read1_path), _fastq_records(read2_path)), start=1
    ):
        if _pair_stem(read1[0]) != _pair_stem(read2[0]):
            raise ValidationError(
                "paired FASTQs are not synchronized at record {}: {} versus {}".format(
                    index, read1[0], read2[0]
                )
            )
    return {"pairs": read1_stats["reads"], "read1": read1_stats, "read2": read2_stats}


def inspect_truth_support(
    truth_path,
    expected_events,
    expected_samples,
    minimum_valid_support=MINIMUM_VALID_SUPPORT,
):
    """Enforce per-event/sample included and skipped support in read truth."""
    support = Counter()
    truth_records = 0
    anchor_pairs = set()
    with _open_text(truth_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "read_id",
            "sample",
            "kind",
            "event_id",
            "ME",
            "isoform",
            "anchor_up",
            "anchor_down",
            "valid_evidence",
            "error_positions_r1",
            "error_positions_r2",
        }
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValidationError(
                "{} is missing truth columns: {}".format(
                    truth_path, ", ".join(sorted(missing))
                )
            )
        for row in reader:
            truth_records += 1
            if row["kind"] != "target":
                continue
            event = row["event_id"]
            sample = row["sample"]
            if event not in expected_events:
                raise ValidationError("truth contains unknown event {}".format(event))
            if sample not in expected_samples:
                raise ValidationError("truth contains unknown sample {}".format(sample))
            if row["isoform"] not in ("included", "skipped"):
                raise ValidationError(
                    "truth has invalid isoform {}".format(row["isoform"])
                )
            try:
                anchors = (int(row["anchor_up"]), int(row["anchor_down"]))
            except ValueError as error:
                raise ValidationError("truth anchors must be integers") from error
            anchor_pairs.add(anchors)
            if row["valid_evidence"] == "1":
                support[(event, sample, row["isoform"])] += 1
    problems = []
    for event in sorted(expected_events):
        for sample in sorted(expected_samples):
            for isoform in ("included", "skipped"):
                observed = support[(event, sample, isoform)]
                if observed < minimum_valid_support:
                    problems.append(
                        "{} {} {} has {} valid reads; minimum is {}".format(
                            event, sample, isoform, observed, minimum_valid_support
                        )
                    )
    if problems:
        raise ValidationError("; ".join(problems))
    return {
        "truth_records": truth_records,
        "valid_support_records": sum(support.values()),
        "distinct_anchor_pairs": len(anchor_pairs),
    }


def _read_event_table(events_path):
    events = {}
    with open(events_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "event_id" not in reader.fieldnames or "ME" not in reader.fieldnames:
            raise ValidationError("events.tsv must contain event_id and ME columns")
        for row in reader:
            if not row["event_id"] or not row["ME"]:
                raise ValidationError("events.tsv contains an empty event_id or ME")
            if row["event_id"] in events:
                raise ValidationError("duplicate event_id {}".format(row["event_id"]))
            events[row["event_id"]] = row["ME"]
    return events


def _read_bulk_samples(samples_path):
    samples = []
    with open(samples_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "condition"}.issubset(reader.fieldnames):
            raise ValidationError("bulk_samples.tsv needs sample and condition columns")
        for row in reader:
            samples.append(row["sample"])
    if len(samples) != len(set(samples)):
        raise ValidationError("bulk_samples.tsv contains duplicate samples")
    return samples


def _read_sample_groups(samples_path):
    with open(samples_path) as handle:
        return {
            row["sample"]: row["condition"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _read_target_psi(events_path):
    """Return {ME: {group: target PSI}} from psi_<group> columns, if present."""
    targets = {}
    with open(events_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        groups = [
            field[len("psi_"):]
            for field in reader.fieldnames or ()
            if field.startswith("psi_")
        ]
        for row in reader:
            targets[row["ME"]] = {
                group: float(row["psi_" + group]) for group in groups
            }
    return targets


def _parse_positions(value):
    if not value:
        return ()
    try:
        positions = tuple(int(position) for position in value.split(","))
    except ValueError as error:
        raise ValidationError("truth error positions must be comma-separated integers") from error
    if any(position < 1 or position > READ_LENGTH for position in positions):
        raise ValidationError("truth error position falls outside a 100 nt read")
    return positions


def _validate_events(events_path):
    expected_bins = {"1-6": 7, "7-12": 27, "13-18": 28, "19-24": 22, "25-27": 16}
    expected_profiles = {
        "stable_medium",
        "stable_high",
        "A_enriched",
        "B_enriched",
        "C_enriched",
    }
    with open(events_path) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 100:
        raise ValidationError("events.tsv has {} events; expected 100".format(len(rows)))
    if len({row["event_id"] for row in rows}) != 100:
        raise ValidationError("events.tsv event IDs are not unique")
    if len({row["ME"] for row in rows}) != 100:
        raise ValidationError("events.tsv microexon coordinates are not unique")
    if len({row["gene_id"] for row in rows}) != 100:
        raise ValidationError("events.tsv contains more than one event per gene")
    for row in rows:
        try:
            length = int(row["length"])
            start = int(row["start"])
            end = int(row["end"])
        except ValueError as error:
            raise ValidationError("events.tsv coordinates and lengths must be integers") from error
        if end - start != length or not 1 <= length <= 27:
            raise ValidationError("{} has inconsistent coordinates or length".format(row["event_id"]))
        if row["motif"] != "AGGT":
            raise ValidationError("{} does not have canonical AG-GT boundaries".format(row["event_id"]))
    status_counts = Counter(row["annotation_status"] for row in rows)
    if not status_counts["annotated"] or not status_counts["non_annotated"]:
        raise ValidationError("events.tsv must contain both annotated and non-annotated events")
    bins = Counter(row["length_bin"] for row in rows)
    if bins != Counter(expected_bins):
        raise ValidationError("combined length-bin counts are incorrect: {}".format(dict(bins)))
    profiles = Counter(row["profile"] for row in rows)
    if profiles != Counter({profile: 20 for profile in expected_profiles}):
        raise ValidationError("combined profile counts are incorrect: {}".format(dict(profiles)))
    return rows


def _validate_checksums(fixture_root):
    checksum_path = pathlib.Path(fixture_root) / "SHA256SUMS"
    if not checksum_path.is_file():
        raise ValidationError("fixture is missing SHA256SUMS")
    checked = 0
    with open(checksum_path) as handle:
        for line in handle:
            if not line.strip():
                continue
            expected, relative = line.rstrip("\n").split("  ", 1)
            artifact = pathlib.Path(fixture_root) / relative
            if not artifact.is_file():
                raise ValidationError("checksum target is missing: {}".format(relative))
            digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
            if digest != expected:
                raise ValidationError("checksum mismatch for {}".format(relative))
            checked += 1
    if checked == 0:
        raise ValidationError("SHA256SUMS contains no artifacts")
    return checked


def _read_truth_details(truth_path, events_by_id, samples):
    rows = []
    seen = set()
    target_count = 0
    background_count = 0
    weak_count = 0
    anchor_pairs = set()
    errors = {1: {}, 2: {}}
    inserts = []
    informative_mates = set()
    orientations = set()
    with gzip.open(truth_path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            read_id = row["read_id"]
            if read_id in seen:
                raise ValidationError("truth contains duplicate read ID {}".format(read_id))
            seen.add(read_id)
            if row["sample"] not in samples:
                raise ValidationError("truth contains unknown sample {}".format(row["sample"]))
            errors[1][read_id] = _parse_positions(row["error_positions_r1"])
            errors[2][read_id] = _parse_positions(row["error_positions_r2"])
            orientations.add(row["orientation"])
            informative_mates.add(row["informative_mate"])
            if row["insert_size"]:
                insert_size = int(row["insert_size"])
                if not 140 <= insert_size <= 260:
                    raise ValidationError("paired insert size is outside 140-260 nt")
                inserts.append(insert_size)
            if row["kind"] == "background":
                background_count += 1
            elif row["kind"] == "target":
                target_count += 1
                if row["event_id"] not in events_by_id:
                    raise ValidationError("truth contains unknown event {}".format(row["event_id"]))
                if row["ME"] != events_by_id[row["event_id"]]["ME"]:
                    raise ValidationError("truth event {} has the wrong coordinate".format(row["event_id"]))
                anchors = (int(row["anchor_up"]), int(row["anchor_down"]))
                anchor_pairs.add(anchors)
                weak = row["weak_anchor"] == "1"
                valid = row["valid_evidence"] == "1"
                if weak:
                    weak_count += 1
                    if valid or min(anchors) not in range(3, 8):
                        raise ValidationError("weak-anchor truth row is marked as valid or lacks a 3-7 nt anchor")
                elif not valid or min(anchors) < 8:
                    raise ValidationError("strong-anchor truth row lacks valid >=8 nt anchors")
            else:
                raise ValidationError("truth kind must be target or background")
            rows.append(row)
    background_fraction = background_count / float(target_count)
    weak_fraction = weak_count / float(target_count)
    if not 0.09 <= background_fraction <= 0.11:
        raise ValidationError("background reads are not approximately 10% of target reads")
    if not 0.13 <= weak_fraction <= 0.17:
        raise ValidationError("weak-anchor reads are not approximately 15% of target reads")
    return {
        "rows": rows,
        "read_ids": seen,
        "errors": errors,
        "target_reads": target_count,
        "background_reads": background_count,
        "distinct_anchor_pairs": len(anchor_pairs),
        "inserts": inserts,
        "informative_mates": informative_mates,
        "orientations": orientations,
    }


def _quality_error_statistics(fastq_paths, error_positions_by_mate, paired=False):
    observed_ids = set()
    quality_sum = 0
    quality_count = 0
    error_qualities = []
    for mate, paths in fastq_paths.items():
        for path in paths:
            for read_id, _, qualities in _fastq_records(path):
                stem = _pair_stem(read_id) if paired else read_id
                if stem in observed_ids and not paired:
                    raise ValidationError("read ID {} occurs in multiple FASTQs".format(stem))
                observed_ids.add(stem)
                numeric = [ord(character) - 33 for character in qualities]
                quality_sum += sum(numeric)
                quality_count += len(numeric)
                for position in error_positions_by_mate[mate].get(stem, ()):
                    error_qualities.append(numeric[position - 1])
    if not error_qualities:
        raise ValidationError("fixture contains no simulated substitution errors")
    overall_mean = quality_sum / float(quality_count)
    error_mean = sum(error_qualities) / float(len(error_qualities))
    if error_mean >= overall_mean:
        raise ValidationError("substitution errors are not associated with lower-quality positions")
    return observed_ids, overall_mean, error_mean, len(error_qualities)


def validate_fixture(fixture_root, layout, verify_checksums=True):
    """Validate all committed artifacts for one fixture layout."""
    if layout not in ("single_end", "paired_end"):
        raise ValidationError("layout must be single_end or paired_end")
    fixture_root = pathlib.Path(fixture_root)
    events = _validate_events(fixture_root / "truth" / "events.tsv")
    events_by_id = {row["event_id"]: row for row in events}
    layout_root = fixture_root / layout
    samples = _read_bulk_samples(layout_root / "bulk_samples.tsv")
    if len(samples) != 9:
        raise ValidationError("{} bulk manifest has {} samples; expected 9".format(layout, len(samples)))
    truth_path = fixture_root / "truth" / "{}_reads.tsv.gz".format(layout)
    truth = _read_truth_details(truth_path, events_by_id, set(samples))
    support = inspect_truth_support(truth_path, set(events_by_id), set(samples))

    fastq_root = layout_root / "fastq"
    fastq_reads = 0
    fastq_pairs = 0
    fastq_paths = defaultdict(list)
    if layout == "single_end":
        for sample in samples:
            path = fastq_root / "{}.fastq.gz".format(sample)
            stats = inspect_fastq(path)
            fastq_reads += stats["reads"]
            fastq_paths[1].append(path)
    else:
        pair_path = layout_root / "paired_samples.tsv"
        pairs = []
        with open(pair_path) as handle:
            for row in csv.reader(handle, delimiter="\t"):
                if len(row) != 2:
                    raise ValidationError("paired_samples.tsv must be headerless with two columns")
                pairs.append(tuple(row))
        if len(pairs) != 9 or [pair[0] for pair in pairs] != samples:
            raise ValidationError("paired_samples.tsv does not match the nine biological samples")
        for read1, read2 in pairs:
            read1_path = fastq_root / "{}.fastq.gz".format(read1)
            read2_path = fastq_root / "{}.fastq.gz".format(read2)
            stats = inspect_paired_fastqs(read1_path, read2_path)
            fastq_pairs += stats["pairs"]
            fastq_paths[1].append(read1_path)
            fastq_paths[2].append(read2_path)

    observed_ids, overall_q, error_q, substitutions = _quality_error_statistics(
        fastq_paths, truth["errors"], paired=layout == "paired_end"
    )
    if observed_ids != truth["read_ids"]:
        raise ValidationError("FASTQ and read-truth identifiers do not match")
    if layout == "single_end" and fastq_reads != len(truth["rows"]):
        raise ValidationError("single-end FASTQ and truth record counts differ")
    if layout == "paired_end" and fastq_pairs != len(truth["rows"]):
        raise ValidationError("paired FASTQ and truth record counts differ")
    if truth["orientations"] != {"+", "-"}:
        raise ValidationError("both transcript orientations must occur")
    if layout == "paired_end":
        if not {"R1", "R2", "both"}.issubset(truth["informative_mates"]):
            raise ValidationError("paired truth must include R1, R2, and overlapping evidence")
        if not 180 <= sum(truth["inserts"]) / float(len(truth["inserts"])) <= 205:
            raise ValidationError("paired insert-size mean is outside the expected range")

    checksums = _validate_checksums(fixture_root) if verify_checksums else 0
    return {
        "events": len(events),
        "samples": len(samples),
        "truth_records": len(truth["rows"]),
        "fastq_reads": fastq_reads,
        "fastq_pairs": fastq_pairs,
        "distinct_anchor_pairs": support["distinct_anchor_pairs"],
        "valid_support_records": support["valid_support_records"],
        "substitutions": substitutions,
        "overall_mean_quality": overall_q,
        "error_mean_quality": error_q,
        "checksums": checksums,
    }


def _numeric(row, field, sample, microexon):
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as error:
        raise ValidationError(
            "sample {} event {} field {} must be numeric".format(
                sample, microexon, field
            )
        ) from error
    if not math.isfinite(value):
        raise ValidationError(
            "sample {} event {} field {} must be numeric and finite".format(
                sample, microexon, field
            )
        )
    return value


def validate_pipeline_outputs(
    layout,
    events_path,
    samples_path,
    run_directory,
    expected_missing=(),
    psi_tolerance=PSI_TOLERANCE,
):
    """Validate discovery and quantification except declared collapsed events.

    When events.tsv carries psi_<group> target columns and psi_tolerance is not
    None, the mean corrected PSI of every truth event in every sample group must
    lie within psi_tolerance of its simulated target.
    """
    if layout not in ("single_end", "paired_end"):
        raise ValidationError("layout must be single_end or paired_end")
    run_directory = pathlib.Path(run_directory)
    events = _read_event_table(events_path)
    truth_microexons = set(events.values())
    samples = _read_bulk_samples(samples_path)
    expected_missing = set(expected_missing)
    unknown_expected = expected_missing - truth_microexons
    if unknown_expected:
        raise ValidationError(
            "expected missing events are absent from truth: {}".format(
                ", ".join(sorted(unknown_expected))
            )
        )

    robust_path = run_directory / "Report" / "out.robustly_detected.txt"
    with open(robust_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        discovered = {row["ME"] for row in reader}
    missing_discoveries = truth_microexons - discovered
    if missing_discoveries != expected_missing:
        raise ValidationError(
            "robust output missing {} truth events (expected only {}): {}".format(
                len(missing_discoveries),
                ", ".join(sorted(expected_missing)) or "none",
                ", ".join(sorted(missing_discoveries)[:10]) or "none",
            )
        )

    layout_code = "se" if layout == "single_end" else "pe"
    required_microexons = truth_microexons - expected_missing
    quantified = 0
    psi_deviations = []
    sample_groups = _read_sample_groups(samples_path)
    observed_psi = defaultdict(list)
    for sample in samples:
        quant_path = (
            run_directory
            / "Report"
            / "quant"
            / "corrected"
            / "PSI_sparse"
            / "bulk"
            / layout_code
            / "{}.corrected.PSI.gz".format(sample)
        )
        rows = {}
        with gzip.open(quant_path, "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                rows[row["ME"]] = row
            missing_quantifications = sorted(required_microexons - set(rows))
        if missing_quantifications:
            raise ValidationError(
                "sample {} is missing {} truth quantifications: {}".format(
                    sample,
                    len(missing_quantifications),
                    ", ".join(missing_quantifications[:10]),
                )
            )
        for microexon in required_microexons:
            row = rows[microexon]
            psi = _numeric(row, "PSI", sample, microexon)
            ci_low = _numeric(row, "CI_Lo", sample, microexon)
            ci_high = _numeric(row, "CI_Hi", sample, microexon)
            if not (0.0 <= ci_low <= psi <= ci_high <= 1.0):
                raise ValidationError(
                    "sample {} event {} has invalid PSI/confidence bounds".format(
                        sample, microexon
                    )
                )
            quantified += 1
            observed_psi[(microexon, sample_groups[sample])].append(psi)

    targets = _read_target_psi(events_path)
    for (microexon, group), values in sorted(observed_psi.items()):
        target = targets.get(microexon, {}).get(group)
        if target is None:
            continue
        observed = sum(values) / len(values)
        psi_deviations.append(
            {
                "ME": microexon,
                "group": group,
                "observed": observed,
                "target": target,
                "deviation": observed - target,
            }
        )
    if psi_tolerance is not None:
        failures = [
            item for item in psi_deviations
            if abs(item["deviation"]) > psi_tolerance
        ]
        if failures:
            failures.sort(key=lambda item: -abs(item["deviation"]))
            failing_events = sorted({item["ME"] for item in failures})
            raise ValidationError(
                "{} truth events deviate from target PSI by more than {} "
                "in at least one group; largest: {}".format(
                    len(failing_events),
                    psi_tolerance,
                    ", ".join(
                        "{ME} {group} observed {observed:.2f} target {target:.2f}".format(**item)
                        for item in failures[:10]
                    ),
                )
            )

    return {
        "truth_events": len(truth_microexons),
        "reported_truth_events": len(truth_microexons - missing_discoveries),
        "expected_missing_events": sorted(expected_missing),
        "samples": len(samples),
        "quantified_event_samples": quantified,
        "additional_discoveries": sorted(discovered - truth_microexons),
        "psi_deviations": psi_deviations,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout", required=True, choices=("single_end", "paired_end"))
    parser.add_argument("--fixture-dir", type=pathlib.Path, default=pathlib.Path(__file__).parent)
    parser.add_argument("--run-dir", type=pathlib.Path, required=True)
    parser.add_argument(
        "--psi-tolerance",
        type=float,
        default=PSI_TOLERANCE,
        help="largest accepted |mean group PSI - target PSI| (default: %(default)s)",
    )
    parser.add_argument(
        "--skip-psi-accuracy",
        action="store_true",
        help="check PSI bounds only, not agreement with the simulated targets",
    )
    args = parser.parse_args(argv)
    report = validate_pipeline_outputs(
        args.layout,
        args.fixture_dir / "truth" / "events.tsv",
        args.fixture_dir / args.layout / "bulk_samples.tsv",
        args.run_dir,
        expected_missing=EXPECTED_COLLAPSED_TRUTH_EVENTS,
        psi_tolerance=None if args.skip_psi_accuracy else args.psi_tolerance,
    )
    print(
        "Validated {reported_truth_events}/{truth_events} reported truth events across {samples} samples "
        "({quantified_event_samples} event-sample quantifications).".format(**report)
    )
    print("Additional discoveries (non-failing): {}".format(len(report["additional_discoveries"])))
    if report["psi_deviations"]:
        largest = max(abs(item["deviation"]) for item in report["psi_deviations"])
        print("Largest group PSI deviation from target: {:.3f}".format(largest))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValidationError) as error:
        print("ERROR: {}".format(error), file=sys.stderr)
        sys.exit(1)
