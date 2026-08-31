"""Pure functions implementing MicroExonator robustness filtering."""

import csv
import gzip
from collections import defaultdict


ROBUST_OUTPUT_COLUMNS = (
    "ME",
    "Transcript",
    "Total_coverage",
    "Total_SJs",
    "ME_coverages",
    "ME_length",
    "ME_seq",
    "ME_matches",
    "U2_score",
    "Mean_conservation",
    "P_MEs",
    "Total_ME",
)


def detect_group_microexons(
    psi_files,
    spanning_read_files,
    min_ci_lower=0.1,
    min_detected_fraction=0.5,
    min_spanning_reads=3,
):
    """Return microexons reproducibly quantified within one sample group."""

    ci_lower_by_microexon = defaultdict(list)
    for file_name in psi_files:
        with gzip.open(file_name, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                ci_lower_by_microexon[row["ME"]].append(float(row["CI_Lo"]))

    spanning_reads = defaultdict(int)
    for file_name in spanning_read_files:
        with open(file_name) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                spanning_reads[row["ME"]] += int(row["Spanning_cov"])

    detected = {}
    for microexon, ci_lower_values in ci_lower_by_microexon.items():
        total_measurements = len(ci_lower_values)
        detected_samples = sum(
            value >= min_ci_lower for value in ci_lower_values
        )
        detected_fraction = detected_samples / float(total_measurements)
        if (
            detected_fraction > min_detected_fraction
            and spanning_reads[microexon] >= min_spanning_reads
        ):
            detected[microexon] = {
                "total_measurements": total_measurements,
                "detected_samples": detected_samples,
                "spanning_reads": spanning_reads[microexon],
            }
    return detected


def collect_robust_microexons(detection_files, min_detected_samples=1):
    """Collect microexons meeting the within-group sample threshold."""

    robust_microexons = set()
    for file_name in detection_files:
        with open(file_name) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if int(row["detected_samples"]) >= min_detected_samples:
                    robust_microexons.add(row["ME"])
    return robust_microexons


def write_robust_microexons(me_centric_file, output_file, robust_microexons):
    """Write selected ME-centric rows with their complete twelve-column schema."""

    with open(me_centric_file) as source, open(output_file, "w") as output:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(ROBUST_OUTPUT_COLUMNS)
        for row in reader:
            if row and row[0] in robust_microexons:
                writer.writerow(row)
