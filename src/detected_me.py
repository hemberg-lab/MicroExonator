"""Snakemake adapter for within-group robustness filtering."""

import csv
import os
import sys

# Snakemake executes ``script:`` rules from a temporary copy under
# ``.snakemake/scripts``.  Resolve the project's staged source directory from
# the workflow working directory instead of relying on that temporary copy's
# import path.
sys.path.insert(0, os.path.join(os.getcwd(), "src"))

from robustness_filter import detect_group_microexons


detected = detect_group_microexons(
    snakemake.input["PSI_files"],
    snakemake.input["ME_reads"],
)

with open(snakemake.output["detected"], "w") as output:
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(
        ("ME", "total_measurements", "detected_samples", "spanning_reads")
    )
    for microexon in sorted(detected):
        evidence = detected[microexon]
        writer.writerow(
            (
                microexon,
                evidence["total_measurements"],
                evidence["detected_samples"],
                evidence["spanning_reads"],
            )
        )
