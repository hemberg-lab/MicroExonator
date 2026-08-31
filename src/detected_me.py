"""Snakemake adapter for within-group robustness filtering."""

import csv

from src.robustness_filter import detect_group_microexons


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
