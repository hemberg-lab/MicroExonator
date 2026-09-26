"""Differential microexon inclusion between two groups of samples, without Whippet.

A Python port of the model in Whippet's whippet-delta.jl (src/diff.jl),
applied directly to MicroExonator's per-sample PSI and read counts:

- each sample with a PSI and at least min_reads reads gives the posterior
  Beta(PSI * N + 1, (1 - PSI) * N + 1), from which `size` values are drawn;
- the draws of a group's replicates are pooled, a Beta is fitted to them by
  the method of moments (Distributions.jl fit(Beta, x)), and `size` new values
  are drawn from the fitted Beta;
- Psi_A and Psi_B are the means of the fitted Betas, DeltaPsi = Psi_A - Psi_B,
  and Probability = max(P(A - B > 0), P(B - A > 0)) over the paired draws.

N is the sample's corrected inclusion plus exclusion reads. The Whippet route
used Whippet's own Total_Reads for the node instead.

Each event gets its own random stream (seed plus a hash of the microexon), so
its result does not depend on which other events are in the run.
"""

import argparse
import csv
import gzip
import re
import sys
import zlib

import numpy


HEADER = ["exon_ID", "Gene", "Node", "Coord", "Strand", "Type",
          "Psi_A", "Psi_B", "DeltaPsi", "Probability", "Complexity", "Entropy",
          "Samples_A", "Samples_B", "Reads_A", "Reads_B"]


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_sample(path):
    """Return {ME: (PSI, N)} for the microexons with a PSI in one sample."""

    sample = dict()
    with open_text(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["PSI"] in ("NA", "nan", ""):
                continue
            if "excluding_covs" in row:
                ME = row["ME"]
                N = float(row["ME_coverages"]) + float(row["excluding_covs"])
            else:
                # uncorrected table: one coverage per junction
                ME = row["ME_coords"]
                N = sum(map(float, row["ME_coverages"].split(","))) + sum(map(float, row["SJ_coverages"].split(",")))
            sample[ME] = (float(row["PSI"]), N)
    return sample


def posterior_draws(psi, N, rng, size):
    return rng.beta(psi * N + 1.0, (1.0 - psi) * N + 1.0, size)


def fit_beta(x):
    """Method-of-moments Beta fit, as Distributions.jl fit(Beta, x)."""

    mean = x.mean()
    var = x.var(ddof=1)
    temp = mean * (1.0 - mean) / var - 1.0
    return mean * temp, (1.0 - mean) * temp


def group_posterior(draws, rng, size):
    alpha, beta = fit_beta(numpy.concatenate(draws))
    return alpha / (alpha + beta), rng.beta(alpha, beta, size)


def probability(a, b):
    return max(numpy.mean(a - b > 0), numpy.mean(b - a > 0))


def compare(values_A, values_B, rng, min_reads=5, min_samples=1, size=1000):
    """values_X: list of (PSI, N), one per sample with a PSI.

    Returns None when a group has fewer than min_samples usable samples.
    """

    used_A = [(psi, N) for psi, N in values_A if psi >= 0 and N >= min_reads]
    used_B = [(psi, N) for psi, N in values_B if psi >= 0 and N >= min_reads]
    if len(used_A) < min_samples or len(used_B) < min_samples:
        return None

    # same draw order as whippet-delta: A replicates, B replicates, then the fitted A and B
    draws_A = [posterior_draws(psi, N, rng, size) for psi, N in used_A]
    draws_B = [posterior_draws(psi, N, rng, size) for psi, N in used_B]
    psi_A, fitted_A = group_posterior(draws_A, rng, size)
    psi_B, fitted_B = group_posterior(draws_B, rng, size)

    return {"Psi_A": psi_A, "Psi_B": psi_B, "DeltaPsi": psi_A - psi_B,
            "Probability": probability(fitted_A, fitted_B),
            "Samples_A": len(used_A), "Samples_B": len(used_B),
            "Reads_A": sum(N for psi, N in used_A), "Reads_B": sum(N for psi, N in used_B)}


def event_rng(seed, ME):
    return numpy.random.default_rng([seed, zlib.crc32(ME.encode())])


def transcript_genes(gtf):
    genes = dict()
    with open_text(gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            match = re.search('gene_id "([^"]+)".*transcript_id "([^"]+)"', line)
            if match:
                genes[match.group(2)] = match.group(1)
                genes[match.group(2).split(".")[0]] = match.group(1)
    return genes


def read_microexons(path, genes):
    """Microexons to report, with a gene ID from their transcript when known."""

    microexons = dict()
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            transcript = row.get("Transcript", "")
            microexons[row["ME"]] = genes.get(transcript, genes.get(transcript.split(".")[0], "NA"))
    return microexons


def fmt(x):
    return "%.5g" % x


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-a", required=True, help="comma-separated PSI tables for group A")
    parser.add_argument("-b", required=True, help="comma-separated PSI tables for group B")
    parser.add_argument("--microexons", required=True, help="microexons to report (out.robustly_detected.txt or out.high_quality.txt)")
    parser.add_argument("--gtf", help="annotation GTF, to add gene IDs")
    parser.add_argument("--min-reads", type=float, default=5)
    parser.add_argument("--min-samples", type=int, default=1)
    parser.add_argument("--size", type=int, default=1000, help="empirical distribution size")
    parser.add_argument("--seed", type=int, default=123456)
    args = parser.parse_args(argv)

    samples_A = [read_sample(path) for path in args.a.split(",")]
    samples_B = [read_sample(path) for path in args.b.split(",")]
    genes = transcript_genes(args.gtf) if args.gtf else dict()
    microexons = read_microexons(args.microexons, genes)

    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(HEADER)
    for ME in sorted(microexons):
        result = compare([s[ME] for s in samples_A if ME in s], [s[ME] for s in samples_B if ME in s],
                         event_rng(args.seed, ME), args.min_reads, args.min_samples, args.size)
        if result is None:
            continue
        chrom = "_".join(ME.split("_")[:-3])
        strand, start, end = ME.split("_")[-3:]
        out.writerow([ME, microexons[ME], "NA", "%s:%d-%s" % (chrom, int(start) + 1, end), strand, "CE",
                      fmt(result["Psi_A"]), fmt(result["Psi_B"]), fmt(result["DeltaPsi"]),
                      "%.3f" % result["Probability"], "NA", "NA",
                      result["Samples_A"], result["Samples_B"], fmt(result["Reads_A"]), fmt(result["Reads_B"])])


if __name__ == "__main__":
    main()
