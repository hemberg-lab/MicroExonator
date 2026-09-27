"""Reduce Salmon transcript quantifications into immutable group shards.

Gene level (for tximport/DESeq2): estimated counts, TPM and the TPM-weighted
effective length. A gene with zero TPM in a run gets the mean effective length
of its transcripts, as tximport does, so its offset is never log(0).
Transcript level (for SUPPA2): TPM and estimated counts.
"""

import argparse
import csv
import gzip
from collections import defaultdict

if __package__:
    from src.shard_guard import write_immutable_bundle
else:
    from shard_guard import write_immutable_bundle


def _number(value):
    return format(value, ".12g")


def _gzip_text(content):
    return gzip.compress(content.encode(), mtime=0)


def _tx2gene(path):
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not {"TXNAME", "GENEID"}.issubset(reader.fieldnames):
            raise ValueError("tx2gene requires TXNAME and GENEID columns")
        result = {}
        for row in reader:
            transcript, gene = row["TXNAME"], row["GENEID"]
            if transcript in result and result[transcript] != gene:
                raise ValueError("conflicting tx2gene assignment: {}".format(transcript))
            result[transcript] = gene
        return result


def _quant(path, mapping):
    genes = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0, 0])
    transcripts = {}
    with open(path, newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames or not {"Name", "EffectiveLength", "TPM", "NumReads"}.issubset(reader.fieldnames):
            raise ValueError("invalid Salmon quant.sf: {}".format(path))
        for row in reader:
            transcript = row["Name"]
            if transcript not in mapping:
                raise ValueError("unmapped transcript: {}".format(transcript))
            count = float(row["NumReads"])
            tpm = float(row["TPM"])
            effective_length = float(row["EffectiveLength"])
            if min(count, tpm, effective_length) < 0:
                raise ValueError("negative Salmon quantity for {}".format(transcript))
            transcripts[transcript] = (count, tpm)
            values = genes[mapping[transcript]]
            values[0] += count
            values[1] += tpm
            values[2] += tpm * effective_length
            values[3] += effective_length
            values[4] += 1
    gene_values = {gene: (count, tpm, weighted / tpm if tpm else length_sum / n)
                   for gene, (count, tpm, weighted, length_sum, n) in genes.items()}
    return gene_values, transcripts


def _matrix(label, rows, runs, lookup, index):
    lines = ["{}\t{}\n".format(label, "\t".join(runs))]
    for row in rows:
        lines.append("{}\t{}\n".format(row, "\t".join(
            _number(lookup[run].get(row, (0.0, 0.0, 0.0))[index]) for run in runs)))
    return _gzip_text("".join(lines))


def reduce_group(quant_by_run, tx2gene_path, prefix, reference_id, manifest_sha256):
    if not quant_by_run:
        raise ValueError("Salmon shard requires at least one run")
    mapping = _tx2gene(tx2gene_path)
    runs = sorted(quant_by_run)
    parsed = {run: _quant(quant_by_run[run], mapping) for run in runs}
    genes_by_run = {run: parsed[run][0] for run in runs}
    tx_by_run = {run: parsed[run][1] for run in runs}
    genes = sorted(set().union(*(set(value) for value in genes_by_run.values())))
    transcripts = sorted(set().union(*(set(value) for value in tx_by_run.values())))
    outputs = {}
    for kind, index in (("counts", 0), ("tpm", 1), ("length", 2)):
        outputs[prefix + ".salmon_" + kind + ".tsv.gz"] = _matrix(
            "gene_id", genes, runs, genes_by_run, index)
    for kind, index in (("tx_counts", 0), ("tx_tpm", 1)):
        outputs[prefix + ".salmon_" + kind + ".tsv.gz"] = _matrix(
            "transcript_id", transcripts, runs, tx_by_run, index)
    write_immutable_bundle(outputs, prefix + ".salmon_checksums.json", {
        "reference_id": reference_id, "manifest_sha256": manifest_sha256,
        "runs": runs,
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tx2gene", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--quant", action="append", required=True,
                        help="run_id=path/to/quant.sf")
    args = parser.parse_args()
    quant_by_run = dict(value.split("=", 1) for value in args.quant)
    if len(quant_by_run) != len(args.quant):
        parser.error("duplicate run ID")
    reduce_group(quant_by_run, args.tx2gene, args.prefix,
                 args.reference_id, args.manifest_sha256)


if __name__ == "__main__":
    main()
