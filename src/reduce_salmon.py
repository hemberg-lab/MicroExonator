"""Reduce Salmon transcript quantifications into immutable gene shards."""

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
    genes = defaultdict(lambda: [0.0, 0.0, 0.0])
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
            values = genes[mapping[transcript]]
            values[0] += count
            values[1] += tpm
            values[2] += tpm * effective_length
    return {gene: (count, tpm, weighted / tpm if tpm else 0.0)
            for gene, (count, tpm, weighted) in genes.items()}


def reduce_group(quant_by_run, tx2gene_path, prefix, reference_id, manifest_sha256):
    if not quant_by_run:
        raise ValueError("Salmon shard requires at least one run")
    mapping = _tx2gene(tx2gene_path)
    runs = sorted(quant_by_run)
    values = {run: _quant(quant_by_run[run], mapping) for run in runs}
    genes = sorted(set().union(*(set(value) for value in values.values())))
    outputs = {}
    for kind, index in (("counts", 0), ("tpm", 1), ("length", 2)):
        lines = ["gene_id\t{}\n".format("\t".join(runs))]
        for gene in genes:
            lines.append("{}\t{}\n".format(gene, "\t".join(
                _number(values[run].get(gene, (0.0, 0.0, 0.0))[index]) for run in runs)))
        outputs[prefix + ".salmon_" + kind + ".tsv.gz"] = _gzip_text("".join(lines))
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
