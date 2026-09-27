"""One frozen reference bundle shared by all umbrella quantifiers."""

import gzip
import json
import re
from pathlib import Path

from src.shard_guard import (hisat2_members, reference_manifest,
                             validate_fixed_microexons, write_immutable)


UMBRELLA_REFERENCE = config.get("umbrella_reference")
if not isinstance(UMBRELLA_REFERENCE, dict):
    raise WorkflowError("umbrella_reference mapping is required with umbrella_manifest")
UMBRELLA_REFERENCE_IDS = {run.reference_id for run in UMBRELLA_MANIFEST.included_runs()}
if len(UMBRELLA_REFERENCE_IDS) != 1:
    raise WorkflowError("umbrella_manifest must use one reference_id per invocation")
UMBRELLA_REFERENCE_ID = next(iter(UMBRELLA_REFERENCE_IDS))
UMBRELLA_REFERENCE_ROOT = "umbrella/reference/" + UMBRELLA_REFERENCE_ID
for field in ("genome_fasta", "annotation_gtf", "whippet_gtf", "me_db",
              "transcriptome_fasta", "decoys", "salmon_gtf", "splice_sites"):
    if not UMBRELLA_REFERENCE.get(field):
        raise WorkflowError("umbrella_reference requires {}".format(field))

UMBRELLA_HISAT_LARGE = UMBRELLA_REFERENCE.get("hisat2_index_type", "small") == "large"
if UMBRELLA_REFERENCE.get("hisat2_index_type", "small") not in ("small", "large"):
    raise WorkflowError("hisat2_index_type must be small or large")
UMBRELLA_HISAT_PREFIX = UMBRELLA_REFERENCE.get(
    "hisat2_index_prefix", UMBRELLA_REFERENCE_ROOT + "/hisat2/genome")
UMBRELLA_HISAT_MEMBERS = hisat2_members(UMBRELLA_HISAT_PREFIX, UMBRELLA_HISAT_LARGE)
UMBRELLA_WHIPPET_INDEX = UMBRELLA_REFERENCE.get(
    "whippet_index", UMBRELLA_REFERENCE_ROOT + "/whippet/whippet.jls")
UMBRELLA_WHIPPET_MEMBERS = [UMBRELLA_WHIPPET_INDEX,
                           UMBRELLA_WHIPPET_INDEX + ".exons.tab.gz"]
UMBRELLA_SALMON_INDEX = UMBRELLA_REFERENCE.get(
    "salmon_index", UMBRELLA_REFERENCE_ROOT + "/salmon/index")
UMBRELLA_TX2GENE = UMBRELLA_REFERENCE_ROOT + "/tx2gene.tsv"
UMBRELLA_FIXED_UNIVERSE = UMBRELLA_REFERENCE_ROOT + "/fixed_universe.json"
UMBRELLA_REFERENCE_MANIFEST = UMBRELLA_REFERENCE_ROOT + "/manifest.json"


rule umbrella_fixed_universe:
    input:
        fixed_gtf=UMBRELLA_REFERENCE["whippet_gtf"],
        me_db=UMBRELLA_REFERENCE["me_db"]
    output:
        protected(UMBRELLA_FIXED_UNIVERSE)
    run:
        count = validate_fixed_microexons(input.fixed_gtf, input.me_db)
        write_immutable(output[0], (json.dumps({"microexons": count}) + "\n").encode())


if "hisat2_index_prefix" not in UMBRELLA_REFERENCE:
    rule umbrella_hisat2_index:
        input:
            UMBRELLA_REFERENCE["genome_fasta"]
        output:
            UMBRELLA_HISAT_MEMBERS
        params:
            prefix=UMBRELLA_HISAT_PREFIX,
            large="--large-index" if UMBRELLA_HISAT_LARGE else ""
        threads: 8
        conda:
            "../envs/umbrella-quant.yaml"
        shell:
            "hisat2-build {params.large} -p {threads} {input} {params.prefix}"


if "whippet_index" not in UMBRELLA_REFERENCE:
    rule umbrella_whippet_index:
        input:
            genome=UMBRELLA_REFERENCE["genome_fasta"],
            fixed_gtf=UMBRELLA_REFERENCE["whippet_gtf"],
            me_db=UMBRELLA_REFERENCE["me_db"],
            validated=UMBRELLA_FIXED_UNIVERSE
        output:
            UMBRELLA_WHIPPET_MEMBERS
        params:
            julia=config.get("julia", "julia"),
            whippet_bin=config.get("whippet_bin_folder", "")
        log:
            UMBRELLA_REFERENCE_ROOT + "/whippet_index.log"
        shell:
            "{params.julia} {params.whippet_bin}/whippet-index.jl "
            "--fasta {input.genome} --gtf {input.fixed_gtf} "
            "--index {output[0]} 2> {log}"


rule umbrella_tx2gene:
    input:
        UMBRELLA_REFERENCE["salmon_gtf"]
    output:
        protected(UMBRELLA_TX2GENE)
    run:
        opener = gzip.open if str(input[0]).endswith(".gz") else open
        mapping = {}
        with opener(input[0], "rt") as stream:
            for line in stream:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9 or fields[2] != "transcript":
                    continue
                transcript = re.search(r'(?:^|;\s*)transcript_id "([^"]+)"', fields[8])
                gene = re.search(r'(?:^|;\s*)gene_id "([^"]+)"', fields[8])
                if not transcript or not gene:
                    raise WorkflowError("Salmon GTF transcript lacks gene_id/transcript_id")
                transcript, gene = transcript.group(1), gene.group(1)
                if transcript in mapping and mapping[transcript] != gene:
                    raise WorkflowError("conflicting gene for transcript {}".format(transcript))
                mapping[transcript] = gene
        if not mapping:
            raise WorkflowError("Salmon GTF contains no transcript features")
        content = "TXNAME\tGENEID\n" + "".join(
            "{}\t{}\n".format(tx, mapping[tx]) for tx in sorted(mapping))
        write_immutable(output[0], content.encode())


if "salmon_index" not in UMBRELLA_REFERENCE:
    rule umbrella_gentrome:
        input:
            transcripts=UMBRELLA_REFERENCE["transcriptome_fasta"],
            genome=UMBRELLA_REFERENCE["genome_fasta"]
        output:
            temp(UMBRELLA_REFERENCE_ROOT + "/salmon/gentrome.fa")
        run:
            Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
            with open(output[0], "wb") as destination:
                for source in (input.transcripts, input.genome):
                    opener = gzip.open if str(source).endswith(".gz") else open
                    with opener(source, "rb") as stream:
                        last = b"\n"
                        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                            destination.write(chunk)
                            last = chunk[-1:]
                        if last != b"\n":
                            destination.write(b"\n")

    rule umbrella_salmon_index:
        input:
            gentrome=UMBRELLA_REFERENCE_ROOT + "/salmon/gentrome.fa",
            decoys=UMBRELLA_REFERENCE["decoys"]
        output:
            directory(UMBRELLA_SALMON_INDEX)
        threads: 8
        conda:
            "../envs/umbrella-quant.yaml"
        shell:
            "salmon index -t {input.gentrome} -d {input.decoys} "
            "-i {output} -p {threads}"


rule umbrella_reference_manifest:
    input:
        genome=UMBRELLA_REFERENCE["genome_fasta"],
        annotation=UMBRELLA_REFERENCE["annotation_gtf"],
        whippet_annotation=UMBRELLA_REFERENCE["whippet_gtf"],
        me_db=UMBRELLA_REFERENCE["me_db"],
        salmon_gtf=UMBRELLA_REFERENCE["salmon_gtf"],
        transcripts=UMBRELLA_REFERENCE["transcriptome_fasta"],
        decoys=UMBRELLA_REFERENCE["decoys"],
        splice_sites=UMBRELLA_REFERENCE["splice_sites"],
        annotation_bed12=[UMBRELLA_REFERENCE["annotation_bed12"]] if UMBRELLA_REFERENCE.get("annotation_bed12") else [],
        tx2gene=UMBRELLA_TX2GENE,
        fixed_universe=UMBRELLA_FIXED_UNIVERSE,
        hisat2=UMBRELLA_HISAT_MEMBERS,
        whippet=UMBRELLA_WHIPPET_MEMBERS,
        salmon=UMBRELLA_SALMON_INDEX
    output:
        protected(UMBRELLA_REFERENCE_MANIFEST)
    run:
        paths = {"genome": str(input.genome), "annotation": str(input.annotation),
                 "whippet_annotation": str(input.whippet_annotation),
                 "me_db": str(input.me_db), "salmon_gtf": str(input.salmon_gtf),
                 "transcripts": str(input.transcripts), "decoys": str(input.decoys),
                 "tx2gene": str(input.tx2gene), "hisat2": list(input.hisat2),
                 "whippet": list(input.whippet), "salmon": str(input.salmon),
                 "fixed_universe": str(input.fixed_universe)}
        paths["splice_sites"] = str(input.splice_sites)
        if input.annotation_bed12:
            paths["annotation_bed12"] = str(input.annotation_bed12[0])
        manifest = reference_manifest(paths, expected_id=UMBRELLA_REFERENCE_ID,
                                      versions=UMBRELLA_REFERENCE.get("versions", {}))
        write_immutable(output[0], (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode())
