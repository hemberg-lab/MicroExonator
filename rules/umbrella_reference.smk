"""One frozen reference bundle shared by all umbrella quantifiers."""

import gzip
import json
import re
from pathlib import Path

from src.shard_guard import hisat2_members, reference_identity, reference_manifest, write_immutable


UMBRELLA_REFERENCE = config.get("umbrella_reference")
if not isinstance(UMBRELLA_REFERENCE, dict):
    raise WorkflowError("umbrella_reference mapping is required with umbrella_manifest")
UMBRELLA_REFERENCE_IDS = {run.reference_id for run in UMBRELLA_MANIFEST.included_runs()}
if len(UMBRELLA_REFERENCE_IDS) != 1:
    raise WorkflowError("umbrella_manifest must use one reference_id per invocation")
UMBRELLA_REFERENCE_ID = next(iter(UMBRELLA_REFERENCE_IDS))
UMBRELLA_REFERENCE_ROOT = "umbrella/reference/" + UMBRELLA_REFERENCE_ID
for field in ("genome_fasta", "annotation_gtf", "whippet_gtf", "me_db"):
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
UMBRELLA_FIXED_MICROEXONS = UMBRELLA_REFERENCE_ROOT + "/fixed_microexons.tsv"
UMBRELLA_REFERENCE_MANIFEST = UMBRELLA_REFERENCE_ROOT + "/manifest.json"
UMBRELLA_SALMON_GTF = UMBRELLA_REFERENCE.get("salmon_gtf", UMBRELLA_REFERENCE["annotation_gtf"])
UMBRELLA_TRANSCRIPTS = UMBRELLA_REFERENCE.get(
    "transcriptome_fasta", UMBRELLA_REFERENCE_ROOT + "/salmon/transcripts.fa")
UMBRELLA_DECOYS = UMBRELLA_REFERENCE.get(
    "decoys", UMBRELLA_REFERENCE_ROOT + "/salmon/decoys.txt")

# ME_DB microexons missing from a GTF are added to a derived copy built here
# (src/insert_microexons_gtf.py); the configured GTFs are never changed.
# Whippet, the junction labels, HISAT2 splice hints, featureCounts and rMATS
# use the derived copies; Salmon and SUPPA2 keep salmon_gtf, whose transcripts
# must match UMBRELLA_TRANSCRIPTS. `umbrella_insert_microexons: false` uses
# the configured GTFs as they are.
UMBRELLA_INSERT_MICROEXONS = str(config.get("umbrella_insert_microexons", True)).lower() not in ("false", "f", "0", "no")
UMBRELLA_BASE_GTF = {"whippet": UMBRELLA_REFERENCE["whippet_gtf"],
                     "annotation": UMBRELLA_REFERENCE["annotation_gtf"]}
if UMBRELLA_INSERT_MICROEXONS:
    UMBRELLA_GTF = {base: UMBRELLA_REFERENCE_ROOT + "/{}.microexons.gtf.gz".format(base)
                    for base in UMBRELLA_BASE_GTF}
else:
    UMBRELLA_GTF = dict(UMBRELLA_BASE_GTF)
UMBRELLA_SPLICE_HINTS = UMBRELLA_REFERENCE_ROOT + "/hisat2_splice_sites.txt"


UMBRELLA_IDENTITY = UMBRELLA_REFERENCE_ROOT + "/identity.json"
_IDENTITY_PATHS, _IDENTITY_SETTINGS = reference_identity(UMBRELLA_REFERENCE, UMBRELLA_INSERT_MICROEXONS)


# Checks the manifest's reference_id against the configured inputs before any
# index is built, so a wrong ID fails in minutes rather than after the builds.
rule umbrella_reference_identity:
    input:
        [member for value in _IDENTITY_PATHS.values()
         for member in (value if isinstance(value, list) else [value])]
    output:
        UMBRELLA_IDENTITY
    run:
        manifest = reference_manifest(_IDENTITY_PATHS, expected_id=UMBRELLA_REFERENCE_ID,
                                      versions=UMBRELLA_REFERENCE.get("versions", {}),
                                      settings=_IDENTITY_SETTINGS)
        Path(output[0]).write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")


rule umbrella_microexon_gtf:
    input:
        gtf=lambda w: UMBRELLA_BASE_GTF[w.base],
        me_db=UMBRELLA_REFERENCE["me_db"],
        identity=UMBRELLA_IDENTITY
    output:
        gtf=protected(UMBRELLA_REFERENCE_ROOT + "/{base}.microexons.gtf.gz"),
        report=protected(UMBRELLA_REFERENCE_ROOT + "/{base}.microexons.report.tsv")
    wildcard_constraints:
        base="whippet|annotation"
    shell:
        "python3 src/insert_microexons_gtf.py --gtf {input.gtf} --me-db {input.me_db} "
        "--output {output.gtf} --report {output.report}"


rule umbrella_splice_hints:
    # All introns of the annotation GTF in use, including inserted microexons;
    # optionally union a supplied external hint set for backwards compatibility.
    input:
        gtf=UMBRELLA_GTF["annotation"],
        sites=[UMBRELLA_REFERENCE["splice_sites"]] if UMBRELLA_REFERENCE.get("splice_sites") else []
    output:
        protected(UMBRELLA_SPLICE_HINTS)
    run:
        from src.derive_umbrella_reference import splice_sites
        if input.sites:
            from src.summarize_junctions import _gtf_introns, _hint_introns
            introns = _gtf_introns(input.gtf) | _hint_introns(input.sites[0])
            Path(output[0]).write_text("".join(
                "{}\t{}\t{}\t{}\n".format(chrom, start - 2, end, strand)
                for chrom, start, end, strand in sorted(introns)))
        else:
            Path(output[0]).write_text(splice_sites(input.gtf))


rule umbrella_fixed_universe:
    # How the fixed microexon universe reached the Whippet annotation:
    # present already, inserted into a host transcript, or no host found.
    input:
        gtf=UMBRELLA_GTF["whippet"],
        me_db=UMBRELLA_REFERENCE["me_db"]
    output:
        protected(UMBRELLA_FIXED_UNIVERSE)
    run:
        from collections import Counter
        from src.insert_microexons_gtf import insert, read_me_db
        _, report = insert(input.gtf, read_me_db(input.me_db))
        counts = Counter(line.split("\t")[4] for line in report.splitlines()[1:])
        if not sum(counts.values()):
            raise WorkflowError("fixed microexon universe is empty")
        write_immutable(output[0], (json.dumps(
            {"microexons": sum(counts.values()), "in_whippet_gtf": counts["present"],
             "without_host": counts["no_host"]}, sort_keys=True) + "\n").encode())


rule umbrella_fixed_microexons:
    input:
        me_centric="Round2/TOTAL.ME_centric.txt",
        identity=UMBRELLA_IDENTITY
    output:
        protected(UMBRELLA_FIXED_MICROEXONS)
    run:
        from src.robustness_filter import write_fixed_microexons
        write_fixed_microexons(input.me_centric, output[0])


if "hisat2_index_prefix" not in UMBRELLA_REFERENCE:
    rule umbrella_hisat2_index:
        input:
            genome=UMBRELLA_REFERENCE["genome_fasta"],
            identity=UMBRELLA_IDENTITY
        output:
            UMBRELLA_HISAT_MEMBERS
        params:
            prefix=UMBRELLA_HISAT_PREFIX,
            large="--large-index" if UMBRELLA_HISAT_LARGE else ""
        threads: 8
        conda:
            "../envs/umbrella-quant.yaml"
        shell:
            "hisat2-build {params.large} -p {threads} {input.genome} {params.prefix}"


if "whippet_index" not in UMBRELLA_REFERENCE:
    rule umbrella_whippet_index:
        input:
            genome=UMBRELLA_REFERENCE["genome_fasta"],
            fixed_gtf=UMBRELLA_GTF["whippet"],
            validated=UMBRELLA_FIXED_UNIVERSE,
            identity=UMBRELLA_IDENTITY
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
        UMBRELLA_SALMON_GTF
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


if "transcriptome_fasta" not in UMBRELLA_REFERENCE:
    rule umbrella_transcriptome_fasta:
        input:
            gtf=UMBRELLA_SALMON_GTF,
            genome=UMBRELLA_REFERENCE["genome_fasta"],
            tx2gene=UMBRELLA_TX2GENE,
            identity=UMBRELLA_IDENTITY
        output:
            protected(UMBRELLA_TRANSCRIPTS)
        conda:
            "../envs/umbrella-quant.yaml"
        shell:
            "gffread -w {output} -g {input.genome} {input.gtf} && "
            "python3 src/derive_umbrella_reference.py validate-transcripts "
            "--fasta {output} --tx2gene {input.tx2gene}"


if "decoys" not in UMBRELLA_REFERENCE:
    rule umbrella_decoy_names:
        input:
            genome=UMBRELLA_REFERENCE["genome_fasta"],
            identity=UMBRELLA_IDENTITY
        output:
            protected(UMBRELLA_DECOYS)
        shell:
            "python3 src/derive_umbrella_reference.py decoys "
            "--input {input.genome} --output {output}"


if "salmon_index" not in UMBRELLA_REFERENCE:
    rule umbrella_gentrome:
        input:
            transcripts=UMBRELLA_TRANSCRIPTS,
            genome=UMBRELLA_REFERENCE["genome_fasta"],
            identity=UMBRELLA_IDENTITY
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
            decoys=UMBRELLA_DECOYS
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
        salmon_gtf=UMBRELLA_SALMON_GTF,
        transcripts=UMBRELLA_TRANSCRIPTS,
        decoys=UMBRELLA_DECOYS,
        hisat2_hints=UMBRELLA_SPLICE_HINTS,
        configured_splice_sites=[UMBRELLA_REFERENCE["splice_sites"]] if UMBRELLA_REFERENCE.get("splice_sites") else [],
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
        paths["hisat2_hints"] = str(input.hisat2_hints)
        if input.configured_splice_sites:
            paths["splice_sites"] = str(input.configured_splice_sites[0])
        if input.annotation_bed12:
            paths["annotation_bed12"] = str(input.annotation_bed12[0])
        # the ID covers configured inputs only, so `python3 src/shard_guard.py
        # reference-id --configfile <config>` gives it before anything is built
        identity, settings = reference_identity(UMBRELLA_REFERENCE, UMBRELLA_INSERT_MICROEXONS)
        manifest = reference_manifest(paths, expected_id=UMBRELLA_REFERENCE_ID,
                                      versions=UMBRELLA_REFERENCE.get("versions", {}),
                                      identity=identity, settings=settings)
        write_immutable(output[0], (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode())
