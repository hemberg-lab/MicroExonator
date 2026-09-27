"""One shared HISAT2 alignment per run, and the analyses that read it.

Every consumer of the temporary BAM (junctions, featureCounts, rMATS prep,
coverage) is a declared rule, so Snakemake deletes the BAM only after all of
them have run. Durable outputs are group-by-batch shards written through
shard_guard; nothing per run is kept except the rMATS prep files, which
`--task post` needs later.

Optional analyses (config `umbrella_optional`, all on by default):
  rmats       rMATS-turbo prep per run, kept for later per-comparison post
  coverage    summed CPM coverage per group and batch (bigWig) plus n
  leafcutter  per-comparison LeafCutter (also needs config leafcutter_dir)
  suppa2      per-comparison SUPPA2 on Salmon transcript TPM
"""

import json
from pathlib import Path

from src.shard_guard import sha256_file


UMBRELLA_OPTIONAL = {"rmats": True, "coverage": True, "leafcutter": True, "suppa2": True}
UMBRELLA_OPTIONAL.update(config.get("umbrella_optional", {}) or {})
UMBRELLA_HISAT2_FLAGS = config.get("umbrella_hisat2_flags", "")
UMBRELLA_INTRON_CATALOG = UMBRELLA_REFERENCE_ROOT + "/introns.tsv.gz"
UMBRELLA_CHROM_SIZES = UMBRELLA_REFERENCE_ROOT + "/chrom.sizes"
UMBRELLA_WORK = "umbrella/work/{reference_id}/{project_id}/{batch_id}/{run_id}"
UMBRELLA_SHARD = "{reference_id}/{project_id}/{group}/{batch_id}"


def umbrella_align_run(wildcards):
    record = UMBRELLA_MANIFEST.by_run.get(wildcards.run_id)
    if record is None or not record.include:
        raise WorkflowError("umbrella alignment run is absent or excluded")
    if wildcards.reference_id != UMBRELLA_REFERENCE_ID:
        raise WorkflowError("umbrella alignment reference_id conflicts with frozen bundle")
    for field in ("project_id", "batch_id", "group"):
        if hasattr(wildcards, field) and getattr(wildcards, field) != getattr(record, field):
            raise WorkflowError("umbrella alignment path conflicts with manifest")
    return record


def umbrella_shard_runs(wildcards):
    if wildcards.reference_id != UMBRELLA_REFERENCE_ID:
        raise WorkflowError("umbrella shard reference_id conflicts with frozen bundle")
    records = UMBRELLA_MANIFEST.runs_for(wildcards.project_id, wildcards.group, wildcards.batch_id)
    if not records:
        raise WorkflowError("umbrella shard has no included runs")
    return records


def umbrella_hisat2_reads(wildcards):
    reads = UMBRELLA_MANIFEST.native_reads(umbrella_align_run(wildcards).run_id)
    return "-1 {} -2 {}".format(*reads) if len(reads) == 2 else "-U {}".format(reads[0])


def umbrella_hisat2_strand(wildcards):
    record = umbrella_align_run(wildcards)
    if record.strandedness in ("unstranded", "auto"):
        return ""
    return "--rna-strandness " + {"firststrand": "R", "secondstrand": "F"}[
        record.strandedness] * (2 if record.layout == "PE" else 1)


# ------------------------------------------------------------ reference level

rule umbrella_intron_catalog:
    input:
        whippet_gtf=UMBRELLA_REFERENCE["whippet_gtf"],
        annotation_gtf=UMBRELLA_REFERENCE["annotation_gtf"],
        splice_sites=UMBRELLA_REFERENCE["splice_sites"]
    output:
        protected(UMBRELLA_INTRON_CATALOG)
    shell:
        "python3 src/summarize_junctions.py catalog --whippet-gtf {input.whippet_gtf} "
        "--annotation-gtf {input.annotation_gtf} --splice-sites {input.splice_sites} "
        "--output {output}"


rule umbrella_chrom_sizes:
    input:
        UMBRELLA_REFERENCE["genome_fasta"]
    output:
        protected(UMBRELLA_CHROM_SIZES)
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "samtools faidx {input} --fai-idx - | cut -f1,2 > {output}"


# ------------------------------------------------------------ per run

rule umbrella_hisat2:
    input:
        reads=lambda w: UMBRELLA_MANIFEST.native_reads(umbrella_align_run(w).run_id),
        valid=UMBRELLA_WORK + "/reads.valid",
        index=UMBRELLA_HISAT_MEMBERS,
        splice_sites=UMBRELLA_REFERENCE["splice_sites"],
        reference=UMBRELLA_REFERENCE_MANIFEST
    output:
        bam=temp(UMBRELLA_WORK + "/align/aligned.bam"),
        bai=temp(UMBRELLA_WORK + "/align/aligned.bam.bai"),
        summary=temp(UMBRELLA_WORK + "/align/hisat2.summary.txt")
    params:
        prefix=UMBRELLA_HISAT_PREFIX,
        reads=umbrella_hisat2_reads,
        strand=umbrella_hisat2_strand,
        flags=UMBRELLA_HISAT2_FLAGS
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.hisat2.log"
    threads: 8
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "hisat2 -p {threads} -x {params.prefix} {params.reads} "
        "--known-splicesite-infile {input.splice_sites} {params.strand} "
        "--new-summary --summary-file {output.summary} {params.flags} 2> {log} "
        "| samtools sort -@ {threads} -o {output.bam} - 2>> {log} "
        "&& samtools index {output.bam}"


rule umbrella_junctions:
    input:
        bam=UMBRELLA_WORK + "/align/aligned.bam",
        bai=UMBRELLA_WORK + "/align/aligned.bam.bai"
    output:
        table=temp(UMBRELLA_WORK + "/align/junctions.tsv.gz"),
        summary=temp(UMBRELLA_WORK + "/align/junctions.summary.json")
    params:
        min_anchor=config.get("umbrella_min_anchor", 8)
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "samtools view {input.bam} | python3 src/summarize_junctions.py extract "
        "--table {output.table} --summary {output.summary} --min-anchor {params.min_anchor}"


rule umbrella_featurecounts:
    input:
        bam=UMBRELLA_WORK + "/align/aligned.bam",
        bai=UMBRELLA_WORK + "/align/aligned.bam.bai",
        summary=UMBRELLA_WORK + "/align/junctions.summary.json",
        gtf=UMBRELLA_REFERENCE["annotation_gtf"]
    output:
        counts=temp(UMBRELLA_WORK + "/align/featurecounts.txt"),
        summary=temp(UMBRELLA_WORK + "/align/featurecounts.txt.summary")
    params:
        declared=lambda w: umbrella_align_run(w).strandedness,
        paired=lambda w: "-p --countReadPairs" if umbrella_align_run(w).layout == "PE" else ""
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.featurecounts.log"
    threads: 4
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "featureCounts -T {threads} -a {input.gtf} -t exon -g gene_id {params.paired} "
        "-s $(python3 src/summarize_junctions.py strandedness --tool featurecounts "
        "--declared {params.declared} --summary {input.summary}) "
        "-o {output.counts} {input.bam} 2> {log}"


rule umbrella_read_length:
    input:
        lambda w: UMBRELLA_MANIFEST.native_reads(umbrella_align_run(w).run_id)[0]
    output:
        temp(UMBRELLA_WORK + "/align/read_length.txt")
    run:
        import gzip as _gzip
        longest = 0
        with _gzip.open(input[0], "rt") as stream:
            for number, line in enumerate(stream):
                if number % 4 == 1:
                    longest = max(longest, len(line.rstrip("\n")))
                if number >= 40000:
                    break
        if not longest:
            raise WorkflowError("no reads in {}".format(input[0]))
        Path(output[0]).write_text("{}\n".format(longest))


rule umbrella_rmats_prep:
    input:
        bam=lambda w: umbrella_align_run(w).work_dir + "/align/aligned.bam",
        bai=lambda w: umbrella_align_run(w).work_dir + "/align/aligned.bam.bai",
        read_length=lambda w: umbrella_align_run(w).work_dir + "/align/read_length.txt",
        gtf=UMBRELLA_REFERENCE["annotation_gtf"]
    output:
        prep=protected("rmats/" + UMBRELLA_SHARD + "/{run_id}.rmats"),
        record=protected("rmats/" + UMBRELLA_SHARD + "/{run_id}.rmats.json")
    params:
        work=lambda w: umbrella_align_run(w).work_dir + "/rmats",
        layout=lambda w: "paired" if umbrella_align_run(w).layout == "PE" else "single",
        manifest_layout=lambda w: umbrella_align_run(w).layout
    log:
        "umbrella/logs/{reference_id}/{project_id}/{batch_id}/{run_id}.{group}.rmats_prep.log"
    threads: 4
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "rm -rf {params.work} && mkdir -p {params.work}/tmp {params.work}/od "
        "&& echo {input.bam} > {params.work}/b1.txt "
        "&& rmats.py --b1 {params.work}/b1.txt --gtf {input.gtf} --od {params.work}/od "
        "--tmp {params.work}/tmp -t {params.layout} --readLength $(cat {input.read_length}) "
        "--variable-read-length --novelSS --nthread {threads} --task prep > {log} 2>&1 "
        "&& python3 src/inventory_rmats_prep.py keep --tmp {params.work}/tmp "
        "--destination {output.prep} --run-id {wildcards.run_id} --bam {input.bam} "
        "--layout {params.manifest_layout} --read-length $(cat {input.read_length}) "
        "&& rm -rf {params.work}"


rule umbrella_coverage_run:
    input:
        bam=UMBRELLA_WORK + "/align/aligned.bam",
        bai=UMBRELLA_WORK + "/align/aligned.bam.bai"
    output:
        temp(UMBRELLA_WORK + "/align/coverage.cpm.bedGraph")
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        "mapped=$(samtools view -c -F 0x904 {input.bam}) "
        "&& bedtools genomecov -ibam {input.bam} -split -bg "
        "-scale $(python3 -c \"print(1e6 / max(1, $mapped))\") "
        "| LC_ALL=C sort -k1,1 -k2,2n > {output}"


# ------------------------------------------------------------ per group and batch

def umbrella_run_files(relative):
    def paths(wildcards):
        return [record.work_dir + relative for record in umbrella_shard_runs(wildcards)]
    return paths


def umbrella_manifest_sha():
    return sha256_file(str(UMBRELLA_MANIFEST.path))


rule umbrella_junction_shard:
    input:
        tables=umbrella_run_files("/align/junctions.tsv.gz"),
        catalog=UMBRELLA_INTRON_CATALOG,
        reference=UMBRELLA_REFERENCE_MANIFEST,
        manifest=str(UMBRELLA_MANIFEST.path)
    output:
        shard=protected("junctions/" + UMBRELLA_SHARD + ".junctions.tsv.gz"),
        checksums=protected("junctions/" + UMBRELLA_SHARD + ".junctions_checksums.json")
    params:
        tables=lambda w: " ".join("--table {}={}".format(record.run_id, record.work_dir + "/align/junctions.tsv.gz")
                                  for record in umbrella_shard_runs(w)),
        prefix="junctions/" + UMBRELLA_SHARD
    run:
        shell("python3 src/summarize_junctions.py shard --catalog {input.catalog} "
              "--prefix {params.prefix} --reference-id {wildcards.reference_id} "
              "--manifest-sha256 " + umbrella_manifest_sha() + " {params.tables}")


rule umbrella_capture_shard:
    input:
        shard="junctions/" + UMBRELLA_SHARD + ".junctions.tsv.gz",
        catalog=UMBRELLA_INTRON_CATALOG,
        manifest=str(UMBRELLA_MANIFEST.path)
    output:
        capture=protected("junctions/" + UMBRELLA_SHARD + ".capture.tsv.gz"),
        rates=protected("junctions/" + UMBRELLA_SHARD + ".capture_rates.json"),
        checksums=protected("junctions/" + UMBRELLA_SHARD + ".capture_checksums.json")
    params:
        prefix="junctions/" + UMBRELLA_SHARD
    run:
        shell("python3 src/summarize_junctions.py capture --catalog {input.catalog} "
              "--shard {input.shard} --prefix {params.prefix} "
              "--reference-id {wildcards.reference_id} --manifest-sha256 " + umbrella_manifest_sha())


rule umbrella_featurecounts_shard:
    input:
        counts=umbrella_run_files("/align/featurecounts.txt"),
        reference=UMBRELLA_REFERENCE_MANIFEST,
        manifest=str(UMBRELLA_MANIFEST.path)
    output:
        counts=protected("genes/" + UMBRELLA_SHARD + ".featurecounts.tsv.gz"),
        checksums=protected("genes/" + UMBRELLA_SHARD + ".featurecounts_checksums.json")
    run:
        from src.reduce_featurecounts import reduce_group as reduce_featurecounts
        reduce_featurecounts(
            {record.run_id: record.work_dir + "/align/featurecounts.txt"
             for record in umbrella_shard_runs(wildcards)},
            "genes/{}/{}/{}/{}".format(wildcards.reference_id, wildcards.project_id,
                                       wildcards.group, wildcards.batch_id),
            wildcards.reference_id, umbrella_manifest_sha())


rule umbrella_qc_shard:
    input:
        hisat2=umbrella_run_files("/align/hisat2.summary.txt"),
        junctions=umbrella_run_files("/align/junctions.summary.json"),
        featurecounts=umbrella_run_files("/align/featurecounts.txt.summary"),
        manifest=str(UMBRELLA_MANIFEST.path)
    output:
        qc=protected("qc/" + UMBRELLA_SHARD + ".qc.tsv.gz"),
        checksums=protected("qc/" + UMBRELLA_SHARD + ".qc_checksums.json")
    run:
        from src.reduce_qc import reduce_group as reduce_qc
        runs = [{"run_id": record.run_id, "layout": record.layout,
                 "strandedness": record.strandedness,
                 "hisat2_summary": record.work_dir + "/align/hisat2.summary.txt",
                 "junction_summary": record.work_dir + "/align/junctions.summary.json",
                 "featurecounts_summary": record.work_dir + "/align/featurecounts.txt.summary"}
                for record in umbrella_shard_runs(wildcards)]
        reduce_qc(runs, "qc/{}/{}/{}/{}".format(wildcards.reference_id, wildcards.project_id,
                                                 wildcards.group, wildcards.batch_id),
                  wildcards.reference_id, umbrella_manifest_sha())


rule umbrella_rmats_inventory:
    input:
        records=lambda w: ["rmats/{}/{}/{}/{}/{}.rmats.json".format(
            w.reference_id, w.project_id, w.group, w.batch_id, record.run_id)
            for record in umbrella_shard_runs(w)],
        manifest=str(UMBRELLA_MANIFEST.path)
    output:
        inventory=protected("rmats/" + UMBRELLA_SHARD + ".rmats_inventory.tsv"),
        checksums=protected("rmats/" + UMBRELLA_SHARD + ".rmats_checksums.json")
    params:
        prefix="rmats/" + UMBRELLA_SHARD
    run:
        shell("python3 src/inventory_rmats_prep.py inventory --prefix {params.prefix} "
              "--reference-id {wildcards.reference_id} --manifest-sha256 "
              + umbrella_manifest_sha() + " {input.records}")


rule umbrella_coverage_shard:
    input:
        tracks=umbrella_run_files("/align/coverage.cpm.bedGraph"),
        sizes=UMBRELLA_CHROM_SIZES,
        manifest=str(UMBRELLA_MANIFEST.path)
    output:
        bigwig=protected("coverage/" + UMBRELLA_SHARD + ".sum_cpm.bw"),
        tally=protected("coverage/" + UMBRELLA_SHARD + ".sum_cpm.json")
    params:
        runs=lambda w: json.dumps(sorted(record.run_id for record in umbrella_shard_runs(w)))
    conda:
        "../envs/umbrella-alignment.yaml"
    shell:
        # One track: summed CPM coverage. The group mean is the sum over all of
        # a group's batch shards divided by the total n recorded next to each.
        "tmp=$(mktemp {output.bigwig}.XXXXXX.bedGraph) && "
        "if [ $(echo {input.tracks} | wc -w) -gt 1 ]; then "
        "bedtools unionbedg -i {input.tracks} "
        "| awk 'BEGIN{{OFS=\"\\t\"}} {{s=0; for (i=4; i<=NF; i++) s+=$i; if (s>0) print $1,$2,$3,s}}' > $tmp; "
        "else cp {input.tracks} $tmp; fi "
        "&& bedGraphToBigWig $tmp {input.sizes} {output.bigwig} && rm -f $tmp "
        "&& echo '{{\"n\": '$(echo {input.tracks} | wc -w)', \"runs\": {params.runs}}}' > {output.tally}"


UMBRELLA_ALIGNMENT_TARGETS = []
for reference_id, project_id, group, batch_id in sorted({
        (run.reference_id, run.project_id, run.group, run.batch_id)
        for run in UMBRELLA_MANIFEST.included_runs()}):
    shard = "{}/{}/{}/{}".format(reference_id, project_id, group, batch_id)
    UMBRELLA_ALIGNMENT_TARGETS += [
        "junctions/{}.junctions.tsv.gz".format(shard),
        "junctions/{}.capture.tsv.gz".format(shard),
        "genes/{}.featurecounts.tsv.gz".format(shard),
        "qc/{}.qc.tsv.gz".format(shard),
    ]
    if UMBRELLA_OPTIONAL.get("rmats"):
        UMBRELLA_ALIGNMENT_TARGETS.append("rmats/{}.rmats_inventory.tsv".format(shard))
    if UMBRELLA_OPTIONAL.get("coverage"):
        UMBRELLA_ALIGNMENT_TARGETS.append("coverage/{}.sum_cpm.bw".format(shard))
