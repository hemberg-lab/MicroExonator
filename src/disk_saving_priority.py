"""Opt-in Snakemake 7 scheduling policy; no producer or file-lifecycle changes."""

# Explicit roles rather than automatic thread inspection: thread overrides must
# not accidentally promote acquisition or demote essential single-core gates.
PRIORITIES = {
    # Acquisition opens new fronts; keep it below already-ready consumers.
    'umbrella_stage_reads': 0,
    'umbrella_stage_single': 0,
    # Cheap gates can pin an entire library's reads if they are postponed.
    'umbrella_validate_reads': 110,
    'umbrella_legacy_bridge': 110,
    'umbrella_read_length': 110,
    'Round2_alingment_pre_processing': 110,
    'ME_reads': 110,
    # Parallel direct FASTQ readers and expensive BAM consumers.
    'umbrella_hisat2': 100,
    'umbrella_salmon_quant': 100,
    'umbrella_qapa_quant': 100,
    'umbrella_fastqc': 100,
    'umbrella_legacy_fastq': 100,
    'Round2_bowtie_to_tags': 100,
    'bowtie_to_genome': 100,
    'umbrella_featurecounts': 100,
    'umbrella_rmats_prep': 100,
    'umbrella_majiq_sj': 100,
    # Other direct readers: all must finish, not just the alignment.
    'umbrella_whippet_quant': 90,
    'umbrella_junctions': 90,
    'umbrella_coverage_run': 90,
    'umbrella_dapars2_coverage': 90,
    # Per-library cleanup chains for temporary ME intermediates.
    'Round2_filter': 80,
    'ME_SJ_coverage': 80,
}


def apply_disk_saving_priorities(workflow, config):
    """Assign numeric rule priorities after includes, preserving producer code.

    Absence/false preserves every existing priority. Apply only to umbrella
    manifests, including core-only tool selection. No custom dynamic scheduler.
    """
    value = config.get('umbrella_disk_saving_priority', False)
    text = str(value).strip().lower()
    if text not in {'true', 't', '1', 'yes', 'false', 'f', '0', 'no'}:
        raise ValueError('umbrella_disk_saving_priority must be true or false')
    if text in {'false', 'f', '0', 'no'} or 'umbrella_manifest' not in config:
        return
    for rule in workflow.rules:
        if rule.name in PRIORITIES:
            rule.priority = PRIORITIES[rule.name]
