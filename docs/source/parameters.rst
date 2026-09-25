.. _parameters:

=======================
Configuration reference
=======================

All parameters go in ``config.yaml`` inside the ``MicroExonator/`` folder. Boolean values are written as ``T`` or ``F``. Paths should be absolute. The genome and the annotation must come from the same source (for example UCSC, GENCODE, Ensembl or FlyBase), so that every chromosome in the annotation exists in the genome.

Required
========

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``Genome_fasta``
     - Multi-FASTA file with the genome sequence.
   * - ``Gene_anontation_bed12`` or ``Gene_anontation_GTF``
     - Transcript annotation, either as a 12-column BED file or as a GTF with gene, transcript and exon features (for example from GENCODE or Ensembl). At least one is required; the GTF is also needed for the Whippet modules.
   * - ``working_directory``
     - Path to the cloned ``MicroExonator/`` folder.
   * - ``min_number_files_detected``
     - Minimum number of samples in which a putative novel microexon must be found during discovery to be kept for quantification. We recommend at least 2 for single-end data and 3 when paired-end files are present.
   * - ``bulk_samples``
     - Tab-separated file assigning each bulk sample to a condition (columns ``sample`` and ``condition``; see :doc:`setup`). Required for bulk RNA-seq runs; pure single-cell runs use the cluster metadata instead.

Splice sites, conservation and microexon length
===============================================

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Default
     - Description
   * - ``ME_len``
     - 30
     - Maximum microexon length in nt.
   * - ``GT_AG_U2_5``, ``GT_AG_U2_3``
     - NA
     - Position weight matrices of U2 GT-AG 5' and 3' splice sites, used to score candidate microexons. Examples for human and mouse are in the ``PWM/`` folder; `SpliceRack <http://katahdin.cshl.edu/splice/>`_ has them for other species. With ``NA``, MicroExonator builds the matrices from annotated splice sites.
   * - ``conservation_bigwig``
     - NA
     - BigWig file of genome-wide conservation scores (phyloP or PhastCons, for example from the UCSC Genome Browser). Highly conserved candidate microexons are protected from being discarded. With ``NA``, conservation is not used.
   * - ``min_conservation``
     - 2
     - Mean conservation score above which a candidate microexon is considered conserved, when ``conservation_bigwig`` is given.

Input data and disk use
=======================

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Default
     - Description
   * - ``paired_samples``
     - none
     - Tab-separated file with two columns and no header, pairing the two files of each paired-end sample. Counts from both files are pooled and reported under the name in the first column.
   * - ``Keep_fastq_gz``
     - F
     - Keep the downloaded or copied FASTQ files in ``FASTQ/`` after they are processed. By default they are deleted.
   * - ``Optimize_hard_drive``
     - F
     - Keep separate temporary FASTQ copies for discovery and for quantification and delete each as soon as it is used. This minimises disk use at the cost of downloading or copying each file twice. See "Running large datasets" in :doc:`discovery_and_quantification`.
   * - ``max_read_len``
     - 150
     - Length of exon sequence on each side of every splice-junction tag, and the length reads are trimmed to before they are aligned to the tags. Any read up to ``max_read_len`` + 8 nt is counted without bias, so one value serves all shorter read lengths. See :ref:`read_length`.

Quantification
==============

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Default
     - Description
   * - ``ME_DB``
     - none
     - File of known microexons to quantify in addition to the annotated and discovered ones, for example from `VastDB <https://vastdb.crg.eu/>`_. BED12 rows, BED6 rows and single-column ``chrom_strand_start_end`` IDs can be mixed in one file (see :doc:`discovery_and_quantification`).
   * - ``min_reads_PSI``
     - 5
     - Minimum number of junction reads (inclusion plus exclusion) needed to report a PSI value; below it PSI is ``NA``.
   * - ``synthetic_skip_tags``
     - T
     - Add a skipping tag for every intron that contains a microexon but has no annotated exon-exon junction. Without it, microexons that are only annotated as included have no exclusion evidence and PSI stays at 1. See :ref:`skipping_tags`.
   * - ``consecutive_microexons``
     - T
     - Correct inclusion and exclusion counts for microexons next to another microexon. See :ref:`consecutive`.
   * - ``shared_splice_site_correction``
     - T
     - Count microexons that share one splice site with a longer annotated exon from their unshared side, and count the longer exon once as exclusion. See :ref:`shared_site`.
   * - ``mismatch_aware_blacklist``
     - T
     - Remove a read that also aligns to the genome without splicing only if the genome alignment is at least as good as the tag alignment. With ``F``, every read with a genome alignment is removed. See :ref:`blacklist`.
   * - ``skip_discovery``
     - F
     - Skip the discovery module and quantify only annotated microexons and those in ``ME_DB``.
   * - ``only_db``
     - F
     - Build the quantification tags from the annotation and ``ME_DB`` only. Usually combined with ``skip_discovery : T``.
   * - ``split_DB``
     - F
     - With ``only_db``, split ``ME_DB`` into chunks that are processed in parallel. Useful for very large microexon databases.
   * - ``skip_genome_alignment``
     - F
     - Do not remove reads that align to the genome without splicing. Not recommended except for testing.

Confidence filtering
====================

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Default
     - Description
   * - ``filter_method``
     - robustness
     - ``robustness`` keeps microexons with a confident PSI in most samples of at least one group; ``legacy_mixture`` uses the Gaussian mixture model of the original publication. See :doc:`discovery_and_quantification`.
   * - ``min_detected_samples``
     - 1
     - For ``robustness``: number of measurements within one group that must meet the criteria.

Differential inclusion (Whippet)
================================

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Default
     - Description
   * - ``whippet_bin_folder``
     - none
     - Path to the ``bin/`` folder of the Whippet installation.
   * - ``julia``
     - none
     - Path to the Julia executable compatible with Whippet.
   * - ``whippet_delta``
     - none
     - YAML file listing the comparisons between sample groups (see :doc:`differential_inclusion_analysis`).
   * - ``downstream_only``
     - F
     - Skip discovery and quantification and run Whippet only on the annotation.
   * - ``use_uncorrected_PSI``
     - F
     - Give Whippet the uncorrected MicroExonator PSI values instead of the corrected ones (for reproducing older analyses).

Single-cell analysis
====================

.. list-table::
   :header-rows: 1
   :widths: 30 12 58

   * - Parameter
     - Default
     - Description
   * - ``Single_Cell``
     - F
     - Enable the single-cell modules.
   * - ``cluster_metadata``
     - none
     - Tab-separated file with a header, assigning each single-cell sample to a cluster (cell type).
   * - ``cluster_name``, ``file_basename``
     - none
     - Names of the columns of ``cluster_metadata`` that hold the cluster and the sample name.
   * - ``cells_pseudobulks``
     - 15
     - Target number of cells per pseudo-bulk; the number of pseudo-bulks per cluster is the number of cells divided by this value.
   * - ``n_pseudobulks``
     - from ``cells_pseudobulks``
     - Fixed number of pseudo-bulks per cluster, overriding ``cells_pseudobulks``.
   * - ``min_number_of_reads_single_cell``
     - none
     - Minimum number of reads to compute PSI in pseudo-bulk comparisons (the single-cell counterpart of ``min_reads_PSI``).
   * - ``min_number_of_samples_single_cell``
     - none
     - Minimum number of pseudo-bulks per group in which a node must be quantified.
   * - ``run_metadata``
     - none
     - Tab-separated file describing the comparisons between cell types (see :doc:`single_cell_analysis`).
   * - ``cdf_t``
     - none (0.8 recommended)
     - Probability threshold for the Beta-distribution test of differential inclusion across computational replicates.
   * - ``min_p_mean``
     - none (0.9 recommended)
     - Minimum mean probability of differential inclusion.
   * - ``min_delta``
     - none (0.1 recommended)
     - Minimum mean ΔPSI.
   * - ``min_rep``
     - none (25 recommended)
     - Minimum number of computational replicates in which a node could be tested.
   * - ``seed``
     - 123
     - Seed for the random assignment of cells to pseudo-bulks.
   * - ``Only_snakepool``
     - F
     - Skip discovery and quantification and test all Whippet splicing nodes.
   * - ``Get_Bamfiles``
     - F
     - Have ``whippet-quant`` write alignments so that BAM files per cell type can be generated for visualisation.
