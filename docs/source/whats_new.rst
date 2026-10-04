.. _whats_new:

==========
What's new
==========

MicroExonator 2.0.0-beta.1
==========================

This version ports MicroExonator to Python 3 and fixes several quantification problems that affected a sizeable share of microexons. Every change is covered by tests, including an end-to-end test on simulated reads with known PSI. The ``MiMB`` branch, which accompanies the *Methods in Molecular Biology* chapter, is unchanged.

Python 3 and defaults
---------------------

* All pipeline scripts run under Python 3, with pinned conda environments.
* The robustness filter is the default confidence filter; it needs a ``bulk_samples`` file with the condition of each bulk sample. The original Gaussian mixture filter is available as ``filter_method : legacy_mixture`` and needs no ``bulk_samples`` file.
* Whippet receives corrected MicroExonator PSI values by default.

Quantification fixes
--------------------

Each fix can be switched off to reproduce the counts of earlier versions. See :doc:`how_it_works` for details.

* **Read length.** Tag flanks were fixed at 100 nt and reads were not trimmed, so reads longer than 100 nt over-counted inclusion (increasingly with microexon length) and most of their junction reads were lost. Tag flanks now follow ``max_read_len`` (default 150) and longer reads are trimmed to it. On simulated 150 nt reads, the old version overestimated PSI by +0.02 for microexons under 10 nt, +0.06 for 10 to 19 nt and +0.085 for 20 to 30 nt.
* **Microexons without an annotated skipping junction.** Many microexons are annotated only as included, so their exclusion reads could not be counted and PSI stayed at 1. A skipping tag is now added for each such intron (``synthetic_skip_tags``).
* **Consecutive microexons.** Skipping reads of one microexon that land on its neighbour's inclusion junction, and reads on the junction between two microexons, are now counted for both (``consecutive_microexons``).
* **Microexons sharing a splice site with a longer exon.** Inclusion is now counted from the microexon's unshared side, and the longer exon counts once as exclusion, through its far junction. Previously its reads were counted over both of its introns and once per annotated exon sharing them, and were also subtracted from inclusion, so PSI was underestimated (by 0.12 on average in simulations, up to 0.37) (``shared_splice_site_correction``).
* **Junctions with a copy in the genome.** Reads that also align to the genome are now removed only if they fit the genome at least as well as the junction (``mismatch_aware_blacklist``).
* **Mixed** ``ME_DB`` **files.** BED6 rows are read correctly, and single-column rows no longer change the length limit for the rows after them.

Differential inclusion without Whippet
--------------------------------------

* Differential inclusion uses MicroExonator's own test by default; ``delta_method : whippet`` selects the Whippet route. ``downstream_only : T`` on its own now runs Whippet on the annotation alone, as the *Methods in Molecular Biology* chapter describes; it used to also need ``Only_whippet : T``. ``delta_method : microexonator`` tests differential inclusion of microexons with MicroExonator's own PSI and read counts, using Whippet's statistical model, so Julia and Whippet do not need to be installed. With the same inputs it reproduces Whippet's ``DeltaPsi`` to within 0.02. See :doc:`differential_inclusion_analysis`.

New reports
-----------

* ``Report/read_lengths.tsv``: read lengths per sample, with warnings for samples that look like another library type.
* ``Report/ME_junction_genomic_copies.txt``: microexons whose junctions have an unspliced copy in the genome.
* ``Report/ME_ambiguous_positions.txt``: short microexons whose position in the intron cannot be resolved from reads.
* ``Report/ME_consecutive_runs.txt``: microexons adjacent to other microexons.

Input and cluster setup
-----------------------

* Local FASTQ files are copied with a bare ``+`` separator line and read names limited to printable ASCII, which Whippet and Bowtie need; clean ``.fastq.gz`` files are copied as they are. See :doc:`setup`.
* Paired-end mates joined from SRA downloads also get a bare ``+`` separator line.
* A PBS cluster example (``Examples/Cluster_config/pbs/cluster.PBS.json``) with resources measured on hg38 runs, and a fixed LSF example. ``Round2_filter`` memory grows with read depth; see :ref:`cluster_resources`.
* Start-up no longer prints Python ``DeprecationWarning`` messages from other packages; deprecated configuration keys are reported as one ``Note:`` line.
* Faster quantification: reading the microexon reads from each FASTQ (``ME_reads``) runs about 10 times faster, finding the introns that contain each microexon (``Get_ME_from_annotation``) about 10 times faster, and the two Bowtie steps of quantification use 8 threads. Outputs are unchanged.
* ``Optimize_hard_drive : T`` works again: quantification fetches its own temporary copy of each FASTQ after discovery and every quantification step reads it, so nothing is downloaded twice. ``validate_fastq_list`` works again for quantification; a duplicate helper in ``Round2.smk`` pointed it at files no rule produced.
* The corrected PSI files use ``min_reads_PSI`` as their read cut-off; it was fixed at 5 before.
* The single-cell targets of the *Methods in Molecular Biology* chapter (``snakepool``, ``quant_unpool_single_cell``, ``collapse_whippet`` and ``cluster_bams``) work again; they had been switched off by mistake. Each comparison now keeps its own number of pseudo-bulks and repeats, and differentially included nodes are matched against the microexons of the configured confidence filter.
* The experimental Google Cloud Storage input (``google_path``) has been removed, so the Google Cloud packages are no longer needed.

Testing
-------

* Unit tests cover the counting rules, tag generation and configuration handling.
* An mm10 smoke test runs the whole workflow on simulated reads with known PSI, including loci with two and three consecutive microexons, and checks that the mean PSI of every event in every sample group is within 0.15 of its simulated value.
