.. _whats_new:

==========
What's new
==========

This version ports MicroExonator to Python 3 and fixes several quantification problems that affected a sizeable share of microexons. Every change is covered by tests, including an end-to-end test on simulated reads with known PSI. The ``MiMB`` branch, which accompanies the *Methods in Molecular Biology* chapter, is unchanged.

Python 3 and defaults
=====================

* All pipeline scripts run under Python 3, with pinned conda environments.
* The robustness filter is the default confidence filter; the original Gaussian mixture filter is available as ``filter_method : legacy_mixture``.
* Whippet receives corrected MicroExonator PSI values by default.

Quantification fixes
====================

Each fix can be switched off to reproduce the counts of earlier versions. See :doc:`how_it_works` for details.

* **Read length.** Tag flanks were fixed at 100 nt and reads were not trimmed, so reads longer than 100 nt over-counted inclusion (increasingly with microexon length) and most of their junction reads were lost. Tag flanks now follow ``max_read_len`` (default 150) and longer reads are trimmed to it. On simulated 150 nt reads, the old version overestimated PSI by +0.02 for microexons under 10 nt, +0.06 for 10 to 19 nt and +0.085 for 20 to 30 nt.
* **Microexons without an annotated skipping junction.** Many microexons are annotated only as included, so their exclusion reads could not be counted and PSI stayed at 1. A skipping tag is now added for each such intron (``synthetic_skip_tags``).
* **Consecutive microexons.** Skipping reads of one microexon that land on its neighbour's inclusion junction, and reads on the junction between two microexons, are now counted for both (``consecutive_microexons``).
* **Junctions with a copy in the genome.** Reads that also align to the genome are now removed only if they fit the genome at least as well as the junction (``mismatch_aware_blacklist``).
* **Mixed** ``ME_DB`` **files.** BED6 rows are read correctly, and single-column rows no longer change the length limit for the rows after them.

New reports
===========

* ``Report/read_lengths.tsv``: read lengths per sample, with warnings for samples that look like another library type.
* ``Report/ME_junction_genomic_copies.txt``: microexons whose junctions have an unspliced copy in the genome.
* ``Report/ME_ambiguous_positions.txt``: short microexons whose position in the intron cannot be resolved from reads.
* ``Report/ME_consecutive_runs.txt``: microexons adjacent to other microexons.

Testing
=======

* Unit tests cover the counting rules, tag generation and configuration handling.
* An mm10 smoke test runs the whole workflow on simulated reads with known PSI, including loci with two and three consecutive microexons, and checks that the mean PSI of every event in every sample group is within 0.15 of its simulated value.
