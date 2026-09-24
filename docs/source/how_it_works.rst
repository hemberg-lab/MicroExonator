.. _how_it_works:

============================
How quantification works
============================

This page explains how MicroExonator turns reads into PSI values, and what each quantification option changes. It is not needed to run the workflow, but it helps to interpret the results.

Two rounds of alignment
=======================

Genome aligners often fail to map reads over microexons, because a microexon leaves only a few nucleotides of the read on its own exon. MicroExonator avoids this by aligning reads to short sequences, called tags, that represent splice junctions directly.

**Discovery (Round1).** Reads are aligned with BWA-MEM to tags that join annotated exon ends (exon-exon junction tags). A read that comes from a transcript with an unannotated microexon aligns to the junction tag with a short insertion at the junction. These insertions are collected, placed back in the intron where they are flanked by canonical splice sites (AG...GT), scored with the U2 splice-site matrices and conservation, and filtered to give a list of putative novel microexons.

**Quantification (Round2).** For every annotated, database and discovered microexon, MicroExonator builds an *inclusion tag* (upstream exon, microexon, downstream exon) for each intron that contains it. Reads are aligned end to end with Bowtie to these tags together with the *skipping tags* (upstream exon joined directly to downstream exon). Reads that also align to the genome without splicing are then removed (see :ref:`blacklist`).

For each microexon and sample, reads crossing an inclusion junction with at least 8 nt on each side count as inclusion, and reads on the skipping junctions of its containing introns count as exclusion. Reads on alternative 5' or 3' splice sites of the microexon count as exclusion too. PSI is inclusion / (inclusion + exclusion), with an exact binomial 95% confidence interval, and is reported only when at least ``min_reads_PSI`` reads are available.

**Corrected PSI.** An inclusion read can cross one or both junctions of a microexon, while a skipping read crosses one junction. To keep the two comparable, reads that span the whole microexon count as one, and reads that cross only one of its junctions count as half. These corrected values (``Report/quant/corrected/``) are the ones used by the confidence filter and passed to Whippet.

.. _read_length:

Tag length and read length
==========================

Reads are aligned end to end, so a read is counted only if it fits entirely inside a tag. Each tag carries ``max_read_len`` nt of exon on each side of the junction (default 150). A read of length *L* that crosses a junction with the minimum 8 nt anchors ends at most *L* - 8 nt from the junction, so every read up to ``max_read_len`` + 8 nt fits in every position where it crosses a junction, and inclusion and exclusion are sampled equally.

If reads are longer than the flanks, they fit only in some positions, and an inclusion tag has *m* more positions than a skipping tag, where *m* is the microexon length. Inclusion is then over-counted by a factor that grows with microexon length (up to 1.6x for a 30 nt microexon with 150 nt reads on 100 nt flanks), and most junction reads are lost. To prevent this, reads longer than ``max_read_len`` are trimmed to that length before they are aligned to the tags.

Tag flanks are cut from spliced transcripts, so when an exon is shorter than the flank, a tag's flank also contains the next junction, and a read from that junction fits equally well inside the neighbouring tag's flank, where it crosses no junction. To keep such reads, Bowtie reports up to 10 alignments per read, and MicroExonator keeps one alignment per read that crosses its tag's junction with at least 8 nt on each side (fewest mismatches first, ties broken by a fixed seed). This makes the counts independent of the flank length.

In practice, set ``max_read_len`` to the length of the longest reads in your data. Shorter reads in the same run are handled correctly. ``Report/read_lengths.tsv`` shows the read lengths of each sample and flags samples that look like a different library type (for example single-cell data included by mistake).

.. _skipping_tags:

Microexons without an annotated skipping junction
=================================================

Skipping tags come from annotated exon-exon junctions. Many microexons are annotated only in transcripts that include them, so the junction that skips them is not in the annotation and there is no tag to count exclusion reads on; their PSI would then always be 1.

With ``synthetic_skip_tags : T`` (default), MicroExonator adds one skipping tag for each intron that contains a microexon but has no skipping tag, built from the same exon flanks as the inclusion tag. No other exon combinations are generated; in human and mouse GENCODE this adds under 1% to the tag library. The added tags are listed in ``data/ME_synthetic_skip_tags.txt``.

.. _consecutive:

Consecutive microexons
======================

Some genes contain two or more microexons in a row, for example E1-m1-m2-E2. The junction that skips only m1 (E1 to m2) is m2's own inclusion junction, so reads on it align to m2's inclusion tag and were not counted as exclusion for m1. Likewise, a read on the m1-m2 junction fits both m1's and m2's inclusion tags, but Bowtie reports only one alignment, so one of the two microexons loses the read.

With ``consecutive_microexons : T`` (default), a read that crosses a junction of an inclusion tag also counts:

* as exclusion for every other microexon lying inside that junction, and
* as inclusion for the adjacent microexon on the other side of that junction (a full read if it spans that microexon, half a read otherwise).

Only the junction the read actually crosses is used, so runs of any length are handled. In runs of three or more, reads that reach a microexon two positions away are not credited to it, and its PSI can be slightly underestimated. ``Report/ME_consecutive_runs.txt`` lists every microexon that is adjacent to another one, with the members of its run.

.. _blacklist:

Reads that also fit the genome
==============================

After tag alignment, reads are aligned to the genome with Bowtie. A read that also aligns without splicing is likely genomic rather than a junction read, and is removed. With ``mismatch_aware_blacklist : T`` (default), a read is removed only if its genome alignment has no more mismatches than its tag alignment.

Some junctions have a copy elsewhere in the genome, typically in a processed pseudogene (a retrocopy of one isoform's mRNA). Reads from such a junction also fit the copy. When the copy differs from the junction, reads that cover the difference are kept; reads that fit both equally cannot be assigned and are still removed. A copy of the skipping junction therefore pushes PSI up, and a copy of the inclusion junction pushes it down. ``Report/ME_junction_genomic_copies.txt`` lists the microexons whose skipping or inclusion junction has an unspliced copy with up to 2 mismatches; treat their PSI with care.

Microexons whose position cannot be resolved
=============================================

A junction read shows the upstream exon, the microexon sequence and the downstream exon, but not where in the intron the microexon sits. A very short microexon (for example a 3 nt CAG) can match several positions flanked by canonical splice sites in the same intron. MicroExonator reports the position with the best U2 splice-site score, and ``Report/ME_ambiguous_positions.txt`` lists the alternative positions for each such event, so an expected coordinate that is missing from the output can be traced to the reported one.

Whippet and microexons
======================

Whippet quantifies splicing on a contiguous splice graph whose k-mers are long relative to microexons, so Whippet's own PSI for microexon nodes is not reliable. The ``.diff.ME.microexons`` results, which replace Whippet's PSI with MicroExonator's for microexon nodes, should be preferred for microexons. When the two disagree for a microexon, MicroExonator's value is the one to trust.
