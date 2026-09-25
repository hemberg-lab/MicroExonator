.. discovery_and_quantification

============================
Discovery and Quantification
============================

Discovery (Round1) and quantification (Round2) are the core modules of MicroExonator. A normal run executes both, but each can also be run on its own. During discovery, reads are aligned with BWA-MEM to exon-exon junction tags to find putative novel microexons. During quantification, reads are aligned with Bowtie to tags representing the inclusion and the skipping of every annotated, database and discovered microexon, and PSI values are computed per sample. :doc:`how_it_works` describes both steps in detail.

The code for each module is in the ``rules/`` folder, and its intermediate files are written to ``Round1/`` and ``Round2/``. Snakemake can list every command before it runs (a dry-run) and draw the graph of all jobs, which makes each step easy to inspect.

Configuration
=============

Minimal configuration
---------------------

These parameters are needed for every run, even when only downstream modules are used. A minimal ``config.yaml`` looks like this:

.. code-block:: yaml

    Genome_fasta : /path/to/Genome.fa
    Gene_anontation_bed12 : /path/to/annotation.bed12
    GT_AG_U2_5 : /path/to/GT_AG_U2_5.good.matrix
    GT_AG_U2_3 : /path/to/GT_AG_U2_3.good.matrix
    conservation_bigwig : /path/to/conservation.bw
    working_directory : /path/to/MicroExonator/
    ME_len : 30
    max_read_len : 150
    min_number_files_detected : 3
    bulk_samples : /path/to/bulk_samples.tsv

* ``Genome_fasta``: a multi-FASTA file with the chromosomes.
* ``Gene_anontation_bed12`` (or ``Gene_anontation_GTF``): the transcript annotation. BED12 files can be exported from the `UCSC Table Browser <http://genome.ucsc.edu/cgi-bin/hgTables>`_. The genome and the annotation must come from the same source, so that every annotated chromosome exists in the genome.
* ``GT_AG_U2_5`` and ``GT_AG_U2_3``: splice-site position weight matrices, for example from `SpliceRack <http://katahdin.cshl.edu/splice/>`_. Human and mouse examples are in the ``PWM/`` folder. With ``NA``, MicroExonator builds them from annotated splice sites.
* ``conservation_bigwig``: genome-wide phyloP or PhastCons scores, available from the `UCSC Genome Browser <http://hgdownload.cse.ucsc.edu/downloads.html>`_ for many species, or ``NA``.
* ``working_directory``: the path of the cloned ``MicroExonator/`` folder.
* ``ME_len``: maximum microexon length (default 30).
* ``max_read_len``: tag flank length and read trimming length (default 150). Set it to the length of the longest reads in your data; see :ref:`read_length`.
* ``min_number_files_detected``: the number of samples in which a putative novel microexon must be found during discovery. At least 2 is recommended for single-end data and 3 when paired-end files are present.
* ``bulk_samples``: a tab-separated file assigning every bulk sample to a condition (columns ``sample`` and ``condition``; see :doc:`setup`). Sample names must match the input names exactly. For paired-end data, list the first-read sample declared in ``paired_samples``. Condition names are used in output paths: spaces become underscores, and other characters that are not valid in paths are rejected. Pure single-cell runs do not need this file (see :doc:`single_cell_analysis`).

Optional configuration
----------------------

Every other parameter is listed, with its default, in :doc:`parameters`. The ones most often changed for discovery and quantification are:

* ``ME_DB``: a file of known microexons to quantify in addition to the annotated and discovered ones, for example from VastDB. Each line is read on its own, so one file can mix BED12 rows (internal blocks of at most ``ME_len`` nt with canonical AG/GT splice sites), BED6 rows (same filters) and single-column ``chrom_strand_start_end`` IDs (used as given). Events that do not lie inside an annotated intron cannot be quantified and are listed in ``data/DB.ME_centric.non_overlaping.txt``.
* ``min_reads_PSI``: minimum number of junction reads to report a PSI value (default 5).
* ``skip_discovery : T`` together with ``only_db : T``: quantify only the annotated microexons and those in ``ME_DB``, without discovery.
* ``paired_samples``: pool the two files of each paired-end sample.
* ``synthetic_skip_tags``, ``consecutive_microexons``, ``shared_splice_site_correction`` and ``mismatch_aware_blacklist`` (all ``T`` by default): corrections for microexons without an annotated skipping junction, for adjacent microexons, for microexons that share a splice site with a longer exon, and for junctions with a copy in the genome. They are explained in :doc:`how_it_works`. Set them to ``F`` to reproduce results from earlier versions.

.. note::

    VastDB microexons can be exported as BED12 by adding the `VastDB track hub <http://vastdb.crg.eu/tracks/VastDBhub/hub.txt>`_ to the UCSC Genome Browser and downloading the track with the `Table Browser <http://genome-euro.ucsc.edu/cgi-bin/hgTables>`_.

Confidence filtering
--------------------

``filter_method`` selects the filter that decides which microexons are reported as detected. Its default value is ``robustness``:

.. code-block:: yaml

    filter_method : robustness
    min_detected_samples : 1

The robustness filter keeps a microexon within a sample group when the lower bound of its PSI confidence interval is at least 0.1 in more than half of its measurements and it has at least three reads spanning the whole microexon across that group. A microexon must meet these criteria in at least ``min_detected_samples`` measurements within one group. The result is written to ``Report/out.robustly_detected.txt``.

The Gaussian mixture filter of the original publication is kept for reproducibility, but it is no longer the default and must be requested explicitly:

.. code-block:: yaml

    filter_method : legacy_mixture

This route writes ``Report/out.high_quality.txt`` and keeps the historical posterior calculation unchanged. The description of which mixture component is used in the published Methods does not fully match the implementation; the implementation is kept as it was so that old analyses can be reproduced. This is tracked in `Issue #37 <https://github.com/hemberg-lab/MicroExonator/issues/37>`_.

The mixture fit needs enough quantified microexons with variation in their U2 scores. When there are too few, the pipeline stops with a message that recommends ``filter_method: robustness``; it does not switch methods silently. The old ``filter_mode`` and ``skip_mixture_model_filter`` settings are still accepted, with a deprecation warning.

Run
===

A run can take hours or days, depending on the amount of data and the hardware, so on a remote machine start it inside a ``screen`` (or ``tmux``) session, which keeps running if the connection is lost:

.. code-block:: bash

    screen -S session_name

Activate the environment with Snakemake:

.. code-block:: bash

    conda activate snakemake_env

and run MicroExonator from inside the ``MicroExonator/`` folder:

.. code-block:: bash

    snakemake -s MicroExonator.smk --cluster-config cluster.json --cluster "{cluster system params}" --use-conda -k -j {number of parallel jobs}

* ``-s MicroExonator.smk`` runs the main workflow file.
* ``--cluster-config`` and ``--cluster`` tell Snakemake how to submit jobs to the scheduler. ``{cluster system params}`` is the submission command, with resources read from ``cluster.json`` through ``{cluster.*}`` fields. On a workstation, leave both flags out.
* ``--use-conda`` creates the conda environments each step needs.
* ``-j`` is the maximum number of jobs running at the same time; between 5 and 50 suits most users.

For example, on LSF with up to 50 simultaneous jobs:

.. code-block:: bash

    snakemake -s MicroExonator.smk --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}" --use-conda -k -j 50

To detach from the screen session press ``Ctrl-a d``; ``screen -ls`` lists sessions and ``screen -r session_name`` reattaches one.

Useful Snakemake flags
----------------------

* ``-n`` (dry-run): build the job plan and check the configuration without running anything. Always worth doing before a large run. If it fails, check that you are inside the ``MicroExonator/`` folder and that ``config.yaml`` is valid.
* ``-p``: print the shell command of every job.
* ``-k``: keep running independent jobs when one fails. Failed steps can be fixed and the run resumed; finished files are not recomputed.
* ``--notemp``: keep the intermediate files that are normally deleted. This uses a lot of disk space, but helps with troubleshooting.
* ``--resources get_data=N``: allow at most N downloads or copies of input files at the same time.

Targets
-------

A target at the end of the command selects what to produce. Without one, MicroExonator uses ``quant``, which runs discovery and quantification.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Target
     - Result
   * - ``quant`` (default)
     - Discovery and quantification.
   * - ``discovery``
     - Discovery only. Produces ``Round1/*.sam.row_ME.filter1`` files with the putative microexons of each sample. These files are protected: Snakemake will not overwrite or delete them.
   * - ``Report/out.robustly_detected.gtf``
     - The annotation with the detected microexons added, for use in other splicing tools such as DEXSeq, rMATS or Whippet (``Report/out.high_quality.gtf`` with ``filter_method: legacy_mixture``).
   * - ``differential_inclusion``
     - Differential inclusion between sample groups with Whippet (see :doc:`differential_inclusion_analysis`).
   * - ``snakepool``, ``quant_unpool_single_cell``, ``collapse_pseudo_pools``, ``cluster_bams``
     - Single-cell modules (see :doc:`single_cell_analysis`).

Running large datasets
----------------------

Because MicroExonator is a Snakemake workflow, it scales to large datasets. A few things help when processing many samples.

**Limit simultaneous downloads.** Downloads from NCBI can fail when too many run at once. Limit them with the ``get_data`` resource, for example to 50:

.. code-block:: bash

    snakemake -s MicroExonator.smk --cluster-config cluster.json --cluster "{cluster system params}" --use-conda -k -j {number of parallel jobs} --resources get_data=50

**Save disk space.** If there is not enough space to hold all FASTQ files at once, set ``Optimize_hard_drive : T`` and run the two modules one after the other. First run discovery:

.. code-block:: bash

    snakemake -s MicroExonator.smk --cluster-config cluster.json --cluster "{cluster system params}" --use-conda -k -j {number of parallel jobs} discovery

Each FASTQ file is deleted as soon as discovery has used it. Then run quantification:

.. code-block:: bash

    snakemake -s MicroExonator.smk --cluster-config cluster.json --cluster "{cluster system params}" --use-conda -k -j {number of parallel jobs} quant

The files are downloaded or copied again for quantification and deleted after use.

**Mind the number of files.** Clusters often limit the number of files per user as well as the space. After a run has finished, the ``.snakemake/`` and ``logs/`` folders can be deleted:

.. code-block:: bash

    rm -rf .snakemake logs

**Check read lengths.** When combining public datasets, check ``Report/read_lengths.tsv`` for samples with unusual read lengths, which often turn out to be another library type (for example single-cell data).

.. note::

    If a run stops or fails, submit the same command again. Snakemake reruns only the jobs whose output is missing. To process every sample again from the start, delete the ``download/`` folder as well.

Output
======

The main results are in the ``Report/`` folder.

Detected microexons
-------------------

With the default robustness filter, detected microexons are listed in ``Report/out.robustly_detected.txt``:

.. list-table:: **out.robustly_detected.txt**
   :header-rows: 1

   * - Column
     - Description
   * - ME
     - Microexon coordinates (``chrom_strand_start_end``)
   * - Transcript
     - Transcript in which the microexon was found
   * - Total_coverage
     - Total coverage across all splice junctions of the microexon
   * - Total_SJs
     - Splice junctions (introns) in which the microexon was found
   * - ME_coverages
     - Comma-separated coverage of each of those junctions
   * - ME_length
     - Microexon length
   * - ME_seq
     - Microexon sequence
   * - ME_matches
     - Number of matches of the microexon sequence inside the intron
   * - U2_score
     - U2 splice-site score
   * - Mean_conservation
     - Mean conservation score (if ``conservation_bigwig`` was given)
   * - P_MEs
     - Microexon confidence score
   * - Total_ME
     - Coordinates, U2 score and conservation of every position that matches the microexon sequence

The legacy mixture filter adds the columns ``ME_P_value`` (value used by the final filters) and ``ME_type`` (``IN``, ``RESCUED`` or ``OUT``) to ``Report/out.high_quality.txt``, and also writes ``out_shorter_than_3_ME.txt`` and ``out_low_scored_ME.txt`` for microexons that are likely false positives.

Quantification
--------------

PSI values per sample are in ``Report/quant/``. The files in ``Report/quant/corrected/PSI_sparse/`` use corrected counts, in which reads spanning the whole microexon count as one and reads crossing only one of its junctions count as half (see :doc:`how_it_works`). These are the values used by the confidence filter and by the Whippet modules. The per-sample quantification tables have these columns:

.. list-table:: **Quantification output**
   :header-rows: 1

   * - Column
     - Description
   * - File
     - Sample name
   * - ME_coords
     - Microexon coordinates
   * - SJ_coords
     - Splice junctions whose reads were counted for this microexon
   * - ME_coverages
     - Comma-separated number of reads supporting inclusion, per splice junction
   * - SJ_coverages
     - Comma-separated number of reads supporting exclusion (skipping), per splice junction
   * - PSI
     - Percent spliced-in
   * - CI_Lo
     - Lower bound of the 95% confidence interval of PSI
   * - CI_Hi
     - Upper bound of the 95% confidence interval of PSI
   * - Alt5
     - Alternative 5' splice sites of the microexon
   * - Alt3
     - Alternative 3' splice sites of the microexon
   * - Alt5_coverages
     - Reads supporting the alternative 5' splice sites
   * - Alt3_coverages
     - Reads supporting the alternative 3' splice sites
   * - Unique_ME_reads
     - Number of distinct read sequences supporting inclusion
   * - sum_ME_SJ_coverage_up
     - Reads covering the upstream splice junction of the microexon
   * - sum_ME_SJ_coverage_down
     - Reads covering the downstream splice junction of the microexon

Diagnostic reports
------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File
     - Content
   * - ``Report/read_lengths.tsv``
     - Reads, median and maximum read length, and percentage of reads trimmed to ``max_read_len`` for each sample, with a note on samples that look like another library type.
   * - ``Report/ME_consecutive_runs.txt``
     - Microexons adjacent to another microexon, with the members of their run. In runs of three or more, PSI can be slightly underestimated.
   * - ``Report/ME_junction_genomic_copies.txt``
     - Microexons whose skipping or inclusion junction has an unspliced copy (up to 2 mismatches) elsewhere in the genome. A skipping copy biases PSI upwards, an inclusion copy downwards.
   * - ``Report/ME_ambiguous_positions.txt``
     - Short microexons whose sequence fits several positions in the same intron, with the alternative positions.
   * - ``data/ME_synthetic_skip_tags.txt``
     - Skipping tags added for introns without an annotated skipping junction.
