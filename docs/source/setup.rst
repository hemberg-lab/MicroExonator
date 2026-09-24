.. _input_files:

===========
Setup
===========

Before running MicroExonator, a few files need to be created inside the ``MicroExonator/`` folder: one or more lists of input samples, a ``config.yaml`` file and, on a cluster, a ``cluster.json`` file. Examples of all of them are in the ``Examples/`` folder.

RNA-seq samples
===============

Samples can be stored locally, downloaded from a URL, or downloaded from NCBI's Sequence Read Archive (SRA). Each source has its own file, the files must sit in the ``MicroExonator/`` folder under these exact names, and at least one of them must exist. Sources can be combined in one run.

**Local files: local_samples.tsv**

A tab-separated file with a header and the columns ``path`` and ``sample``. ``path`` points to a FASTQ file (``.fastq.gz``, ``.fastq.bz2`` or plain ``.fastq``) and ``sample`` is the name the sample will have in every output:

.. code-block:: text

    path	sample
    /data/condition_A.rep1.fastq.gz	A1
    /data/condition_A.rep2.fastq.gz	A2
    /data/condition_A.rep3.fastq.gz	A3
    /data/condition_B.rep1.fastq.gz	B1
    /data/condition_B.rep2.fastq.gz	B2
    /data/condition_B.rep3.fastq.gz	B3

**Files at a URL: sample_url.tsv**

The same structure, with a column named ``url`` instead of ``path``. Each URL must point directly to a gzip-compressed FASTQ file.

**SRA runs: NCBI_accession_list.txt**

One SRA run accession per line, with no header:

.. code-block:: text

    SRR1805814
    SRR1805815
    SRR1805816

The run accession becomes the sample name. For paired-end runs both mates are downloaded and joined into one FASTQ file, and each mate is then treated as a separate read. To report paired-end data from local files as a single sample, see ``paired_samples`` in :doc:`parameters`.

By default, downloaded or copied FASTQ files are kept only while they are needed (see ``Keep_fastq_gz`` and ``Optimize_hard_drive`` in :doc:`parameters`).

Sample groups
=============

Bulk RNA-seq runs also need a ``bulk_samples`` file that assigns every sample to a biological condition. It is used to decide which microexons are detected robustly within each group (see :doc:`discovery_and_quantification`). It is a tab-separated file with the columns ``sample`` and ``condition``:

.. code-block:: text

    sample	condition
    A1	control
    A2	control
    A3	control
    B1	treated
    B2	treated
    B3	treated

Sample names must match the names given in the sample lists above. Its path is given with the ``bulk_samples`` key in ``config.yaml``.

Configuration file
==================

All parameters are given in a file called ``config.yaml`` inside the ``MicroExonator/`` folder, in `YAML <https://yaml.org/>`_ format (``key : value``, one per line). The minimal configuration is described in :doc:`discovery_and_quantification` and every parameter is listed in :doc:`parameters`. Complete examples are in ``Examples/Runs/``.

Cluster configuration
=====================

On a cluster, jobs are usually submitted through a scheduler such as LSF or SLURM. Snakemake can submit each step as a separate job; it reads the resources for each step from a JSON file, which we assume is called ``cluster.json`` and sits in the ``MicroExonator/`` folder.

Every step of the workflow is a Snakemake rule (in the ``rules/`` folder). Default resources for all rules are given under ``"__default__"``; for LSF this may look like:

.. code-block:: json

    "__default__" :
    {
        "queue"     : "normal",
        "nCPUs"     : "1",
        "memory"    : 10000,
        "resources" : "\"select[mem>10000] rusage[mem=10000] span[hosts=1]\"",
        "name"      : "JOBNAME.{rule}.{wildcards}",
        "output"    : "logs/{rule}.{wildcards}.out",
        "error"     : "logs/{rule}.{wildcards}.err",
        "Group"     : "your_group",
        "tCPU"      : "99999"
    }

This gives every job one CPU and 10 GB of memory, which is enough for most steps with the human or mouse genome. Rules that need more are given their own entry, which overrides the default. For example, to give the discovery alignment five CPUs:

.. code-block:: json

    "Round1_bwa_mem_to_tags" :
    {
        "nCPUs"    : 5
    }

A complete file for LSF is in ``Examples/Cluster_config/lsf/cluster.json``. For other schedulers, see the `Snakemake documentation <https://snakemake.readthedocs.io/en/v7.32.4/snakefiles/configuration.html#cluster-configuration-deprecated>`_.

.. _cluster_resources:

Recommended resources
---------------------

.. list-table::
   :header-rows: 1

   * - Rule
     - CPUs
     - Memory (GB)
   * - Round1_bwa_mem_to_tags
     - 5
     - default
   * - hisat2_genome_index
     - 5
     - default
   * - Round2_bowtie_to_tags
     - 5
     - default
   * - bowtie_genome_index
     - default
     - 30 for the human genome
   * - bowtie_to_genome
     - 2
     - default
   * - total_hisat2_to_genome
     - 5
     - default
   * - Output
     - 2
     - 30
   * - whippet_quant
     - default
     - 2
