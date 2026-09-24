.. _troubleshooting:

===============
Troubleshooting
===============

Read the error of the failed step
=================================

When a step fails, Snakemake prints the rule and the command that failed. On a cluster, the full error is usually in the scheduler's output files. With the ``cluster.json`` example in :doc:`setup`, these are written to ``logs/{rule}.{wildcards}.out`` and ``.err``, named after the rule and the sample, next to the log files MicroExonator writes itself. Snakemake's own log of each run is in ``.snakemake/log/``.

Many cluster failures are caused by too little memory for a step. Give that rule more memory in ``cluster.json`` (see :ref:`cluster_resources`) and submit the same command again.

Resume a run
============

Submitting the same command again resumes a stopped or failed run: only jobs whose output is missing are run. With ``-k``, independent jobs keep running while one fails, so fixing a failed step does not require starting over.

Protected files
===============

Some final outputs, such as the per-sample discovery results (``Round1/*.sam.row_ME.filter1``), are protected so that they cannot be overwritten or deleted by accident. A run that would replace them stops with an error. If you really want to regenerate them, remove the files (or make them writable) yourself first.

Start again with different settings
====================================

To rerun with a different genome, annotation or parameters that change the tags, start from a fresh clone of the repository, or delete the folders created by the previous run: ``data/``, ``Round1/``, ``Round2/``, ``Report/``, ``Whippet/`` and ``download/``.

The mixture-model filter fails
==============================

The ``legacy_mixture`` confidence filter fits a Gaussian mixture to the U2 scores of the quantified microexons. With few microexons the fit can fail, and the pipeline stops with a message recommending ``filter_method: robustness``, which is the default and does not need a model fit. See :doc:`discovery_and_quantification`.

Unexpected read lengths
=======================

``Report/read_lengths.tsv`` lists read lengths per sample and warns about samples whose median read is shorter than 40 nt or whose reads were mostly longer than ``max_read_len``. Public datasets sometimes include single-cell or other non-standard libraries under a bulk RNA-seq label; these warnings help find them. If most reads of ordinary bulk samples are being trimmed, raise ``max_read_len`` to their read length (see :ref:`read_length`).

Snakemake version
=================

The current version is tested with Snakemake 7.32.4. Snakemake 8 changed how cluster submission is configured (``--cluster`` and ``--cluster-config`` were replaced by executor plugins). If you see errors about these options, create the Snakemake environment with ``snakemake=7.32.4`` as described in :doc:`install`.
