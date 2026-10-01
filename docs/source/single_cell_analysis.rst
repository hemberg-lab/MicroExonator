.. single_cell_analysis


====================
Single cell analysis
====================

This section describes how to analyse single-cell RNA-seq data: how to quantify microexons across populations of cells and test alternative splicing between them. Single-cell experiments usually sequence each cell shallowly, so we developed a pseudo-pooling strategy to assess differential inclusion of microexons and other alternative splicing events between defined groups of cells (normally cell types defined beforehand from gene expression profiles). This module is called ``snakepool``.

.. note::

    This module uses Whippet. Before using it, follow the installation instructions in :doc:`differential_inclusion_analysis`.


Configuration
=============

To run this module, add the following parameters to ``config.yaml``:

.. code-block:: bash

    Single_Cell : T
    cluster_metadata : /path/to/cluster_metadata.tsv
    cluster_name : broad_type
    file_basename : Run_s
    cdf_t : 0.8
    min_p_mean : 0.9
    min_delta : 0.1
    min_rep : 25
    min_number_of_reads_single_cell : 5
    min_number_of_samples_single_cell : 3
    run_metadata : /path/to/run_metadata.tsv


* ``Single_Cell``: set to ``T`` to enable the single-cell modules.
* ``cluster_metadata``: path of a tab-separated file with a header row and at least two columns, giving the cluster and the sample name of each cell.
* ``cluster_name``: name of the column of ``cluster_metadata`` that holds the cluster names.
* ``file_basename``: name of the column of ``cluster_metadata`` that holds the sample names. These must match the sample names in the input files (see :doc:`setup`).
* ``cdf_t``: probability threshold, between 0.5 and 1. The probabilities of differential inclusion that a node gets across computational replicates are fitted with a `Beta distribution <https://en.wikipedia.org/wiki/Beta_distribution>`_, and its `cumulative distribution function <https://en.wikipedia.org/wiki/Cumulative_distribution_function>`_ at ``cdf_t`` gives a p-value for the node being above this threshold.
* ``min_p_mean``: minimum mean probability of differential inclusion across computational replicates for a node to be called differentially included.
* ``min_delta``: minimum mean ΔPSI across computational replicates for a node to be called differentially included.
* ``min_rep``: minimum number of computational replicates in which a node could be tested. Nodes with little coverage can only be tested in a few replicates, which makes their result unreliable. The number of replicates of each comparison is set in ``run_metadata``; we recommend setting ``min_rep`` to at least half of it.
* ``min_number_of_reads_single_cell`` and ``min_number_of_samples_single_cell``: minimum number of reads for a node to be used in a pseudo-bulk, and minimum number of pseudo-bulks per group in which it must be quantified. They are passed to ``whippet-delta`` as ``-r`` and ``-s``.
* ``run_metadata``: path of a tab-separated file describing the comparisons between cell types (see below).


Grouping for confidence filtering
---------------------------------

For a pure single-cell run, ``bulk_samples.tsv`` is not required. MicroExonator
chooses robustness-filter groups as follows:

* If ``cluster_metadata`` is supplied, every value in the configured
  ``cluster_name`` column becomes a sample group. All listed cells must match
  input sample names, and every input cell must occur exactly once. Spaces in
  cluster names are converted to underscores for output paths; labels that
  collide after conversion or contain other path-special characters are
  rejected.
* If ``cluster_metadata`` is omitted, all input cells are treated as one large
  sample group named ``all_cells``.

The cluster metadata fields are therefore optional for discovery and
quantification. They remain required when running cluster-dependent downstream
analyses such as ``snakepool``. A run containing both bulk and single-cell
inputs must provide both ``bulk_samples`` and ``cluster_metadata`` so every
input can be assigned unambiguously.


.. note::

    ``run_metadata`` is needed only for the ``snakepool`` target. The values of
    ``cdf_t``, ``min_p_mean``, ``min_delta``, ``min_rep`` and the two
    ``*_single_cell`` thresholds shown above are their defaults.


run_metadata
------------

``run_metadata`` is a tab-separated file with these columns:

.. list-table:: **run_metadata.tsv**
   :header-rows: 1

   * - Column
     - Description

   * - Compare_ID
     - Name of the comparison

   * - A.cluster_names
     - Comma-separated list of the cell types in group A

   * - A.number_of_pools
     - Number of pseudo-bulks to generate for group A

   * - B.cluster_names
     - Comma-separated list of the cell types in group B

   * - B.number_of_pools
     - Number of pseudo-bulks to generate for group B

   * - Repeat
     - Number of computational replicates: how many times cells are randomly assigned to pseudo-bulks and the comparison is run

.. warning::

    The order of the columns does not matter, but their names must be exactly as above. Other columns are ignored.

.. note::

    We recommend values of ``A.number_of_pools`` and ``B.number_of_pools`` that put at least five cells in each pseudo-bulk, and a ``Repeat`` of at least 10; more replicates give a better fit of the Beta distribution to the probabilities.

Optional configuration
----------------------

These parameters are optional:

.. code-block:: bash

    seed : 123
    Only_snakepool : T
    Get_Bamfiles : T

* ``seed``: seed for the random assignment of cells to pseudo-bulks (default 123). Keeping the same seed makes results reproducible and stops Snakemake from recomputing finished results.
* ``Only_snakepool``: with ``T``, discovery and quantification are skipped and only the splicing nodes of the annotation are tested. Useful when only alternative splicing events from the annotation are of interest.
* ``Get_Bamfiles``: with ``T``, BAM files are generated for visualisation (see below).

Run
===

With the files above in place, run the module with ``snakepool`` as the target:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} snakepool

.. note::

    Run a dry-run with ``-n`` before submitting a large set of jobs.

Pseudo-bulk quantification (optional)
-------------------------------------

Single cells are usually sequenced shallowly, so few reads cross any one splice junction in a cell. Pooling cells of the same type into pseudo-bulks increases the number of splicing nodes that can be quantified. With the ``collapse_pseudo_pools`` target, MicroExonator randomly assigns the cells of each cluster to pseudo-bulks, runs ``whippet-quant`` on each and writes the results to ``Whippet/Quant/Single_Cell/Pseudo_bulks/``, together with ``pseudo_bulk_membership.tsv``, which records which cells went into each pseudo-bulk:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} collapse_pseudo_pools

The number of pseudo-bulks per cluster is the number of cells divided by ``cells_pseudobulks`` (default 15), or a fixed number given by ``n_pseudobulks``. ``seed`` (default 123) fixes the random assignment. See :doc:`parameters` for all single-cell parameters.

Unpooled quantification (optional)
----------------------------------

To quantify each cell on its own instead of pseudo-bulks, use the ``quant_unpool_single_cell`` target. It writes one ``.psi.gz`` file per cell to ``Whippet/Quant/Single_Cell/Unpooled/``:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} quant_unpool_single_cell
    
This allows custom downstream analyses on the quantification of each cell. To avoid writing very many files, the per-cell results can instead be aggregated by cluster (using ``cluster_metadata``) with the ``collapse_whippet`` target; the results are written to ``Whippet/Quant/Collapsed/``.

.. warning::

    Only cells listed in ``cluster_metadata`` are processed.


Output
======

The ``whippet-delta`` results of every comparison in each computational replicate are in ``Whippet/Delta/Single_Cell/``. The combined results of each comparison are in ``Whippet/Delta/Single_Cell/Sig_nodes/``, with these columns:

.. list-table:: **all_nodes.microexons.txt**
   :header-rows: 1

   * - Column
     - Description

   * - Gene
     - Gene ID

   * - Node
     - Node number inside the gene

   * - Coord
     - Node coordinate

   * - Strand
     - Plus or minus strand

   * - Type
     - Node type. For more information visit `Whippet's GitHub page <https://github.com/timbitz/Whippet.jl#output-formats>`_.

   * - Psi_A.mean
     - Mean PSI of group A across computational replicates.

   * - Psi_B.mean
     - Mean PSI of group B across computational replicates.

   * - DeltaPsi.mean
     - Mean ΔPSI across computational replicates.

   * - DeltaPsi.sd
     - Standard deviation of ΔPSI across computational replicates.

   * - Probability.mean
     - Mean probability of differential inclusion across computational replicates.

   * - Probability.var
     - Variance of the probability across computational replicates.

   * - N.detected.reps
     - Number of replicates in which the node could be tested.

   * - cdf.beta
     - p-value of the probability being above ``cdf_t``

   * - is.diff
     - Whether the node is differentially included according to ``min_rep``, ``min_p_mean`` and ``min_delta``

   * - microexon_ID
     - Microexon ID (genomic coordinates)


Visualization
=============

To visualise the results, ``whippet-quant`` must write alignments. Add to ``config.yaml``:

.. code-block:: bash

    Get_Bamfiles : T

The alignments are converted to indexed BAM files with the ``cluster_bams`` target:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} cluster_bams

One BAM file is generated for each cell type in ``cluster_metadata``. With the coordinates of differentially included nodes and these BAM files, sashimi plots can be drawn with tools such as `ggsashimi <https://github.com/guigolab/ggsashimi>`_ or `IGV <http://software.broadinstitute.org/software/igv/>`_.
