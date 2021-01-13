.. single_cell_analysis


====================
Single cell analysis
====================

On this secction we describe how to perform analysis of single-cell RNA-seq data to quantify microexons across populations of cells and alternative splicing events across them. Since single cell experiments usually provides a shallow sequencing depth for each cell, we have developed a pseudo-pooling strategy to assess differential inclusion of microexons and other types of alternative splicing events across defined groups of cells (normally corresponding to cell-types previously defined by gene expresion profiles). We named this single cell analysis module ``snakepool`` an on this secction we describe how to use it.

.. note::

    Before using this module you must follow the same installation instructions decrived in :doc:`differential_inclusion_analysis` secction.


Configuration
=============

To run this secction the following parameters needs to be incorporated at ``config.yaml``.

.. code-block:: bash

    Single_Cell : T
    cluster_metadata : /lustre/scratch117/cellgen/team218/gp7/Micro-exons/Runs/Paper/MicroExonator/Whippet/Tasic_clustering.txt
    cluster_name : broad_type
    file_basename : Run_s
    cdf_t : 0.8
    min_p_mean : 0.9
    min_delta : 0.1
    min_rep : 25
    run_metadata : /lustre/scratch117/cellgen/team218/gp7/Micro-exons/Runs/Paper/MicroExonator/Whippet/Tasic_run.txt


* ``Single_Cell`` correspond to an optiona parameter that needs to be set as ``T`` in order to run ``snakepool``.
* ``cluster_metadata`` must indicate the path of a tabular separated file that contain at least two colums to indicate the cluster and file base names. This file must have the first row as header.
* ``cluster_name`` indicate the name of the column, inside ``cluster_metadata``, which has cluster name information.
* ``file_baseame`` indicate the name of the column, inside ``cluster_metadata``, which has the sample names. These needs to match with sample names defined on the input files (See :doc:`setup`)
* ``cdf_t`` parameter set a theshold to run a `Cumulative distribution function <https://en.wikipedia.org/wiki/Cumulative_distribution_function>`_ over the resultant Probability values obtained for each node across comutational replicates, asuming these fit a `beta distribution <https://en.wikipedia.org/wiki/Beta_distribution>`_. Which in practical terms can be considered as a user-defined threshold (between 0.5 and 1) to calculate a p-value assosiated to node's probability of differential inclusion across computational replicates.
* ``min_p_mean`` corresponds to a threshold of mean probability values across computational samples to define a node as differentially included across the comparing cell-types.
* ``min_p_delta`` corresponds to a threshold of mean delta PSI values across computational samples to define a node as differentially included across the comparing cell-types.
* ``min_rep`` minimun set of computational replicates that can be considered to define a node as differentially included across the comparing cell-types. This parameter is relevant because when nodes have limmited read coverage across cells, only some few computational replictates could enable quantitative alternative splicing analyses, leading to unreliable assesment of its alternative inclusion. The number of compuational repeats is defined for each sample inside `run_metadata`. We recomend to set this value to at least the half of the compuational replicates that are scheduled to run by the user.
* ``run_metadata`` indicates the path of a tabulat separated file which contain information about user-defined comparisons across cell-types. Additional information about this file can be found bellow.


.. warning::

    If ``Single_Cell`` is set to ``T`` all the other parameters listed above will become compulsory. Thus, you should only activate this module if you have the required parameters on ``config.yaml``.  


run_metadata
------------

The file indicated by ``run_metadata`` must be a tabular separated file containing the following columns:

.. list-table:: **run_metadata.tsv**
   :header-rows: 1

   * - Column
     - Description

   * - Compare_ID
     - User-defined name for scheduled comparions

   * - A.cluster_names
     - Comma-separated list of cell-types to be concider as part of `sample group A`

   * - A.number_of_pools
     - Number of pseudo-bulk to be generated for `sample group A`

   * - B.cluster_names
     - Comma-separated list of cell-types to be concider as part of `sample group B`

   * - B.number_of_pools
     - Number of pseudo-bulk to be generated for `sample group B`

   * - Repeat
     - Node type. For more information visit `Whippet's GitHub page <https://github.com/timbitz/Whippet.jl#output-formats>`_.

.. warning::

    The order of the columns is not relevant, however the column names must correspond to the ones indicated above. Additional columns with other names will be ignored.

.. note::

    We recomend to consider ``A.number_of_pools`` and ``B.number_of_pools`` values that ensure that at least five cells are merged within each pseudo-bulk pool. We also suggest 10 as a minimun value of ``Repeat``, higher values will enable better estimation of parameters to fit the resultant probabilities into beta distribution models. 

Optional Configuration
----------------------

The following parameter are optionals to be defined inside ``config.yaml`` file:

.. code-block:: bash

    seed : 123
    Only_snakepool : T
    Get_Bamfiles : T

* ``seed`` define a specific seed for pseudo number geration. This number influence the arrangement cells into the corresponding pseudo-bulks. Mataining the same seed ensures reproducibility of the results and prevent snakemake of overwrite completed results.
* ``Only_snakepool`` is a bolean variable that if its defined as ``T`` it will force MicroExonator to skip Disovery and Quantification modules. This mode is useful for users who are only interested to find alterantive splicing events from splicing nodes that can be extracted from the annotation.
* ``Get_Bamfiles`` correspond to a bolean variable that if its defined as ``T`` enable the generation of BAM files that can be used for visualization purposes.

Run
===

After setting up all the files described above, this single cell analysis module can be run by adding `snakepool` as target for snakemake:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} snakepool

.. note::

    It is allways a good idea to use ``-np`` to execute an snakemake ``dry-run`` before submiting a large set of jobs.

Unpooled quantification (optional)
----------------------------------

In order to generate PSI quantification files at the single cell level (as opposed to pseudo-bulks), you can run MicroExonator with ``quant_unpool_single_cell`` as a target for snakemake. By doing this `.psi.gz` files will be generated at ``Whippet/Quant/Single_Cell/Unpooled/`` folder:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} quant_unpool_single_cell
    
This can enable users do run custom downstream analysis over alternative splicing quantification files generated for every cell by separated. 

.. warning::

    Only FASTQ files from cells annotated on ``cluster_metadata`` file will be processed.


Output
======

Direct results from `whippet-delta` for every comparison across each computational replicate can be found at `Whippet/Delta/Single_Cell/`. Integrated results for each comparion can be found at ``Whippet/Delta/Single_Cell/Sig_nodes``, these resuls are structured as follow:

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
     - Mean PSI values for group ``cluster A`` across computational replicates.

   * - Psi_B.mean
     - Mean PSI values for group ``cluster B`` across computational replicates.

   * - DeltaPsi.mean
     - Mean DeltaPsi values obtained across computationa replicates.

   * - DeltaPsi.sd
     - Standar deviation of DeltaPsi values obtained across computationa replicates.

   * - Probability.mean
     - Mean probability of differential inclusion obtained across computationa replicates.

   * - Probability.var
     - Variance of probability across computationa replicates.

   * - N.detected.reps
     - Number of replicates in which the differential inclusion could be assessed.

   * - cdf.beta
     - p-value of being above the used defined probability threshold ``cdf_t``

   * - is.diff
     - Bolean variable defining wheather the node was differentially included accoding to the user-defined criteria (``min_rep``, ``min_p_mean`` and ``min_delta``)

   * - microexon_ID
     - Microexon ID based on its genomic coodinates.


Visualization
=============

In order to visualize the results, you need to instruct ``whippet-quant`` to generate SAM files by incorporating the following parameter with a ``True`` value inside ``config.yaml``: 

.. code-block:: bash

    Get_Bamfiles : T

Samfiles are further converted to BAM files and corresponding index files are genrated to enable their visualization. To generate these BAMs ``cluster_bams`` needs to be defined as an snakemake target:

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} cluster_bams

A BAM file will be generated for every cell-type defined at ``cluster_metadata`` file. Given the coodinates of differentially included splicing nodes and the correspondig BAM files, sashimi plots can be generated by using tools such as `ggsashimi <https://github.com/guigolab/ggsashimi>`_ or `IGV <http://software.broadinstitute.org/software/igv/>`_.
