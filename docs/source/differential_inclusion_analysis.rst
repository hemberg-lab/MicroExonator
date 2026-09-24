.. differential_inclusion_analysis


===============================
Differential inclusion analysis
===============================

This downstream module tests for differential inclusion of annotated and novel microexons between user-defined groups of samples, using `Whippet <https://github.com/timbitz/Whippet.jl>`_. Whippet builds a contiguous splice graph for each gene from the annotation, where nodes are non-overlapping stretches of exonic sequence and edges are splice junctions or contiguous sequence, and quantifies splicing by mapping reads directly to these graphs. MicroExonator adds the detected microexons to the annotation first, so that they become nodes of the graph.

Whippet's k-mers are long relative to microexons, so its own PSI for microexon nodes is not reliable. For microexons, use the ``.diff.ME.microexons`` results described below, which are computed from MicroExonator's PSI.

Install
=======

Whippet runs on Julia. The current Whippet release works with Julia 1.6, which can be built from its repository:

.. code-block:: bash

    git clone https://github.com/JuliaLang/julia
    cd julia
    git checkout v1.6.0
    make

This produces a ``julia`` executable inside the ``julia/`` folder (prebuilt binaries from `julialang.org <https://julialang.org/downloads/oldreleases/>`_ work too). Then install Whippet with that Julia:

.. code-block:: bash

    git clone https://github.com/timbitz/Whippet.jl.git
    cd Whippet.jl
    /path/to/julia --project -e 'using Pkg; Pkg.instantiate()'

where ``/path/to/julia`` is the full path of the Julia executable.

Configure
=========

Add these keys to ``config.yaml``:

.. code-block:: yaml

    whippet_bin_folder : /path/to/Whippet.jl/bin
    julia : /path/to/julia/julia
    Gene_anontation_GTF : /path/to/gencode.annotation.gtf
    whippet_delta : /path/to/whippet_delta.yaml

* ``whippet_bin_folder``: the ``bin/`` folder of the Whippet installation.
* ``julia``: the Julia executable used to install Whippet.
* ``Gene_anontation_GTF``: the gene annotation as `GTF <https://en.wikipedia.org/wiki/Gene_transfer_format>`_, with gene, transcript and exon features. Use the same annotation as for discovery and quantification.
* ``whippet_delta``: a YAML file describing the comparisons (see below).

MicroExonator passes its corrected PSI values to Whippet by default. Set ``use_uncorrected_PSI : T`` only to reproduce analyses that used the uncorrected values; it affects the ``.diff.ME.microexons`` results, not the PSI values Whippet calculates itself.

whippet_delta YAML file
-----------------------

Each comparison has this structure:

.. code-block:: bash

    comparison_ID:
      A : sample1,sample2,sample3
      B : sample4,sample5,sample6

where ``sample1`` to ``sample6`` are sample names as defined in the input files (see :doc:`setup`; for SRA data, the run accessions) and ``comparison_ID`` is any name for the comparison. There is no limit on the number of comparisons or on group sizes. The :download:`YAML file <../../Examples/Runs/Parada_et_al/whippet_delta.yaml>` used in our publication is an example.

.. warning::

    Inside this YAML file sample groups must be named ``A`` and ``B``.


Optional parameters
-------------------

To skip discovery and quantification and analyse only the splicing events in the GTF annotation, add:

.. code-block:: bash

    downstream_only : T

Run
===

Run the usual MicroExonator command with ``differential_inclusion`` as the target. If discovery and quantification have not been run yet, they are added to the job plan and their results are passed to Whippet automatically (unless ``downstream_only`` is ``T``). A single command can therefore run thousands of jobs; a dry-run (``-n``) first is recommended.

.. code-block:: bash

    snakemake -s MicroExonator.smk  --cluster-config cluster.json --cluster {cluster system params} --use-conda -k  -j {number of parallel jobs} differential_inclusion



Output
======

Results are in the ``Whippet/`` folder.

* ``Whippet/Quant/``: Whippet quantification of every sample. ``.psi.gz`` files contain PSI for all splicing nodes of the annotation (see the `Whippet documentation <https://github.com/timbitz/Whippet.jl#output-formats>`_), and ``.psi.ME.gz`` files the same table with MicroExonator's PSI for microexon nodes.
* ``Whippet/Delta/``: differential inclusion for each comparison, in Whippet's ``.diff.gz`` format. For easier reading, the nodes that correspond to microexons are also written to ``.microexons`` files: ``.diff.ME.microexons`` uses MicroExonator's PSI and should be preferred for microexons; ``.diff.microexons`` uses Whippet's own PSI. The ``.diff.ME.microexons`` files are not produced when ``downstream_only`` is ``T``.
