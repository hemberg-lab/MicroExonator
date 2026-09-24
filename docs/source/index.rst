========
Overview
========

Microexons are exons of 30 nucleotides or less. They form one of the most conserved programs of alternative splicing described so far, are mostly included in neurons, and are misregulated in the brains of people with autism. They are also easy to miss: RNA-seq aligners that map reads to the genome often fail to place the very short anchor a microexon leaves in a read.

MicroExonator is a Snakemake workflow built for these events. It discovers novel microexons de novo from raw RNA-seq data of any organism with a genome and a gene annotation, quantifies annotated and novel microexons as percent spliced-in (PSI) values, and can test for differential inclusion between groups of samples or cell types with `Whippet <https://github.com/timbitz/Whippet.jl>`_ (`Sterne-Weiler et al. 2018 <https://doi.org/10.1016/j.molcel.2018.08.018>`_). It runs on a workstation or on a high-performance computing cluster, and installs its own software dependencies through conda.

The workflow is divided into modules, each of which can be run on its own:

* **Discovery** (Round1): reads are aligned to exon-exon junction tags to find putative novel microexons.
* **Quantification** (Round2): reads are aligned to tags that represent the inclusion and the skipping of every annotated and discovered microexon, and PSI values are computed per sample.
* **Differential inclusion**: microexon PSI values are compared between user-defined sample groups with Whippet.
* **Single-cell analysis**: microexons are quantified in single cells, in cell-type pseudo-bulks, and compared between cell types (Snakepool).

See :doc:`how_it_works` for a description of how reads are counted, and :doc:`whats_new` for the changes in the current version.

**Citation**

    Parada GE, Munita R, Georgakopoulos-Soares I, Fernandes HJR, Kedlian VR, Metzakopian E, Andres ME, Miska EA, Hemberg M (2021). MicroExonator enables systematic discovery and quantification of microexons across mouse embryonic development. *Genome Biology* 22:43. https://doi.org/10.1186/s13059-020-02246-2

    A step-by-step protocol is also available as a book chapter: Parada GE, Hemberg M. *Identification and quantification of microexons using bulk and single-cell RNA-seq data*, in the *Methods in Molecular Biology* series. The commands in that chapter correspond to the ``MiMB`` branch of the repository.

**Support**

    For questions, ideas, feature requests and bug reports, please open an issue on our `GitHub page <https://github.com/hemberg-lab/MicroExonator/issues>`_ or write to gp7@sanger.ac.uk.

.. toctree::
    :name: MicroExonator-install
    :maxdepth: 1
    :hidden:

    install

.. toctree::
    :name: MicroExonator-setup
    :maxdepth: 1
    :hidden:

    setup

.. toctree::
    :name: MicroExonator-discovery-and-quantification
    :maxdepth: 3
    :hidden:

    discovery_and_quantification

.. toctree::
    :name: MicroExonator-parameters
    :maxdepth: 2
    :hidden:

    parameters

.. toctree::
    :name: MicroExonator-how-it-works
    :maxdepth: 2
    :hidden:

    how_it_works

.. toctree::
    :name: MicroExonator-differential_inclusion_analysis
    :maxdepth: 3
    :hidden:

    differential_inclusion_analysis

.. toctree::
    :name: MicroExonator-single_cell_analysis
    :maxdepth: 3
    :hidden:

    single_cell_analysis

.. toctree::
    :name: MicroExonator-troubleshooting
    :maxdepth: 1
    :hidden:

    troubleshooting

.. toctree::
    :name: MicroExonator-whats-new
    :maxdepth: 1
    :hidden:

    whats_new

.. toctree::
    :name: MicroExonator-Licence
    :maxdepth: 1
    :hidden:

    licence

.. toctree::
    :name: MicroExonator-Support
    :maxdepth: 1
    :hidden:

    support
