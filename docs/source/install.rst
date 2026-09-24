.. _Installation:

=====================
Installation
=====================

Hardware
========

MicroExonator runs on any Unix-like system (Linux or macOS), either on a standalone workstation or on a high-performance computing (HPC) cluster. A cluster is recommended for anything beyond a handful of samples. A single machine needs at least 16 GB of RAM and 4 CPUs; building the Bowtie index of the human genome needs about 30 GB. Recommended resources for individual steps are listed in :ref:`cluster_resources`.

Conda and Snakemake
===================

All software dependencies are handled by conda. Snakemake creates the environments each step needs the first time it runs, so the only thing you install yourself is conda and an environment that contains Snakemake.

Install Miniconda (or any conda distribution) following the `conda documentation <https://docs.conda.io/en/latest/miniconda.html>`_, for example:

.. code-block:: bash

   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh

Open a new terminal, then add the channels that provide bioinformatics software:

.. code-block:: bash

   conda config --add channels defaults
   conda config --add channels bioconda
   conda config --add channels conda-forge

Install ``mamba``, which resolves environments much faster than conda, and use it to create an environment with Snakemake:

.. code-block:: bash

   conda install -n base -c conda-forge mamba
   mamba create -n snakemake_env -c conda-forge -c bioconda snakemake=7.32.4

The current version of MicroExonator is tested with Snakemake 7.32.4. Snakemake 8 changed how cluster execution is configured, so pinning version 7 is the safest choice.

Get MicroExonator
=================

Clone the repository:

.. code-block:: bash

  git clone https://github.com/hemberg-lab/MicroExonator

This creates a folder called ``MicroExonator``. The workflow runs inside this folder and writes all of its intermediate files and results there, so we recommend a fresh clone for every project or run. This keeps intermediate files from different runs apart, avoids overwriting finished results, and keeps the configuration files of each run together with its output.

.. note::

    The commands in the *Methods in Molecular Biology* chapter correspond to the ``MiMB`` branch, which is kept unchanged for reproducibility. To follow the chapter exactly, run ``git checkout MiMB`` inside the cloned folder. The default branch contains the current version described in these pages.
