.. _Installation:
  
=====================
Installation
=====================

To install MicroExonator follow these instructions:

Clone repository
=================
Clone the github repository

.. code-block:: bash

  git clone https://github.com/hemberg-lab/MicroExonator

Install Miniconda 3

.. code-block:: bash

   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   chmod +x Miniconda3-latest-Linux-x86_64.sh ./Miniconda3-latest-Linux-x86_64.sh
  


Set up a master virtual environment
===================================

Create a conda virtual enviroment with the necesary dependencies

.. code-block:: bash

  conda create -n snakemake_env -c bioconda -c conda-forge snakemake


