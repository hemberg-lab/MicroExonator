.. _input_files:

===========  
Setup
===========

Before runnung MicroExonator there are several files that needs to be created inside ``MicroExonator/`` root folder:

RNA-seq samples
===============

Input RNA-seq data either a ``local_samples.tsv``, ``NCBI_accession_list.txt`` or ``sample_url.tsv`` needs to be defined.
If you want to run MicroExonator over RNA-seq samples that are locally stored, they need to be defined inside ``local_samples.tsv``.
MicroExonator can also download and run samples from NCBI if the corresponding SRA accession names are defined inside of ``NCBI_accession_list.txt``,
in addition any ``fastq.gz`` that can be directly download from a URL can be included into the aalysis by defining them inside a ``sample_url.tsv``.
You can find examples of these files inside the ``Examples/`` folder.
Is posible to combine different types of input sources, but at least one of these files needs to be defined inside ``MicroExonator/`` root folder. 

Cluster configuration
=====================

If you are working on a high performace cluster, then it is very likely that you need to submit jobs to queueing systems such as lsf, qsub, SLURM, etc.
To make MicroExonator work with these queueing systems, you need to create a `cluster.json` file. 
We currently provide in the Examples folder a ``cluster.json`` file to run MicroExonator with `lsf <https://www.ibm.com/support/knowledgecenter/en/SSETD4/product_welcome_platform_lsf.html>`_.
To adapt MicroExonator to other quequing systems please see the `SnakeMake documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html?highlight=cluster.json#cluster-configuration>`_.

Config file
===========

Each MicroExonator's module has certain compulsory and optional parameters that needs to be defined inside a ``config.yaml`` file.
The necesary content of ``config.yaml`` is described on each moudle section and examples can be found at the ``Examples/`` folder.
