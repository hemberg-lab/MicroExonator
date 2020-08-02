.. differential_inclusion_analysis


===============================
Differential inclusion analysis
===============================


On this secction we descrive the a downstream module that was developed to perform alternative splicing analysis between sample groups. To quantify and assess differential inclusion of novel and annotated microexons, on this moudle we have integrated `Whippet <https://github.com/timbitz/Whippet.jl>`_, which enables a fast and accurate assesment of alterntive splicing events across user-defined sample groups.

Install
=======

To run this downstream module for the first time you need to create a environment that has `snakemake` and the version of `julia` that is compatible with `Whipet v0.11`. To creat this enviroment execute the following command inside ``MicroExonator/`` folder:

.. code-block:: bash

    conda env create -f Whippet/julia_0.6.1.yaml

Then, activate the newly created enviroment:

.. code-block:: bash

    source activate julia_0.6.1

Enter julia's interactive mode:

.. code-block:: bash

    julia

Install Whippet by excecuting the following command on the interactive session:

.. code-block:: bash

    Pkg.add("Whippet")

.. note::

    To exit julia interactive session press ``control + d``.


Configure
=========

Here there is an list of the additonal keys that need to be incorporated as a part of config.yaml:

.. code-block:: bash
    
    whippet_bin_folder : /path/to/miniconda/envs/julia_0.6.1/share/julia/site/v0.6/Whippet/bin
    Gene_anontation_GTF : /path/to/gene.annotation.gtf
    whippet_delta : /path/to/whippet_delta.yaml

* ``whippet_bin_folder`` correspodn t the path of whippet binary folder (``Whippet/bin``) that is located inside ``julia_0.6.1`` virtual enviroment folder. The specific routh to ``Whippet/bin`` may variate, so it is important that you manually identify the correct path.

* ``Gene_anontation_GTF`` corresponds the path of a gene annotation file as Gene Transfer Format (`GTF <https://en.wikipedia.org/wiki/Gene_transfer_format#:~:text=The%20Gene%20transfer%20format%20(GTF,conventions%20specific%20to%20gene%20information.>`_). Working with the same annotation data base than the one used on the previous steps is recommended. 

* ``whippet_delta`` indicate the path of a `YAML <https://en.wikipedia.org/wiki/YAML#:~:text=Open%20format%3F&text=YAML%20(a%20recursive%20acronym%20for,is%20being%20stored%20or%20transmitted.>`_ file you need to create to provide information about the desired comparisons between groups of samples.


whippet_delta YAML file
-----------------------

This file can contain the information to schedule any number of comparison between sample groups of any size. Every comparison should have the following structure inside the YAML file:

.. code-block:: bash

    comparison_ID:
      A : sample1,sample2,sample3
      B : sample4,sample5,sample6

Where ``sample1 ... sample6`` correspond to base names given to each RNA-seq samples at the corresponding input files (See :doc:`setup`) and `comparison_ID` to any given name for the sheduled comparison. As an example see the :download:`YAML file <../../Examples/Runs/Parada_et_al/whippet_delta.yaml>` we used in our publication. 

.. warning::

    Inside this YAML file sample groups must be named ``A`` and ``B``.


Optional parameters
-------------------

If you just want to skip Discovery and Quantification modules and just asses alternative splicing events annotated at the provided GTF file, then include the following like at the configuratio file:

.. code-block:: bash

    downstream_only : T

Output
======

Quantification files generated per each sample can be found at ``Whipet/Quant``. Differentially included microexon analyses that can be obtained with Whippet, are reported at ``Whippet/Delta`` folder. MicroExonator performs these analyses using both PSI values calculated internally by the pipeline and PSI values directly calculated with Whippet. These results are reported under the same format than the ``diff.gz`` descrived at the `Whippet's GitHub page <https://github.com/timbitz/Whippet.jl#output-formats>`_. However, to provide easier interpretation, we filter the Whippet splicing nodes that correspond to microexon inclusion events, these are reported as ``.microexons`` files, where ``.diff.ME.microexons`` files correspond to the output when MicroExonator PSI values are taken as input and ``.diff.microexons`` when Whippet PSI  values are taken as input.

