.. differential_inclusion_analysis


===============================
Differential inclusion analysis
===============================


Differentially included microexon analyses that can be obtained with Whippet, are reported at ``Whippet/Delta`` folder. MicroExonator performs these analyses using both PSI values calculated internally by the pipeline and PSI values directly calculated with Whippet. These results are reported under the same format than the ``diff.gz`` descrived at the `Whippet's GitHub page <https://github.com/timbitz/Whippet.jl#output-formats>`_. However, to provide easier interpretation, we filter the Whippet splicing nodes that correspond to microexon inclusion events, these are reported as ``.microexons`` files, where ``.diff.ME.microexons`` files correspond to the output when MicroExonator PSI values are taken as input and ``.diff.microexons`` when Whippet PSI  values are taken as input. 

