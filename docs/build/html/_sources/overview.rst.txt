.. overview
  
========
Overview
========

MicroExonator is a fully-integrated computational pipeline that allows for systematic de novo discovery and quantification
of microexons using raw RNA-seq data for any organism with a gene annotation. Compared to other available methods MicroExonator
is more sensitive for discovering smaller microexons and it provides higher specificity for all lengths. Moreover, MicroExonator
provides integrated downstream comparative analysis between cell types or tissues using
`Whippet <https://github.com/timbitz/Whippet.jl>`_. (`Sterne-Weiler et al. 2018 <https://doi.org/10.1016/j.molcel.2018.08.018>`_). 
As a proof of principle MicroExonator  identified X novel microexons in Y RNA-seq samples from mouse early development to provide a systematic characterization 
based on time and tissue specificity.

MicroExonator pipeline is divided in several modules:
    * Discover
    * Quantification
    * Differential Inclusion
    * Single cell analysis

