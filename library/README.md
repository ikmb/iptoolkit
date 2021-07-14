# The immunopeptidomic toolkit library, IPTK #

## Introduction and Project Aim ##

<p>IPTK is a Python library specialized in the analysis of HLA-peptidomes identified through an immunopeptidomic(IP) pipeline.
The library provides a high level API for analyzing and visualizing the identified peptides, integrating transcritomics and protein structure information for a rich analysis of the identified immunopeptidomes. It also provides a toolbox for integrating and comparing different experiments and/or different runs.</p>

## Release 0.6 notice ##

<p> Version 0.6 brings major upgrades to the library and introduce a wide array of function and classes for automating and accelerating IPTK performance </p>
<p> 1- IPTK can now parse and work with mzIdentML files using the function parse_mzIdentML_to_identification_table define in the IO module of the library </p>
<p> 2- IPTK can now process and read mzML files directly using PyOpenMS </p>
<p> 3- IPTK has an improved function executional speed thanks to the AcceleratedFunctions module in the Analysis module which provides an acceleration using Numba</p>
<p> 4- Current release also introduce, the Wrappers module which provide a simple abstraction for creating Experiment and ExperimentSet</p>
<p> 5- Introducing ReplicatedExperiments which provides a simple API for creating experiments obtained from replicates </p>
<p> 6- IPTK, current support concurrent execution, the wrapper submodules, now utilizes multiprocessing for parsing and reading multiple datasets on-parallel </p>
<p> 7- Introducing, chordDiagram for showing overlap among experiments and Proband of experiments </p>

## The road to version 1.0 ##

<p> The major plan is to, first, increase and enhance IPTK scale and execution speed by offloading computational intensive tasks to RUST. Second, increase automation by providing custom analysis recipes for performing commonly used routines. Third, provide an API for integrating other omics layers, namely metabolomics and proteomics. Finally, adding support to PTM modified HLA peptides and proteins</p>

### Planned features for 0.7.* Release ###

<p> 1. Release 0.7.1 will aim at supporting the integration of Proteomic data with the library </p>
<p> 2. Release 0.7.2 will aim at supporting the integration of Metabolomics data with the library</p>
<p> 3. Release 0.7.3 will aim at standardizing all omics API and provide a high-level abstraction for working with them</p>
<p> 4. Release 0.7.4-0.7.7 will aim at re-implement all the class in Rust and provide a python wrapper around these classes, Thus ensuring fast and concurrent execution</p>

### Planned features for 0.8.* ###

<p> 1. Release 0.8.1-0.8.4 will aim at re-implement all IPTK parsers in Rust and provide a python binder to it</p>
<p> 2. Release 0.8.5-0.8.8 will aim at re-implementing all analysis function using Rust </p>

### Planned features for 0.9.* ###

<p> Different minor releases will introduce different analysis Recipes to automate analysis tasks</p>

### Planning for version 1.0.0 ###

<p> IPTK version 1.0 is release on PyPi and on BioConda </p>

## Previous versions Release notice ##

### Release 0.5 notice: ###

<p> 1- Adding a class to query AFND database for allele frequency world-wide. </p>
<p> 2- Adding function for plotting a choropleth for allele frequencies. </p>
<p> 3- Adding classes for working directly with mzML files using pyopenMS framework </p>
<p> 4- An experimental class that act as database interface and provide method for storing and querying immunopeptidomic data</p>

### Release 0.4.11 notice: ###

<p> Adding more control to the function plot_MDS_from_ic_coverage to fine-tune its behavior, for example, by controlling the random seed.</p>

### Release 0.4.10 notice: ###

<p> Corrected a bug in the Experiment class to correctly compute the length of peptides containing parentheses. This bug caused the len function to return the number of characters in the sequence instead of the number of amino acids. </p>

### Release 0.4.8 notice: ###

<p> Corrected a bug in the Peptide class to manage peptides containing parentheses in the sequence. This bug caused the len function to return the number of characters 
in the sequence instead of the number of amino acids. </p>

### Release 0.4.7 notice: ###

<p> Minor corrections in the visualization module</p>

### Release 0.4.6 notice: ###

<p> Minor corrections in the documentation and the default values for some parameters in the visualization functions</p>

### Release 0.4.0 notice: ###

<p> 1- Adding function to compute immunopeptiomic coverage matrix </p>
<p> 2- Introducing MDS plots for comparing the similarities between runs based on immunopeptidomic coverage </p>  


### Tutorials ###

<p>The library have three notebooks that provide a step-by-step guidance to use the library and to utilize its major APIs for interacting with an IPs data.
These tutorial can be found at the Tutorial directory at the project's Github page.</p>

<p> IPTK has been documented using Sphinx, the manual of the library can be found at docs directory and online at <a href= "https://iptk.readthedocs.io/en/latest/index.html"> readthedocs </a> </p> 

### Installation ###

<p>The library can be installed using pip as follows: </p>

```
pip install iptkl --user
```

### Funding ###

The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’). 
