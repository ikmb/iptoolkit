# The immunopeptidomic toolkit library, IPTK # 

### Introduction and Project Aim ###
<p>IPTK is a Python library specialized in the analysis of HLA-peptidomes identified through an immunopeptidomic(IP) pipeline. 
The library provides a high level API for analyzing and visualizing the identified peptides, integrating transcritomics and protein structure information 
for a rich analysis of the identified immunopeptidomes. It also provides a toolbox for integrating and comparing different experiments and/or different runs.</p>

### Release 0.4.0 notice:
<p> 1- Adding function to compute immunopeptiomic coverage matrix </p>
<p> 2- Introducing MDS plots for comparing the similarities between runs based on immunopeptidomic coverage </p>   

### Release 0.4.6 notice:
<p> Minor corrections in the documentation and the default values for some parameters in the visualization functions</p>

### Release 0.4.7 notice:
<p> Minor corrections in the visualization module</p>

### Release 0.4.8 notice:
<p> Corrected a bug in the Peptide class to manage peptides containing parentheses in the sequence. This bug caused the len function to return the number of characters 
in the sequence instead of the number of amino acids. </p>

### Release 0.4.10 notice:
<p> Corrected a bug in the Experiment class to correctly compute the length of peptides containing parentheses. This bug caused the len function to return the number of characters in the sequence instead of the number of amino acids. </p>

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
