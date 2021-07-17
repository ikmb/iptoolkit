![IKMB_LOGO](/Media/IKMB_LOGO.png)

# The immunopeptidomic toolkit library, IPTK # 

## Introduction and Project Aim ##

<p>IPTK is a Python library specialized in the analysis of HLA-peptidomes identified through an Immunopeptidomic(IP) pipeline. 
The library provides a high level API for analyzing and visualizing the identified peptides, integrating transcritomics and protein structure information 
for a rich analysis of the identified immunopeptidomes. It also provides a toolbox for integrating and comparing different experiments and/or different runs.</p>

### Installation ###

<p>The library can be installed using pip as follows: </p> 

```
pip install iptkl --user
```

### Notes and common troubleshooting ###

<p> 1- Please make sure that pip is installed on the system. </p>
<p> 2- For macOS users please make sure Xcode is installed. This can be done using the following command</p>

```
xcode-select --install 
```

<p> For Debian/Ubuntu users, please make sure build-essential is installed. </p>

```
sudo apt install build-essential
```

<p> 3- make sure you run the dashboard with python>=3.6, using conda or pip as follow </p>

```
pip install python==3.6
```

### Dependencies ###

<p> The library requires the following libraries to be installed in order to function properly: </p>
<p> numpy, pandas, biopython, seaborn, matplotlib, plotly, mhcnames, pyteomics, h5py, logomaker, colour, lxml, nglview, sklearn, scipy, statannot</p>
<p> Usually these packages are installed automatically through pip. However, incase this process failed, the dependencies can be installed as follows:</p>

```
pip install -r requirements.txt 
```

### Visualization ###

<p> 1. Incase you are working within a Jupyter Notebooks, you can set the magic command %matplotlib notebook to work interactively with the generated plots.
However, if you are working on an IPython shell, please add the magic command %matplotlib to work. </p>

<p> 2. To save, any of the figures generated using matplotlib or seaborn, use the following command: </p>

```
fig.savefig('my_figure_name.my_extension', dpi=600)
```

<p> 3. To visualize any of the figures generated using plotly library, use either:  </p>

```
fig.show()
```

<p> which will open the figure in the browser for interactive visualization and the generated figure can then be saved from there. A second option is to save the figure directly in python, as follows </p>

```
fig.write_image('my_figure_name.my_extension')
```

<p> 4. To work interactively using Plotly based figures inside Jupyter Notebooks: </p>

<p> A. install chart studio as follows: </p>

```
pip install chart_studio 
```

<p> B. embed the generated plotly figure using the function chart_studio.plotly.iplot as shown in tutorial 2 and 4. </p>


## Get Started! ##

<p>The library has four notebooks that provide a step-by-step guidance to use the library and to utilize its major APIs for interacting with an IPs data.
These tutorials can be found at the Tutorial directory</p>

<p> IPTK has been documented using Sphinx, the manual of the library can be found at the docs directory and online at <a href= "https://iptk.readthedocs.io/en/latest/index.html"> readthedocs </a> </p> 

### Contact ###

<p> Please feel free to write an email to the developer at h.elabd@ikmb.uni-kiel.de or to open an issue here incase of a bug or a required feature. </p>


## Running the tutorials ###

<p>
To run the tutorials locally, run the following steps: 
    <p>1- Download the tutorials, either by cloning the repository or by downloading the tutorials along with the associated datasets only.</p>
    <p>2- Start the notebook by running jupyter-notebook from the terminal. </p>
</p>

## Running the dashboard ##

<p> To start the dashboard: </p>

<p> 1. First, install IPTK, in case you have not already, as follows</p>

```
pip install iptkl --user
```

<p> 2. Install Dash and other dependencies using, as follows  </p>

```
pip install dash dash_bootstrap_components  dash-uploader 
```

<p> 3. Making the app executable, as follows </p>

```
chmod +x Apps/ExperimentUI.py 
```

<p> 4. Launch the App, as follows </p>

```
./ExperimentUI.py  
```

<p> 5. Open the app in the browser by typing the IP: http://127.0.0.1:8050/ </p>

<p> 6. Getting starting !!!</p>
<p> A simple test case can be found at the test_data directory, </p>
<p> 6.a For the identification file drag the file: 0810202_0.5_all_ids_merged_psm_perc_filtered.idXML</p>
<p> 6.b select idXML from the Format drop-down menu</p>
<p> 6.c For the Fasta Database, drag the file: human_proteome.fasta </p>
<p> 6.d Click Create Experiment, wait a few second and enjoy analyzing the data </p>

## FAQs ##

### How does IPTK map HLA types to the identified peptides? ###

<p> The class Experiment is utilized to link or connect different aspects of an experiment together, for example, the transcriptomic layer with the list of identified peptides.
In case of HLA information, it is assumed that all peptides identified in an experiments are coming from the same pool of HLA molecules, e.g. HLA-DRB1*15:01 and HLA-DRB1*13:01, incase HLA-DR
specific antibody has been used for the pulldown of HLA proteins or up to 6 HLA-I alleles, incase HLA pan specific antibodies have been used. The Experiment class is initialized with
an HLASet instance that stores information about HLA types, i.e. HLA types are linked with the list of identified peptides using the experiment class. Finally, the function compute_protein_coverage defined in the AnalysisFunction module will be used to compare protein coverage among different experiments.
 </p>

### What the difference between source code inside library directory and lib_exp_acc? ###

<p> The code inside library uses, contains stable IPTK code, meanwhile, the code inside lib_exp_acc contain experimental code, which one one hand, might contain features not included inside the stable version of the library, on the other hand, it is an experimental code so it might change tremendously between different pushes. For users of IPTK we highly recommend looking up the source code under the library directory, as it contain the source code inside the PyPi and Conda packages. On the other hand, IPTK, developers are highly recommended to work with the code inside lib_exp_acc?</p>

### Where can I get IPTK? ###

<p> IPTK can be downloaded from PyPi and Conda as described below: </p>

1. Pip based installation 

```
pip install iptkl --user
```
2. Conda based installation

```
pip install iptkl --user
```

## Release 0.6 notice ##

<p> Version 0.6 brings major upgrades to the library and introduce a wide array of function and classes for automating and accelerating IPTK performance </p>
<p> 1- IPTK can now parse and work with mzIdentML files using the function parse_mzIdentML_to_identification_table define in the IO module of the library </p>
<p> 2- IPTK can now process and read mzML files directly using PyOpenMS </p>
<p> 3- IPTK has an improved function executional speed thanks to the AcceleratedFunctions module in the Analysis module which provides an acceleration using Numba</p>
<p> 4- Current release also introduce, the Wrappers module which provide a simple abstraction for creating Experiment and ExperimentSet</p>
<p> 5- Introducing ReplicatedExperiments which provides a simple API for creating experiments obtained from replicates </p>
<p> 6- IPTK, current support concurrent execution, the wrapper submodules, now utilizes multiprocessing for parsing and reading multiple datasets on-parallel </p>
<p> 7- Introducing, chordDiagram for showing overlap among experiments and Proband of experiments </p>

### Release 0.6.6 notice ###

<p> 1- Introducing GOEngine class which provides an easy-to-use wrapper around goatools for performing GOEA on the identified proteins.</p>

<p> 2- current release supports Jaccard index as a metric of similarity among experiments</p>

<p> 3- Introducing support for visualizing GOEA results </p>

<p> 4- correction of minor bugs and documentation typos in previous releases </p>


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

### Funding ###

The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’) 

![IKMB_LOGO](/Media/RTG1743.png)