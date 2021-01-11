![IKMB_LOGO](/Media/IKMB_LOGO.png)
# The immunopeptidomic toolkit library, IPTK # 

### Introduction and Project Aim ###
<p>IPTK is a Pythonic library specialized in the analysis of HLA-peptidomes identified through an Immunopeptidomics pipeline. 
The library provides a high-level API for analyzing and visualizing the identified peptides, integrating transcriptomics and protein structure information 
for a rich analysis of the identified immunopeptidomes. It also provides a toolbox for integrating and comparing different experiments and/or different runs. </p>

### Release 0.4.0 notice:
<p> 1- Adding function to compute immunopeptidomics coverage matrix </p>
<p> 2- Introducing MDS plots for measuring the similarities among different runs based on immunopeptidome coverage. </p>   

### Tutorials ### 
<p>The library have three notebooks that provide a step be step guidance to use the library and to utilize its major APIs for handling different aspects of IP identification.
These Tutorial can be found at the <a href= "https://github.com/ikmb/iptoolkit/tree/master/Tutorials"> Tutorials </a></p>

<p> IPTK have be documented using Sphnix, the manual of the library can be found at docs directory and online at <a href= "https://iptk.readthedocs.io/en/latest/index.html"> readthedocs </a> </p> 


### Installation ###
<p>The library can be installed using pip as follows: </p> 

```
pip install iptkl --user
```

### Notes and common troubleshooting ### 
<p> 1- Please make sure the pip is probably installed on the system </p>
<p> 2- For macOS users please make sure Xcode is installed. This can be done using the following command</p>

```
xcode-select --install 
```

<p> For Debian/Ubuntu users, please make sure build-essential installed. </p>

```
sudo apt install build-essential
```

### Running the Tutorials  ####
<p>
To run the tutorial locally, please run the following steps: 
<p>1- Download the tutorials, either by cloning the repository or by downloading individual tutorials along with the associated datasets.</p>
<p>2- Start the notebook by running jupyter-notebook from the terminal.</p>
</p> 

### Dependencies ###
<p> The library requires the following libraries to be installed in order to functional properly: </p>
<p> numPy, pandas, biopython, seaborn, matplotlib, plotly, mhcnames, pyteomics', h5py, logomaker, colour, lxml, nglview, sklearn, scipy, statannot. levenshtein</p>
<p> Usually this packages are installed automatically through pip. However, incase this process failed, the dependcies can be installed using the as follows:</p>

```
pip install -r requirements.txt 
```


### Funding ###
The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’) 

![IKMB_LOGO](/Media/RTG1743.png)