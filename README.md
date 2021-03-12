![IKMB_LOGO](/Media/IKMB_LOGO.png)
# The immunopeptidomic toolkit library, IPTK # 

### Introduction and Project Aim ###
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
<p> numpy, pandas, biopython, seaborn, matplotlib, plotly, mhcnames, pyteomics, h5py, logomaker, colour, lxml, nglview, sklearn, scipy, statannot. levenshtein</p>
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


### Get Started! ### 
<p>The library has four notebooks that provide a step-by-step guidance to use the library and to utilize its major APIs for interacting with an IPs data.
These tutorials can be found at the Tutorial directory</p>

<p> IPTK has been documented using Sphinx, the manual of the library can be found at the docs directory and online at <a href= "https://iptk.readthedocs.io/en/latest/index.html"> readthedocs </a> </p> 

#### Running the tutorials  #####
<p>
To run the tutorials locally, run the following steps: 
<p>1- Download the tutorials, either by cloning the repository or by downloading the tutorials along with the associated datasets only.</p>
<p>2- Start the notebook by running jupyter-notebook from the terminal. </p>
</p> 


#### Running the dashboard #### 
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

### Release 0.4 notice:
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

### Release 0.4.11 notice:
<p> Adding more control to the function plot_MDS_from_ic_coverage to fine-tune its behavior, for example, by controlling the random seed.</p>

### Release 0.5 notice: 
<p> 1- Adding a class to query AFND database for allele frequency world-wide. </p>
<p> 2- Adding function for plotting a choropleth for allele frequencies. </p>
<p> 3- Adding classes for working directly with mzML files using pyopenMS framework </p>

### Funding ###
The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’) 

![IKMB_LOGO](/Media/RTG1743.png)