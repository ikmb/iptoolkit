![IKMB_LOGO](/Media/IKMB_LOGO.png)
# The immunopeptidomic toolkit library, IPTK # 

### Introduction and Project Aim ###
<p>IPTK is a Python library specialized in the analysis of HLA-peptidomes identified through an Immunopeptidomic(IP) pipeline. 
The library provides a high level API for analyzing and visualizing the identified peptides, integrating transcritomics and protein structure information 
for a rich analysis of the identified immunopeptidomes. It also provides a toolbox for integrating and comparing different experiments and/or different runs.</p>

### Release 0.4 notice:
<p> 1- Adding function to compute immunopeptiomic coverage matrix </p>
<p> 2- Introducing MDS plots for comparing the similarities between runs based on immunopeptidomic coverage </p>   

### Tutorials ### 
<p>The library have three notebooks that provide a step be step guidance to use the library and to utilize its major APIs for interacting with an IPs data.
These Tutorial can be found at the Tutorial directory</p>

<p> IPTK have be documented using Sphinx, the manual of the library can be found at docs directory and online at <a href= "https://iptk.readthedocs.io/en/latest/index.html"> readthedocs </a> </p> 


### Installation ###
<p>The library can be installed using pip as follows: </p> 

```
pip install iptkl --user
```

### Notes and common troubleshooting ### 
<p> 1- Please make sure that pip is probably installed on the system. </p>
<p> 2- For macOS users please make sure Xcode is installed. This can be done using the following command</p>

```
xcode-select --install 
```

<p> For Debian/Ubuntu users, please make sure build-essential installed. </p>

```
sudo apt install build-essential
```

### Running the tutorials  ####
<p>
To run the tutorials locally, run the following steps: 
<p>1- Download the tutorials, either by cloning the repository or by downloading the tutorials along with the associated datasets only.</p>
<p>2- Start the notebook by running jupyter-notebook from the terminal. </p>
</p> 

### Dependencies ###
<p> The library requires the following libraries to be installed in order to functional properly: </p>
<p> numpy, pandas, biopython, seaborn, matplotlib, plotly, mhcnames, pyteomics', h5py, logomaker, colour, lxml, nglview, sklearn, scipy, statannot. levenshtein</p>
<p> Usually these packages are installed automatically through pip. However, incase this process failed, the dependencies can be installed as follows:</p>

```
pip install -r requirements.txt 
```

### Visualization ### 
<p> Incase you are working within a Jupyter-notebooks, you can set the magic command %matplotlib notebook to work interactively with the generated plots.
However, if you are working on an IPython shell, please add the magic command %matplotlib to work. </p>

### Starting the dashboard ### 
<p> To start the dashboard: </p>

<p> 1. First, install IPTK, incase you have not already, as follows</p>

```
pip install iptkl --user
```

<p> 2. Install Dash and other dependencies using, as follows  </p>

```
pip install dash dash_bootstrap_components  dash-uploader 
```

<p> 3. Making the app executable, as follows </p>

```
chmod +x ExperimentUI.py 
```

<p> 4. Launch the App, as follows </p>

```
./ExperimentUI.py  
```

<p> 5. Open the app in the browser by typing the IP: http://127.0.0.1:8050/ </p>


### Funding ###
The project was funded by the German Research Foundation (DFG) (Research Training Group 1743, ‘Genes, Environment and Inflammation’) 

![IKMB_LOGO](/Media/RTG1743.png)