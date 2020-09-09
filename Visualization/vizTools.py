#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: visualizaton functions  
@version: 0.0.1
"""
# import the module 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 
import seaborn as sns 
import logomaker as lgm 
from colour import Color 
from Bio.PDB import MMCIFParser
import os 
import nglview as nv 
from nglview.color import ColormakerRegistry 
import pandas as pd 
import numpy as np 
from IPTK.Utils.Types import PlottingKeywards, MappedProteinRepresentation
from typing import List, Dict 
from scipy.stats import pearsonr
# define some helper and formater functions 
@ticker.FuncFormatter
def major_formatter(x,pos):
    """
    @brief: an override for the y-axis major ticke formater that removes the negative sign
    @note: this code was adapted from the great solution posted on stackoverflow 
    https://stackoverflow.com/questions/51086732/how-can-i-remove-the-negative-sign-from-y-tick-labels-in-matplotlib-pyplot-figur
    """
    if x<0:
        return str(-x)
    return str(x)
# define the plotting functions 
def plot_overlap_heatmap(results_df:pd.DataFrame, plotting_kwargs: PlottingKeywards)->sns.matrix.ClusterGrid:
    """
    @brief: plot the peptide/protein overlap clustermap
    @param: results_df: a pandas dataframe table that hold the overlapping results.
    """
    fig=sns.clustermap(results_df, **plotting_kwargs)
    return fig

def plot_motif(pwm_df:pd.DataFrame, plotting_kwargs:PlottingKeywards = {'fade_probabilities':True})->plt.Figure:
    """
    @brief: a generic motif plotter that forward its argument to logomaker for making plots 
    @param: pwm_df: a pandas dataframe containg the position weighted matrix 
    @param: plotting_kwargs: a dictionary of parameter to control the behevier of 
    the logo plotter. 
    @see: logomaker.Logo
    """
    res=lgm.Logo(pwm_df,**plotting_kwargs)
    return res

def plot_paired_represention(protein_one_repr , protein_two_repr, 
                            color_first: str = 'red', color_second: str = 'blue',
                            alpha: float = 0.9 ) -> plt.Figure:  
    """ 
    @brief a defination for a paired representation function
    @param: protein_one_repr: a dict object containing the legand of the first protein along with its mapped array
    @param: protein_two_repr: a dict object containing the legand of the second protein along with its mapped aray  
    @param: color_first: the color of representation for the first protein 
    @param: color_second: the color of the second protein 
    @param: alpha: the transparancy of the figure 
    @example:
    >>> p1={'cond1':np.random.randint(low=0,high=100,size=(100,))}
    >>> p2={'cond2':np.random.randint(low=0,high=100,size=(100,))}
    >>> res=plot_paired_represention(p1,p2)
    >>> res.show() # show the results figure 
    >>> res.savefig('example_one.png',dpi=600) # save the figure with the res.
    """ 
    # get the protein names 
    name_protein_one=list(protein_one_repr.keys())[0] 
    name_protein_two=list(protein_two_repr.keys())[0]
    # assert that both protein have the same length 
    if len(protein_one_repr[name_protein_one]) != len(protein_two_repr[name_protein_two]):
        raise ValueError(f"""The length of the provided proteins does not match!, they have length of: 
    {len(protein_one_repr[name_protein_one])} and {len(protein_two_repr[name_protein_two])}""")
    # create the figure and the axes 
    fig, ax = plt.subplots()   
    # plot the boundaries
    ax.plot(protein_one_repr[name_protein_one],color_first, label=name_protein_one) 
    ax.plot(-1*protein_two_repr[name_protein_two], color_second, label=name_protein_two) 
    # plot the filling 
    protein_pos=range(len(protein_two_repr[name_protein_two])) 
    ax.fill_between(protein_pos, protein_one_repr[name_protein_one], color='red', alpha=0.5) 
    ax.fill_between(protein_pos, -1*protein_two_repr[name_protein_two], color='blue', alpha=0.5) 
    # Adjust the formats of the figure 
    ax.set_ylabel('Coverage') 
    ax.set_xlabel('Amino acid Position')
    ax.yaxis.set_major_formatter(major_formatter) 
    ax.legend() 
    return fig 

def plot_overlay_representation(proteins, alpha: float = 0.5, title: str = None, 
                                legend_pos: int = 2)-> plt.Figure:
    """
    @brief: plot an overlayed representation for the SAME protein or proteins OF EQAUL LENGTH 
    in different conditions, in which the mapped represention of each protein are stacked ontop of each other to generate 
    a representation for the protein representability under different conditions. 
    @param: proteins: a nested dict object containing for each protein a child dict that contain 
    the mapping array and the color in the figure. 
    @param: title: The title of the figure, default is None. 
    @param alpha: the transparancy between proteins 
    @example: 
    >>> simulated_coverage=simulate_protein_representation(4,100,50)
    >>> overlayed_fig=plot_overlay_representation(simulated_coverage,alpha=0.2
    >>> overlayed_fig.savefig('TestFile.png',dpi=1200)
    """
    # create the figure 
    fig=plt.figure()
    # assert that the proteins are of equal length 
    names=list(proteins.keys())
    len_protein_zero= len(proteins[names[0]]['mapped_protein'])
    for name in names: 
        if len(proteins[name]['mapped_protein']) != len_protein_zero: 
            raise ValueError('The provide proteins are not of equal length')
    # plot the mapping to the protein backbone 
    for name in names: 
        # plot the boundary 
        plt.plot(proteins[name]['mapped_protein'], 
        color=proteins[name]['color'], label=name)
        # plot the filling 
        plt.fill_between(range(len_protein_zero),
        proteins[name]['mapped_protein'],alpha=alpha)
    # add the figure labels 
    plt.legend(loc=legend_pos)
    plt.xlabel('Amino acid position')
    plt.ylabel('Coverage')
    if title !=None: plt.title(title)
    return fig 

def plot_protein_coverage(mapped_protein: np.ndarray, col: str = 'r', prot_name: str = None)->plt.Figure:
    """
    @brief: plot the mapped protein array 
    @param:  mapped_protein: a numpy array with shape of 1 by protein length
    @param: col: the color of the coverage respresentation
    @param: prot_name: the default protein name
    """
    fig=plt.figure()
    plt.plot(mapped_protein, col)
    plt.fill_between(range(len(mapped_protein)), mapped_protein, color=col)
    plt.xlabel('Amino acid position')
    plt.ylabel('Coverage')
    if prot_name != None: 
        plt.title(f'Coverage representation for: {prot_name}')
    return fig

def plot_protein_presentation_3D(proteins, plotting_args = {'color':'red'}, title: str = None)-> plt.Figure: 
    """
    @brief plot a 3D surface representation for the SAME protein or proteins OF EQAUL LENGTH 
    in different conditions.  
    @param: proteins: a nested dict object containing for each protein a child dict that contain 
    the mapping array and the color in the figure.  
    @param: plotting_args: a dict that contain further parameter to the plot_surface functions
    """
    # create the figure 
    fig=plt.figure()
    # assert that the proteins are of equal length 
    names=list(proteins.keys())
    len_protein_zero= len(proteins[names[0]]['mapped_protein'])
    for name in names: 
        if len(proteins[name]['mapped_protein']) != len_protein_zero: 
            raise ValueError('The provide proteins are not of equal length')
    # get a list of the proteins names &  their mapped-array 
    presentability_score=[]
    protein_names=list(proteins.keys())
    for protein in protein_names:
        presentability_score.append(proteins[protein]['mapped_protein'].reshape(1,-1))
    # compute the presentability score 
    presentation_matrix=np.concatenate(presentability_score)
    # a 2D matrix of shape (num_of proteins, protein length)
    _num_cond=np.arange(presentation_matrix.shape[1])
    _prot_idx=np.arange(presentation_matrix.shape[0])
    X, Y=np.meshgrid(_num_cond,_prot_idx)
    # plot the results 
    ax=plt.axes(projection='3d')
    ax.plot_surface(X,Y,presentation_matrix,**plotting_args)
    # add the key element of the figure 
    ## add the y-ticks
    ax.set_yticks(list(range(len(protein_names))))
    ax.set_yticklabels(protein_names)
    ## add the y-label
    ax.set_ylabel('Conditions')
    ax.set_xlabel('Amino Acid position')
    ax.set_zlabel('Coverage')
    # set the title 
    if title != None: plt.title(title)
    return fig

def imposed_coverage_on_3D_structure(path2mmCIF: str, mapped_protein: np.ndarray, 
                                    background_color: str ='black', low: str ='red', high: str ='violet')->None: 
    """
    @brief: A function to impose the peptide coverage ontop of a protein 3D structure
    where the color of each position is marked by a color gradient that reflect the number of peptides 
    aligned to this position. 
    @Note: Use the function with Jupyter-note book as it return an NGLWidget object that your can explore 
    and navigate from the browser. 
    @param: path2mmCIF: The path to load the mmCIF file containing the protein structure
    @param: mapped_protein: a numpy array of containg the number of peptides originated from each position in the array
    @param: background_color: the color of the background, default is black 
    @param: low: the color of low covered position, default is red. 
    @param: high: the color of high covered position, default is violet.
    @param: in_notebook: wether or not the function is called from a notebook or not 
    """
    # Check the file exists 
    if not os.path.exists(path2mmCIF): 
        raise ValueError(f'The provided path to the mmCIF file: {path2mmCIF} does not exist!')
    # Load the file 
    parser=MMCIFParser()
    # Get the structure 
    structure=parser.get_structure('VIS_STRUC',path2mmCIF)
    lowest_coverage=np.min(mapped_protein)
    highest_coverage=np.max(mapped_protein)
    difference=highest_coverage-lowest_coverage
    # Generate the gradient
    low_col=Color(low)
    color_gradient=list(low_col.range_to(Color(high),difference))
    # generate the color of each position 
    colors=[]
    aa_counter=0 # coloring head is pointing to position zero 
    for pos in mapped_protein.ravel(): 
        color_gradient.append([str(color_gradient[pos]),aa_counter])
        aa_counter+=1 # set the head to point to thenext amino acids 
    ## visualize the results using NGLview 
    # set the color map
    cm = ColormakerRegistry
    cm.add_selection_scheme('VIS_STRUC', colors)
    ## view the file
    view=nv.show_biopython(structure)
    view.clear_representations()
    view.add_cartoon(color='VIS_STRUC')
    view.backbone=background_color
    view
    return view  

def plot_peptide_length_dist(pep_length: List[int],
                         plotting_kwargs: Dict[str,str]={}, 
                         x_label: str = 'Peptide Length',
                         y_label: str = 'Frequency',
                         title: str = 'Peptide Length distribution'):
    """
    @brief: visualize a histogram of the eluted peptide length 
    @param: pep_length: a list of integer containing the peptides' lengths
    @param: plotting_kwargs: a dict object containing parameters for the function
    seaborn::distplot
    @param: x_label: the label of the x-axis 
    @param: y_label: the label of the y-axis 
    @param: title: the title of the figure
    """
    fig=plt.figure()
    ax=sns.distplot(pep_length,**plotting_kwargs)  
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    return fig

def plot_num_peptides_per_parent(nums_table: pd.DataFrame,
                        num_prot: int = -1, 
                        plotting_kwargs: Dict[str,str]={}, 
                        x_label: str = 'Number of peptides',
                        y_label: str = 'Protein ID',
                        title: str = 'Number of peptides per protein'):
    """
    @brief: visualize a histogram of the eluted peptide length.  
    @param: nums_table: a pandas dataframe containing number of peptides identified from each protein. 
    @param: num_prot, the number of protein to show relative to the first element, for example, the first 10, 20 etc.
    If the default value of -1 is used then all protein will be plotted, however, this might lead to a crowded figure.
    @param: plotting_kwargs: a dict object containing parameters for the function
    seaborn::barplot, enable more finetune control of the function behavior. 
    @param: x_label: the label of the x-axis 
    @param: y_label: the label of the y-axis 
    @param: title: the title of the figure
    """ 
    # check the number of proteins to plot 
    if num_prot !=-1: 
        if num_prot > nums_table.shape[0]: 
            raise ValueError(f'The provided protein number of proteins to plot: {num_prot} is bigger than the number of proteins in the provided table: {nums_table.shape[0]}')
        nums_table=nums_table.iloc[:num_prot,]
    # plot the results 
    fig=plt.figure()
    ax=sns.barplot(x='Number_of_Peptides', y='Proteins', data=nums_table, **plotting_kwargs)
    ax.set_xlabel('Number of peptides')
    ax.set_ylabel('Protein ID')
    ax.set_title(title)
    return fig

def plot_parent_protein_expression_in_tissue(expression_table: pd.DataFrame, 
    tissue_name: str, plotting_kwargs: Dict[str,str]={'orient':'v'}, def_value: float = -1,
    ylabel: str = 'Normalized Expression', title: str = 'Parent proteins gene expression') -> plt.Figure:
    """
    @brief: plot the parent protein expression in tissue 
    @param: expression_table: The protein expression table which contains the expresion value for each parent protein
    @param: tissue_name: The name of the tissue 
    @param: def_value: The default value for proteins that could not be mapped to the expression database 
    @param: plotting_kwargs: a dict object containing parameters for the sns.violinplot function.
    @param: ylabel: the label on the y-axis. 
    @param: title: the title of the figure.
    """
    # First filter the DB for the non-mapped 
    df=expression_table.loc[expression_table.iloc[:,1]!=def_value,]
    # get the num of un-mapped
    num_un_mapped: int = expression_table.shape[0]-df.shape[0]
    # create a figure to plot to it 
    fig= plt.figure()
    ax=sns.violinplot(df.iloc[:,1], **plotting_kwargs)
    # set the axis of the axes, legend, label, etc 
    ax.set_ylabel(ylabel+' in '+tissue_name)
    ax.set_xlabel(f'Number of proteins: {expression_table.shape[0]}, Number of proteins without reference expression value: {num_un_mapped}')
    ax.set_title(title)
    return fig


def plot_num_protein_per_org(counts_table: pd.DataFrame,
                            plotting_kwargs: Dict[str,str]={}
                            )->plt.Figure:
    """
    """
    pass 

def plot_gene_expression_vs_num_peptides(exp_count_table: pd.DataFrame, tissue_name: str,
    def_value: float = -1, plotting_kwargs: Dict[str,str] = {}, 
    xlabel: str = 'Number of peptides', ylabel: str = 'Expression value', 
    title: str= 'Peptides per protein Vs. Expression Level')->plt.Figure:
    """
    @brief: plot the correlation between the gene expression and the num of peptids per protein 
    @param: exp_count_table: A table that contain the number of peptides and the expresion value for each protein in the database. 
    @param: tissue_name: The name of the tissue 
    @param: def_value: The default value for proteins that could not be mapped to the expression database 
    @param: plotting_kwargs: a dict object containing parameters for the sns.scatter function.
    @param: ylabel: the label on the y-axis. 
    @param: title: the title of the figure.
    @param: title: the title of the figure.
    """
    # First filter the DB for the non-mapped 
    df=exp_count_table.loc[exp_count_table.iloc[:,2]!=def_value,]
    # create a figure to plot to it 
    fig= plt.figure()
    ax=sns.scatterplot(df.iloc[:,1],df.iloc[:,2], **plotting_kwargs)
    ax.set_ylabel(ylabel+' in '+tissue_name)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    return fig

def plot_num_protein_per_location(protein_loc: pd.DataFrame, 
    tissue_name: str, plotting_kwargs: Dict[str,str] = {}, 
    drop_unknown: bool = False, xlabel: str = 'Number of Proteins',
    ylabel: str = 'Compartment', title: str= 'Number of proteins per sub-cellular compartment'
    )->plt.Figure: 
    """
    @brief: plot the number of proteins per each sub-cellular compartment
    @param: protein_loc: A table that contain the count of protein from each location.  
    @param: tissue_name: The name of the tissue. 
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: title: the title of the figure.
    @param: title: the title of the figure.
    """
    if drop_unknown:
        protein_loc=protein_loc.loc[protein_loc.iloc[:,0]!='UNK',]
    # create a figure 
    fig=plt.figure()
    ax=sns.barplot(y='Compartment',x='Counts',data=protein_loc, **plotting_kwargs)
    # set the labels 
    ax.set_ylabel(ylabel+' in '+tissue_name)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # return the results 
    return fig










def plot_locations_vs_num_peptides():
    """
    @brief: plot the correlation between the sub-cellular location and the number of observed peptides 
    """
    pass 

def plot_num_parent_per_peptide():
    """
    @brief: plot the number of peptids of parent protein
    """
    pass 

def plot_correlation_in_motife():
    """
    @brief: plot the correlation in motif
    """
    pass 




def plot_cleavage_around_peptide(flank_region, Peptide)->None:
    """
    @brief: return a mapping of the flanking region and the peptide core. 
    n_terminal_cleavage 1         c_terminal_cleavage 1
    n_terminal_cleavage 2         c_terminal_cleavage 2
    n_terminal_cleavage 3         c_terminal_cleavage 3
        .                 Peptide         .
        .                                 .
        .                                 .
    n_terminal_cleavage n        c_terminal_cleavage n
    """
    pass 
    