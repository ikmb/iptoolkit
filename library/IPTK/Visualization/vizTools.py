#!/usr/bin/env Python 
"""The module contain visualization functions the can be used to plot the results obtained from the 
methods of the classes defined in the Class module or from the analysis functions defined in the Analysis Module.   
"""
# import the module 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 
import os 
import math
import seaborn as sns 
import logomaker as lgm 
import nglview as nv 
import pandas as pd 
import numpy as np 
import pyopenms as poms 
import holoviews as hv 
from scipy.stats import pearsonr,ttest_ind
import plotly.express as px 
from plotly.graph_objects import Figure
import plotly.graph_objects as go 
from plotly import tools 
from colour import Color 
from Bio.PDB import MMCIFParser
from nglview.color import ColormakerRegistry 
from IPTK.Utils.Types import PlottingKeywards, MappedProteinRepresentation
from IPTK.Classes.Annotator import Annotator
from IPTK.Classes.Features import Features
from IPTK.Analysis.AnalysisFunction import get_PTMs_modifications_positions
from IPTK.Analysis.AnalysisFunction import get_PTMs_glycosylation_positions
from IPTK.Analysis.AnalysisFunction import get_PTMs_disuldfide_bonds
from IPTK.Analysis.AnalysisFunction import get_sequence_variants_positions
from IPTK.Analysis.AnalysisFunction import get_splice_variants_positions
from IPTK.Classes.MzMLExperiment import MzMLExperiment
from typing import Counter, List, Tuple, Dict, Union 
from sklearn import manifold
from bokeh.io import export_svg 
# define some helper and formater functions 
@ticker.FuncFormatter
def major_formatter(x,pos):
    """an override for the y-axis major ticker formater that removes the negative sign\
    @note: this code was adapted from the great solution posted on stackoverflow \
    https://stackoverflow.com/questions/51086732/how-can-i-remove-the-negative-sign-from-y-tick-labels-in-matplotlib-pyplot-figur
    """
    if x<0:
        return str(-x)
    return str(x)
# define the plotting functions 
def plot_overlap_heatmap(results_df:pd.DataFrame, plotting_kwargs: PlottingKeywards={})->sns.matrix.ClusterGrid:
    """ Plot a user provided dataframe as a cluster heatmap using seaborn library

    :param results_df: A pandas dataframe table that hold the overlapping number
    :type results_df: pd.DataFrame
    :param plotting_kwargs: forward parameter to the clustermap function, defaults to {}
    :type plotting_kwargs: PlottingKeywards, optional
    :return: the generated clustermap 
    :rtype: sns.matrix.ClusterGrid
    """
    fig=sns.clustermap(results_df, **plotting_kwargs)
    return fig

# define the plotly version of the library 
def plotly_overlap_heatmap(results_df: pd.DataFrame)->Figure:
    """Plot a user provided dataframe as a heatmap using plotly library

    :param results_df: A pandas dataframe table that hold the overlapping number.
    :type results_df: pd.DataFrame
    :return: a plotly Figure containing the heatmap 
    :rtype: Figure
    """
    return px.imshow(results_df) 

def plot_motif(pwm_df:pd.DataFrame, plotting_kwargs:PlottingKeywards = {'fade_probabilities':True})->plt.Figure:
    """A generic motif plotter that forward its argument to logomaker for making plots

    :param pwm_df: A pandas dataframe containing the position weighted matrix 
    :type pwm_df: pd.DataFrame
    :param plotting_kwargs: A dictionary of parameter to control the behavior of the logo plotter, defaults to {'fade_probabilities':True}
    :type plotting_kwargs: PlottingKeywards, optional
    :return: a matplotlib figure instance containing the ploted motif 
    :rtype: plt.Figure
    """
    res=lgm.Logo(pwm_df,**plotting_kwargs)
    return res

def plot_paired_represention(protein_one_repr: Dict[str, np.ndarray],
                             protein_two_repr: Dict[str, np.ndarray], 
                            color_first: str = 'red', color_second: str = 'blue',
                            alpha: float = 0.9, title=" Parired protein representation" ) -> plt.Figure:  
    """ visualize the difference in representation for the same protein between two experiments using matplotlib library.

    :param protein_one_repr: a dict object containing the legend of the first condition along with its mapped array
    :type protein_one_repr: Dict[str, np.ndarray]
    :param protein_two_repr: a dict object containing the legend of the second condition along with its mapped array
    :type protein_two_repr: Dict[str, np.ndarray]
    :param: color_first: the color of representation for the first condition 
    :type color_first: str
    :param: color_second: the color of the second condition 
    :type color_second: str
    :param alpha: the transparency of the figure. 
    :type alpha: float. 
    :param: title: the title of the figure. 
    :type title: str
    :return: A matplotlib Figure containing the representation 
    :rtype: plt.Figure
    """ 
    # get the protein names 
    name_protein_one=list(protein_one_repr.keys())[0] 
    name_protein_two=list(protein_two_repr.keys())[0]
    # adjust the correct shape of the mapped tensor 
    ## shaped protein one 
    if len(protein_one_repr[name_protein_one].shape)==2:
        protein_one_repr[name_protein_one]=protein_one_repr[name_protein_one].reshape(-1)
    ## shape protein two 
    if len(protein_two_repr[name_protein_two].shape)==2:
        protein_two_repr[name_protein_two]=protein_two_repr[name_protein_two].reshape(-1)
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
    ax.set_title(title)
    return fig 

def plotly_paired_representation(protein_one_repr: Dict[str, np.ndarray], 
                                 protein_two_repr: Dict[str, np.ndarray], 
                                 title: str =" Parired protein representation") -> Figure:
    """Compare the peptide coverage for the same protein under different conditions using the same protein using plotly library.
   
    :param protein_one_repr: a dict object containing the legend of the first condition along with its mapped array
    :type protein_one_repr: Dict[str, np.ndarray]
    :param protein_two_repr: a dict object containing the legend of the second condition along with its mapped array
    :type protein_two_repr: Dict[str, np.ndarray]
    :return: A plotly Figure containing the representation 
    :rtype: Figure
    """
    # get the protein names 
    name_protein_one=list(protein_one_repr.keys())[0] 
    name_protein_two=list(protein_two_repr.keys())[0]
    # adjust the correct shape of the mapped tensor 
    ## shaped protein one 
    if len(protein_one_repr[name_protein_one].shape)==2:
        protein_one_repr[name_protein_one]=protein_one_repr[name_protein_one].reshape(-1)
    ## shape protein two 
    if len(protein_two_repr[name_protein_two].shape)==2:
        protein_two_repr[name_protein_two]=protein_two_repr[name_protein_two].reshape(-1)
    # assert that both protein have the same length 
    if len(protein_one_repr[name_protein_one]) != len(protein_two_repr[name_protein_two]):
        raise ValueError(f"""The length of the provided proteins does not match!, they have length of: 
    {len(protein_one_repr[name_protein_one])} and {len(protein_two_repr[name_protein_two])}""")
    # create the figure 
    fig= tools.make_subplots(2,1)
    # plot the coverage of the first protein 
    fig.append_trace(go.Scatter(x=np.arange(protein_one_repr[name_protein_one].shape[0]),
    y=protein_one_repr[name_protein_one], fill='tozeroy', name=name_protein_one),1,1) 
    # plot the coverage of the second protein 
    fig.append_trace(go.Scatter(x=np.arange(protein_two_repr[name_protein_two].shape[0]),
    y=protein_two_repr[name_protein_two], fill='tozeroy', name=name_protein_two),2,1) 
    # update the title of the figure
    fig.update_layout(title=title,plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
    # return the figure 
    return fig

def plot_overlay_representation(proteins: Dict[str,Dict[str,np.ndarray]], alpha: float = 0.5, title: str = None, 
                                legend_pos: int = 2)-> plt.Figure:
    """plot an overlayed representation for the SAME protein or proteins OF EQUAL LENGTH in different conditions \
       in which the mapped represention of each protein are stacked on top of each other to generate \
       a representation for the protein representability under different conditions. 

    :param proteins: a nested dict object containing for each protein a child dict that contain the mapping array and the color in the figure. 
    :type proteins: Dict[str,np.ndarray]]
    :param alpha: the transparency between proteins, defaults to 0.5
    :type alpha: float, optional
    :param title: The title of the figure, defaults to None
    :type title: str, optional
    :param legend_pos: the position of the legend, defaults to 2
    :type legend_pos: int, optional
    :raises ValueError: if the provided protein have different lengths 
    :return: a matplotlib figure containing the mapped representation 
    :rtype: plt.Figure
    """
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

def plotly_multi_traced_coverage_representation(proteins: Dict[str,Dict[str,np.ndarray]], 
                title: str= "Protein Coverage Across  ") -> Figure : 
    """plot a multi-traced representation for the same protein across different conditions using plotly library 

    :param proteins: A dict object containing for each protein the corresponding mapped array. 
    :type proteins: Dict[str,Dict[str,np.ndarray]]
    :param title: the title of the figure, defaults to "Protein Coverage Across  "
    :type title: str, optional
    :return: a multitraced traced figure showing the coverage of proteins across different conditions
    :rtype: Figure
    """
     # un roll all the proteins to make sure they are of the correct shape 
    for prot in proteins.keys():
        proteins[prot]= proteins[prot].reshape(-1)
    # get the number of traces 
    num_traces: int = len(proteins)
    # make a figure to hold the number of proteins 
    fig = tools.make_subplots(num_traces,1)
    # adding trances to the figure 
    trace_counter: int = 1
    # loop over all the proteins and plot the results 
    for prot in proteins.keys():
        fig.add_trace(go.Scatter(x=np.arange(proteins[prot].shape[0]),
        y=proteins[prot], fill='tozeroy', name=prot,),
        trace_counter,1)
         # increase the counter 
        trace_counter+=1
    # update the title and the layout of the figure
    fig.update_layout(title=title+" "+str(trace_counter-1)+" Conditions", # correcting an indexing bug 
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)')
    # update the axis 
    fig.update_layout(
        {
            'xaxis'+str(num_traces)+'_title':'Amino Acid Position' ,
            'yaxis'+str(math.ceil(num_traces/2))+'_title':'Coverage'
        })
    # return the figure 
    return fig        

def plot_protein_coverage(mapped_protein: np.ndarray, col: str = 'r', prot_name: str = None)->plt.Figure:
    """plot the peptide coverage for a given protein. 

    :param mapped_protein: a NumPy array with shape of (1, protein length) or shape (protein-length,)
    :type mapped_protein: np.ndarray
    :param col: the color of the coverage respresentation, defaults to 'r'
    :type col: str, optional
    :param prot_name: The default protein name, defaults to None
    :type prot_name: str, optional
    :rtype: plt.Figure
    """
    if len(mapped_protein.shape)==2:
        mapped_protein=mapped_protein.reshape(-1) 
    # create a figure 
    fig=plt.figure()
    # plot the boundary
    plt.plot(mapped_protein, col)
    # plot the coverage 
    plt.fill_between(range(len(mapped_protein)), mapped_protein, color=col)
    plt.xlabel('Amino acid position')
    plt.ylabel('Coverage')
    if prot_name != None: 
        plt.title(f'Coverage representation for: {prot_name}')
    return fig

def plotly_protein_coverage(mapped_protein: np.ndarray, prot_name: str = None)->Figure:
    """Plot the peptide coverage for a given protein

    :param mapped_protein: A NumPy array with shape of 1 by protein length or shape protein-length
    :type mapped_protein: np.ndarray
    :param prot_name: The default protein name, defaults to None
    :type prot_name: str, optional
    :rtype: Figure
    """
    if len(mapped_protein.shape)==2:
        mapped_protein=mapped_protein.reshape(-1) 
    # translate the array into a dataframe 
    df_res= pd.DataFrame(mapped_protein)
    df_res.columns=['Coverage']
    # plot the figure 
    print(df_res)
    fig=px.area(df_res,y="Coverage",
            labels={'index':'Amino Acid Position'},
            title=f'Peptide Coverage of : {prot_name}')
    # update the background colors 
    fig.update_layout(
        {
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )
    # return the figure 
    return fig


def plot_protein_presentation_3D(proteins: Dict[str,Dict[str,np.ndarray]],
            plotting_args = {'color':'red'}, title: str = None)-> plt.Figure: 
    """plot a 3D surface representation for the SAME protein or proteins OF EQAUL LENGTH \
    in different conditions.  
   
    :param proteins: a nested dict object containing for each protein a child dict that contain \
    the mapping array and the color in the figure.  
    :type proteins: [type]
    :param plotting_args: a dict that contain further parameter to the plot_surface functions, defaults to {'color':'red'}
    :type plotting_args: dict, optional
    :param title: the title of the figure, defaults to None
    :type title: str, optional
    :raises ValueError: if the provided proteins are of different length 
    :rtype: plt.Figure
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
                                    background_color: str ='black', low: str ='red', high: str ='blue')->None: 
    """A function to impose the peptide coverage on top of a protein 3D structure \
    where the color of each position is marked by a color gradient that reflect the number of peptides \ 
    aligned to this position. Note: Use the function with Jupyter-note book as it return an NGLWidget object that your can explore \
    and navigate from the browser. 

    :param path2mmCIF: The path to load the mmCIF file containing the protein structure
    :type path2mmCIF: str
    :param mapped_protein: a Numpy array containing the number of peptides originated from each position in the array
    :type mapped_protein: np.ndarray
    :param background_color: the color of the background, defaults to 'black'.
    :type background_color: str, optional
    :param low: the color of low covered position, defaults to 'red'.
    :type low: str, optional
    :param high: the color of high covered position, defaults to 'blue'
    :type high: str, optional
    :raises ValueError: incase the provided path to the structure file does not exists 
    :raises IOError: if the structure file can not be loaded or if more than one file are located in the provided path 
    """
    # Check the file exists 
    if not os.path.exists(path2mmCIF): 
        raise ValueError(f'The provided path to the mmCIF file: {path2mmCIF} does not exist!')
    # check that it is a directory 
    if os.path.isdir(path2mmCIF): 
        hits_in_path: List[str] = [elem for elem in os.listdir(path2mmCIF) if 'cif' in elem]
        if len(hits_in_path) ==0: 
            raise IOError(f'No cif file can be located in the provided path: {path2mmCIF}!') 
        if len(hits_in_path)>1:
            raise IOError(f'{len(hits_in_path)} files have been found in the provided path. The function expect only one function')
        path2mmCIF= os.path.join(path2mmCIF,hits_in_path[0])
    # Load the file 
    parser=MMCIFParser()
    # check that the file is a one-D array 
    if len(mapped_protein.shape)==2: 
        mapped_protein=mapped_protein.reshape(-1)
    # Get the structure 
    structure=parser.get_structure('VIS_STRUC',path2mmCIF)
    lowest_coverage=np.min(mapped_protein)
    highest_coverage=np.max(mapped_protein)
    difference=int(highest_coverage-lowest_coverage)+1
    # Generate the gradient
    low_col=Color(low)
    color_gradient=list(low_col.range_to(Color(high),difference))
    # generate the color of each position 
    colors=[]
    #aa_counter=0 # coloring head is pointing to position zero 
    for idx in range(mapped_protein.shape[0]): 
        amino_acid_col: str = str(color_gradient[int(mapped_protein[idx])])
        colors.append([amino_acid_col,str(idx)])
    ## visualize the results using NGLview 
    # set the color map
    cm = ColormakerRegistry
    cm.add_selection_scheme('VIS_STRUC', colors)
    ## view the file
    view=nv.show_biopython(structure)
    view.clear_representations()
    view.add_cartoon(color='VIS_STRUC')
    view.background=background_color
    view
    return view  

def plot_peptide_length_dist(pep_length: List[int],
                         plotting_kwargs: Dict[str,str]={}, 
                         x_label: str = 'Peptide Length',
                         y_label: str = 'Frequency',
                         title: str = 'Peptide Length distribution'):
    """Visualize a histogram of the eluted peptide length using seaborn library. 
    
    :param pep_length: a list of integer containing the peptides' lengths
    :type pep_length: List[int]
    :param plotting_kwargs: a dict object containing parameters for the function seaborn.histplot, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param x_label: the label of the x-axis , defaults to 'Peptide Length'
    :type x_label: str, optional
    :param y_label: the label of the y-axis , defaults to 'Frequency'
    :type y_label: str, optional
    :param title: the title of the figure, defaults to 'Peptide Length distribution'
    :type title: str, optional
    """
    fig=plt.figure()
    ax=sns.histplot(pep_length,**plotting_kwargs)  
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    return fig

def plotly_peptide_length_dist(pep_length: List[int],
                         x_label: str = 'Peptide Length',
                         y_label: str = 'Counts',
                         title: str = 'Peptide Length distribution'):
    """visualize a histogram of the eluted peptide length using plotly library

    :param pep_length: a list of integer containing the peptides' lengths
    :type pep_length: List[int]
    :param x_label: the label of the x-axis , defaults to 'Peptide Length'
    :type x_label: str, optional
    :param y_label: the label of the y-axis, defaults to 'Counts'
    :type y_label: str, optional
    :param title: the figure title, defaults to 'Peptide Length distribution'
    :type title: str, optional
    """
    # define the figure 
    fig=px.histogram(pd.DataFrame({'Peptide Length':pep_length}),
                    marginal="violin",
                    labels={
                        'value':x_label,
                        'counts':y_label
                        },
                        title=title)
    # update the layout 
    fig.layout.update(showlegend=False)
    # update the background colors 
    fig.update_layout(
        {
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )
    return fig

def plot_num_peptides_per_parent(nums_table: pd.DataFrame,
                        num_prot: int = -1, 
                        plotting_kwargs: Dict[str,str]={},
                        x_label: str = 'Number of peptides',
                        y_label: str = 'Protein ID',
                        hide_y_label: bool = False, 
                        title: str = 'Number of peptides per protein'):
    """Plot the number of peptides belonging to each protein using seaborn library.  

    :param nums_table: a pandas dataframe containing number of peptides identified from each protein. 
    :type nums_table: pd.DataFrame
    :param num_prot: the number of protein to show relative to the first element, for example, the first 10, 20 etc. \
    If the default value of -1 is used then all protein will be plotted, however, this might lead to a crowded figure, \
    defaults to -1. 
    :type num_prot: int, optional
    :param plotting_kwargs: a dict object containing parameters for the function \
    seaborn::histplot, defaults to {}
    :param x_label: the label of the x-axis, defaults to 'Number of peptides'
    :type x_label: str, optional
    :param y_label: the label of the y-axis, defaults to 'Protein ID'
    :type y_label: str, optional
    :param hide_y_label: a boolean flag to hide the text of the y-label, defaults to False. 
    :type hide_y_label: bool, optional
    :param title: the title of the figure, defaults to 'Number of peptides per protein'. 
    :type title: str, optional
    :raises ValueError: if the num_prot is bigger than the number of elements in the provided table. 
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
    if hide_y_label:
        ax.get_yaxis().set_visible(False) 
        ax.set_ylabel('')
    else: 
        ax.set_ylabel('Protein ID')
    ax.set_title(title)
    return fig

def plotly_num_peptides_per_parent(nums_table: pd.DataFrame,
                        num_prot: int = -1, 
                        x_label: str = 'Number of peptides',
                        y_label: str = 'Protein ID',
                        title: str = 'Number of peptides per protein'):
    """Plot the number of peptides belonging to each protein using plotly library. 
    
    :param nums_table: a pandas dataframe containing number of peptides identified from each protein. 
    :type nums_table: pd.DataFrame
    :param num_prot: the number of protein to show relative to the first element, for example, the first 10, 20 etc. \
    If the default value of -1 is used then all protein will be plotted, however, this might lead to a crowded figure., defaults to -1\
    :type num_prot: int, optional
    :param x_label: the label of the x-axis , defaults to 'Number of peptides'
    :type x_label: str, optional
    :param y_label: the label of the y-axis , defaults to 'Protein ID'
    :type y_label: str, optional
    :param title: title, defaults to 'Number of peptides per protein'
    :type title: str, optional
    :raises ValueError: if the num_prot is bigger than the number of elements in the provided table 
    """
    # check the number of proteins to plot 
    if num_prot !=-1: 
        if num_prot > nums_table.shape[0]: 
            raise ValueError(f'The provided protein number of proteins to plot: {num_prot} is bigger than the number of proteins in the provided table: {nums_table.shape[0]}')
        nums_table=nums_table.iloc[:num_prot,]
    # plot the results 
    fig=px.bar(nums_table, 
            x="Number_of_Peptides", y="Proteins",color='Number_of_Peptides',
            color_continuous_scale='blackbody', 
            labels={
             'Number_of_Peptides':x_label,
             'Proteins':y_label   
            }
        )
    # return the background colors 
    fig.update_layout(
        {
            'title':title,
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )
    return fig

def plot_parent_protein_expression_in_tissue(expression_table: pd.DataFrame, 
    ref_expression: pd.DataFrame , tissue_name: str, sampling_num: int = 10,
    plotting_kwargs: Dict[str,str]={'orient':'v'}, def_value: float = -1,
    ylabel: str = 'Normalized Expression') -> plt.Figure:
    """Plot the parent protein expression in tissue relative a sampled collection of non-presented proteins using seaborn library.

    :param expression_table: The protein expression table which contains the expresion value for each parent protein
    :type expression_table: pd.DataFrame
    :param ref_expression: The reference expression of the tissue under investigation. 
    :type ref_expression: pd.DataFrame
    :param tissue_name: The name of the tissue . 
    :type tissue_name: str
    :param sampling_num: The number of times to sample from the non-prsenter, defaults to 10
    :type sampling_num: int, optional
    :param plotting_kwargs: a dict object containing parameters for the sns.violinplot function, defaults to {'orient':'v'}
    :type plotting_kwargs: Dict[str,str], optional
    :param def_value: The default value for proteins that could not be mapped to the expression database, defaults to -1
    :type def_value: float, optional
    :param ylabel: the label on the y-axis, defaults to 'Normalized Expression'
    :type ylabel: str, optional
    :raises ValueError: if the reference gene expression table is smaller than the number of parents 
    """
    # First filter the DB for the non-mapped 
    df=expression_table.loc[expression_table.iloc[:,1]!=def_value,]
    # get the num of un-mapped
    num_un_mapped: int = expression_table.shape[0]-df.shape[0]
    # assert that the number of genes in the database is bigger than the present number of parent proteins
    if ref_expression.shape[0] <= expression_table.shape[0]:
        raise ValueError('The provided reference gene expression table is smaller than the number of parents!')
    # extract the genes that have not-been presented 
    np_df=ref_expression.loc[~ref_expression.iloc[:,0].isin(df.iloc[:,0])].reset_index(drop=True)
    # sample the expression value 
    exp_value: np.ndarray = np.zeros((sampling_num, df.shape[0]))
    # fil the array 
    for idx in range(sampling_num):
        selected_genes: np.ndarray = list(np.random.randint(low=0,high=np_df.shape[0],
                size=(df.shape[0],)))
        # get the expression value
        exp_value[idx,:]=np_df.loc[selected_genes].iloc[:,-1].tolist()
    # compute the average representation of proteins 
    average_np_gen: np.ndarray=np.mean(exp_value,axis=0)
    # construct a dataframe that contain the presented and non-presented protein expression 
    gene_exp_p_np: pd.DataFrame = pd.DataFrame({
        'Presented':df.iloc[:,1].tolist(),
        'Not-Presented':average_np_gen.reshape(-1)
    })
    # compute the t-test between the two population
    ttest_res=ttest_ind(a= gene_exp_p_np['Presented'].tolist(),
                        b= gene_exp_p_np['Not-Presented'].tolist())
    # create a figure to plot to it 
    fig= plt.figure()
    # format the label
    title=f't-test score: {ttest_res.statistic:.3f}, P-val: {ttest_res.pvalue:.3e}'
    # plot the results 
    ax=sns.violinplot(data=gene_exp_p_np, **plotting_kwargs)
    # set the axis of the axes and the title 
    ax.set_ylabel(ylabel+' in '+tissue_name)
    ax.set_xlabel(f'Number of proteins: {expression_table.shape[0]}, Number of proteins without reference expression value: {num_un_mapped}')
    ax.set_title(title)
    return fig

def plotly_parent_protein_expression_in_tissue(expression_table: pd.DataFrame, 
            ref_expression: pd.DataFrame , tissue_name: str, sampling_num: int = 10,
            def_value: float = -1,
            ylabel: str = 'Normalized Expression') -> Figure:
    """plot the parent protein expression in tissue relative a sampled collection of non-presented proteins using plotly library.
    
    :param expression_table: The protein expression table which contains the expresion value for each parent proteins. 
    :type expression_table: pd.DataFrame
    :param ref_expression: The reference expression of the tissue under investigation.
    :type ref_expression: pd.DataFrame
    :param tissue_name: The name of the tissue. 
    :type tissue_name: str
    :param sampling_num: the number of times to sample from the non-presenters, defaults to 10
    :type sampling_num: int, optional
    :param def_value: The default value for proteins that could not be mapped to the expression database, defaults to -1
    :type def_value: float, optional
    :param ylabel: The label on the y-axis, , defaults to 'Normalized Expression'
    :type ylabel: str
    :raises ValueError: If the reference gene expression table is smaller than the number of parents 
    """
    # First filter the DB for the non-mapped 
    df=expression_table.loc[expression_table.iloc[:,1]!=def_value,]
    # get the num of un-mapped
    num_un_mapped: int = expression_table.shape[0]-df.shape[0]
    # assert that the number of genes in the database is bigger than the present number of parent proteins
    if ref_expression.shape[0] <= expression_table.shape[0]:
        raise ValueError('The provided reference gene expression table is smaller than the number of parents!')
    # extract the genes that have not-been presented 
    np_df=ref_expression.loc[~ref_expression.iloc[:,0].isin(df.iloc[:,0])].reset_index(drop=True)
    # sample the expression value 
    exp_value: np.ndarray = np.zeros((sampling_num, df.shape[0]))
    # fil the array 
    for idx in range(sampling_num):
        selected_genes: np.ndarray = list(np.random.randint(low=0,high=np_df.shape[0],
                size=(df.shape[0],)))
        # get the expression value
        exp_value[idx,:]=np_df.loc[selected_genes].iloc[:,-1].tolist()
    # compute the average representation of proteins 
    average_np_gen: np.ndarray=np.mean(exp_value,axis=0)
    # construct a dataframe that contain the presented and non-presented protein expression 
    gene_exp_p_np: pd.DataFrame = pd.DataFrame({
        'Presented':df.iloc[:,1].tolist(),
        'Not-Presented':average_np_gen.reshape(-1)
    })
    # compute the t-test between the two population
    ttest_res=ttest_ind(a= gene_exp_p_np['Presented'].tolist(),
                        b= gene_exp_p_np['Not-Presented'].tolist())
    # create the figure and plot the results 
    fig=px.violin(gene_exp_p_np, box=True) 
    # update the layout
    # return the background colors 
    fig.update_layout(
        {
            'title': f't-test score: {ttest_res.statistic}, P-val: {ttest_res.pvalue }',
            'xaxis_title':f'Number of proteins: {expression_table.shape[0]}, Number of proteins without reference expression value: {num_un_mapped}',
            'yaxis_title':ylabel,
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )
    # return the results 
    return fig


def plot_gene_expression_vs_num_peptides(exp_count_table: pd.DataFrame, tissue_name: str,
    def_value: float = -1, plotting_kwargs: Dict[str,str] = {}, 
    xlabel: str = 'Number of peptides', ylabel: str = 'Expression value', 
    title: str= 'Peptides per protein Vs. Expression Level')->plt.Figure:
    """Plot the correlation between the gene expression and the num of peptids per protein using seaborn library 

    :param exp_count_table: A table that contain the number of peptides and the expresion value for each protein in the database
    :type exp_count_table: pd.DataFrame
    :param tissue_name: The name of the tissue 
    :type tissue_name: str
    :param def_value: The default value for proteins that could not be mapped to the expression database, defaults to -1
    :type def_value: float, optional
    :param plotting_kwargs: a dict object containing parameters for the sns.scatter function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param xlabel: the label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'Expression value'
    :type ylabel: str, optional
    :param title: The title of the figure, defaults to 'Peptides per protein Vs. Expression Level'
    :type title: str, optional
    """
    # First filter the DB for the non-mapped 
    df=exp_count_table.loc[exp_count_table.iloc[:,2]!=def_value,]
    # create a figure to plot to it 
    fig= plt.figure()
    ax=sns.scatterplot(x=df.iloc[:,1],y=df.iloc[:,2], **plotting_kwargs)
    ax.set_ylabel(ylabel+' in '+tissue_name)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    return fig

def plotly_gene_expression_vs_num_peptides(exp_count_table: pd.DataFrame, tissue_name: str,
    def_value: float = -1, xlabel: str = 'Number of peptides', ylabel: str = 'Expression value', 
    title: str= 'Peptides per protein Vs. Protein Expression Level')->plt.Figure:
    """Plot the correlation between the gene expression and the number of peptids per protein using plotly library. 

    :param exp_count_table: A table that contain the number of peptides and the expresion value for each protein in the database
    :type exp_count_table: pd.DataFrame
    :param tissue_name: The name of the tissue 
    :type tissue_name: str
    :param def_value: The default value for proteins that could not be mapped to the expression database, defaults to -1
    :type def_value: float, optional
    :param xlabel: The label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: The label on the y-axis., defaults to 'Expression value'
    :type ylabel: str, optional
    :param title: The title of the figure, defaults to 'Peptides per protein Vs. Protein Expression Level'
    :type title: str, optional
    """
    # First filter the DB for the non-mapped 
    df=exp_count_table.loc[exp_count_table.iloc[:,2]!=def_value,]
    # create a figure to plot to it 
    fig=px.scatter(df, x="Number_of_Peptides",y="Expression",
                labels={
                    'Number_of_Peptides':'Number of Peptides',
                    'Expression':'Expression Level'
                })
    # update the layout of the figure 
    fig.update_layout(
        {
            'title': ylabel+' in '+tissue_name, 
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )    
    return fig

def plot_num_protein_per_location(protein_loc: pd.DataFrame, 
    plotting_kwargs: Dict[str,str] = {}, 
    drop_unknown: bool = False, xlabel: str = 'Number of Proteins',
    ylabel: str = 'Compartment', title: str= 'Number of proteins per sub-cellular compartment'
    )->plt.Figure: 
    """plot the number of proteins per each sub-cellular compartment

    :param protein_loc: A table that contain the count of protein from each location.  
    :type protein_loc: pd.DataFrame
    :param plotting_kwargs: a dict object containing parameters for the sns.barplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param drop_unknown: whether or not to drop protein with unknown location, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: the label on the x-axis, defaults to 'Number of Proteins'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'Compartment'
    :type ylabel: str, optional
    :param title: the title of the figure, defaults to 'Number of proteins per sub-cellular compartment'
    :type title: str, optional

    """
    if drop_unknown:
        protein_loc=protein_loc.loc[protein_loc.iloc[:,0]!='UNK',]
    # create a figure 
    fig=plt.figure()
    ax=sns.barplot(y='Compartment',x='Counts',data=protein_loc, **plotting_kwargs)
    # set the labels 
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # return the results 
    return fig

def plotly_num_protein_per_location(protein_loc: pd.DataFrame, 
    drop_unknown: bool = False, xlabel: str = 'Number of Proteins',
    ylabel: str = 'Compartment', 
    title: str= 'Number of proteins per sub-cellular compartment'
    )->Figure: 
    """plot the number of proteins per each sub-cellular compartment
    
    :param protein_loc: A table that contain the count of protein from each location
    :type protein_loc: pd.DataFrame
    :param drop_unknown: whether or not to drop protein with unknown location, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: the label on the x-axis, defaults to 'Number of Proteins'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'Compartment'
    :type ylabel: str, optional
    :param title: the title of the figure, defaults to 'Number of proteins per sub-cellular compartment'
    :type title: str, optional
    """
    if drop_unknown:
        protein_loc=protein_loc.loc[protein_loc.iloc[:,0]!='UNK',]
    # create a figure 
    fig=px.bar(protein_loc,
        y="Compartment", x="Counts", color="Counts",
        color_continuous_scale='blackbody',
        labels={
            'Counts':xlabel, 
            'Compartment':ylabel 
        })
    # update the layout of the figure 
    fig.update_layout(
        {
            'title': title, 
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )    
    return fig

def plot_num_protein_per_go_term(protein2goTerm: pd.DataFrame, 
    tissue_name: str, plotting_kwargs: Dict[str,str] = {}, 
    drop_unknown: bool = False, xlabel: str = 'Number of Proteins',
    ylabel: str = 'Compartment', title: str= 'Number of proteins per sub-cellular compartment'
    )->plt.Figure: 
    """plot the number of proteins per each GO Term 

    :param protein2goTerm: A table that contain the count of proteins from each GO-Term 
    :type protein2goTerm: pd.DataFrame
    :param tissue_name: a dict object containing parameters for the sns.barplot function.
    :type tissue_name: str
    :param plotting_kwargs: a dict object containing parameters for the sns.barplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param drop_unknown: whether or not to drop protein with unknown location, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: the label on the x-axis, defaults to 'Number of Proteins'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'Compartment'
    :type ylabel: str, optional
    :param title: the title of the figure, defaults to 'Number of proteins per sub-cellular compartment'
    :type title: str, optional
    """
    if drop_unknown:
        protein2goTerm=protein2goTerm.loc[protein2goTerm.iloc[:,0]!='UNK',]
    # create a figure 
    fig=plt.figure()
    ax=sns.barplot(y='GO-Terms',x='Counts',data=protein2goTerm, **plotting_kwargs)
    # set the labels 
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # return the results 
    return fig

def plotly_num_protein_per_go_term(
    protein2goTerm: pd.DataFrame, 
    drop_unknown: bool = False, xlabel: str = 'Number of Proteins',
    ylabel: str = 'GO-Term', title: str= 'Number of proteins per GO-Term'
    )->Figure: 
    """ plot the number of proteins per each GO Term using plotly library 

    :param protein2goTerm: A table that contain the count of proteins from each GO-Term  
    :type protein2goTerm: pd.DataFrame
    :param drop_unknown: whether or not to drop protein with unknown location, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: the label on the x-axis, defaults to 'Number of Proteins'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'GO-Term'
    :type ylabel: str, optional
    :param title: the title of the figure, defaults to 'Number of proteins per GO-Term'
    :type title: str, optional
    """
    if drop_unknown:
        protein2goTerm=protein2goTerm.loc[protein2goTerm.iloc[:,0]!='UNK',]
    # create a figure 
    fig=px.bar(protein2goTerm,
        y="GO-Terms", x="Counts", color="Counts",
        color_continuous_scale='blackbody',
        labels={
            'Counts':xlabel, 
            'GO-Terms':ylabel 
        })
    # update the layout of the figure 
    fig.update_layout(
        {
            'title': title, 
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )    
    return fig

def plot_num_peptides_per_location(pep2loc: pd.DataFrame, plotting_kwargs: Dict[str,str] = {}, 
    drop_unknown: bool = False, xlabel: str = 'Number of peptides',
    ylabel: str = 'Compartment', 
    title: str= 'Number of peptides per sub-cellular compartment' ) -> plt.Figure:
    """plot the number of peptides obtained from each compartment using seaborn library.
    
    :param pep2loc: A table that contain the count of peptides from each compartment 
    :type pep2loc: pd.DataFrame
    :param plotting_kwargs: a dict object containing parameters for the sns.barplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param drop_unknown: whether or not to drop protein with unknown location, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: The label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: The label on the y-axis, defaults to 'Compartment'
    :type ylabel: str, optional
    :param title: The title of the figure, defaults to 'Number of peptides per sub-cellular compartment'
    :type title: str, optional
    """
    if drop_unknown:
        pep2loc=pep2loc.loc[pep2loc.iloc[:,0]!='UNK',]
    # create a figure 
    fig=plt.figure()
    ax=sns.barplot(y='Compartment',x='Counts',data=pep2loc, **plotting_kwargs)
    # set the labels 
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # return the results 
    return fig

def plotly_num_peptides_per_location(pep2loc: pd.DataFrame,
    drop_unknown: bool = False, xlabel: str = 'Number of peptides',
    ylabel: str = 'Compartment', 
    title: str= 'Number of peptides per sub-cellular compartment' ) -> plt.Figure:
    """  plot the number of peptides obtained from each compartment using plotly library

    :param pep2loc: A table that contain the count of peptides from each compartment 
    :type pep2loc: pd.DataFrame
    :param drop_unknown: whether or not to drop protein with unknown location, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: The label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: The label on the y-axis, defaults to 'Compartment'
    :type ylabel: str, optional
    :param title: The title of the figure, defaults to 'Number of peptides per sub-cellular compartment'
    :type title: str, optional
    """
    if drop_unknown:
        pep2loc=pep2loc.loc[pep2loc.iloc[:,0]!='UNK',]
    # create a figure 
    fig=px.bar(pep2loc,
        y="Compartment", x="Counts", color="Counts",
        color_continuous_scale='blackbody',
        labels={
            'Counts':xlabel, 
            'Compartment':ylabel 
        })
    # update the layout of the figure 
    fig.update_layout(
        {
            'title': title, 
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )    
    return fig

def plot_num_peptide_per_go_term(pep2goTerm: pd.DataFrame,
    plotting_kwargs: Dict[str,str] = {}, 
    drop_unknown: bool = False, xlabel: str = 'Number of peptides',
    ylabel: str = 'GO-Term', 
    title: str= 'Number of peptides per GO Term' ) -> plt.Figure:
    """plot the number of peptides obtained per Go-Term using seaborn library.

    :param pep2goTerm: A table that contain the count of peptides from each GO-Term 
    :type pep2goTerm: pd.DataFrame
    :param plotting_kwargs: a dict object containing parameters for the sns.barplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param drop_unknown: whether or not to drop peptide with unknown GO-term, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: the label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'GO-Term'
    :type ylabel: str, optional
    :param title: The title of the figure, defaults to 'Number of peptides per GO Term'
    :type title: str, optional
    """
    if drop_unknown:
        pep2goTerm=pep2goTerm.loc[pep2goTerm.iloc[:,0]!='UNK',]
    # create a figure 
    fig=plt.figure()
    ax=sns.barplot(y='GO-Terms',x='Counts',data=pep2goTerm, **plotting_kwargs)
    # set the labels 
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # return the results 
    return fig

def plotly_num_peptide_per_go_term(pep2goTerm: pd.DataFrame,
    drop_unknown: bool = False, xlabel: str = 'Number of peptides',
    ylabel: str = 'GO-Term', title: str= 'Number of peptides per GO Term' ) -> Figure:
    """ plot the number of peptides obtained per Go-Term  using plotly library.

    :param pep2goTerm: A table that contain the count of peptides from each GO-Term 
    :type pep2goTerm: pd.DataFrame
    :param drop_unknown: whether or not to drop peptide with unknown GO-term, defaults to False
    :type drop_unknown: bool, optional
    :param xlabel: the label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: the label on the y-axis, defaults to 'GO-Term'
    :type ylabel: str, optional
    :param title: the title of the figure, defaults to 'Number of peptides per GO Term'
    :type title: str, optional
    """
    if drop_unknown:
        pep2goTerm=pep2goTerm.loc[pep2goTerm.iloc[:,0]!='UNK',]
    # create a figure 
    fig=px.bar(pep2goTerm,
        y="GO-Terms", x="Counts", color="Counts",
        color_continuous_scale='blackbody',
        labels={
            'Counts':xlabel, 
            'GO-Terms':ylabel 
        })
    # update the layout of the figure 
    fig.update_layout(
        {
            'title': title, 
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )    
    return fig

def plot_num_peptides_per_organism(pep_per_org: pd.DataFrame, log_scale: bool =False, 
    plotting_kwargs: Dict[str,str] = {},  xlabel: str = 'Number of peptides',
    ylabel: str = 'Organism', title: str= 'Number of peptides per organism' ) -> plt.Figure:
    """plot the number of peptides per each organism inferred from the experiment using seaborn and matlotlib.

    :param pep_per_org: A table that contain the number of peptides belonging to each organism
    :type pep_per_org: pd.DataFrame
    :param log_scale: Whether or not to scale the number of peptides using a log scale, default is False, defaults to False
    :type log_scale: bool, optional
    :param plotting_kwargs: a dict object containing parameters for the sns.barplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param xlabel: the label on the x-axis, defaults to 'Number of peptides'
    :type xlabel: str, optional
    :param ylabel: The label on the y-axis, defaults to 'Organism'
    :type ylabel: str, optional
    :param title: The title of the figure, defaults to 'Number of peptides per organism'
    :type title: str, optional
    """
    # create a figure 
    fig=plt.figure()
    if log_scale:
        pep_per_org['Counts']=np.log10(pep_per_org['Counts'])
    ax=sns.barplot(y='Organisms',x='Counts',data=pep_per_org, **plotting_kwargs)
    # set the labels 
    if log_scale:
        xlabel="log10 of "+xlabel
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    # return the results 
    return fig

def plotly_num_peptides_per_organism(pep_per_org: pd.DataFrame, log_scale: bool =False, 
    xlabel: str = 'Number of Peptides', ylabel: str = 'Organism', 
    title: str= 'Number of peptides per organism') -> Figure:
    """plot the number of peptides per each organism inferred from the experiment using plotly library.
    
    :param pep_per_org: A table that contain the count of peptides from each organism
    :type pep_per_org: pd.DataFrame
    :param log_scale: Whether or not to scale the number of peptide using a log scale, defaults to False
    :type log_scale: bool, optional
    :param xlabel: The label on the x-axis, defaults to 'Number of Peptides'
    :type xlabel: str, optional
    :param ylabel: The label on the y-axis, defaults to 'Organism'
    :type ylabel: str, optional
    :param title: the title of the figure , defaults to 'Number of peptides per organism'
    :type title: str, optional
    """
    # normalize the counts by 
    if log_scale:
        pep_per_org['Counts']=np.log10(pep_per_org['Counts'])
    # plot the figure 
    fig=px.bar(pep_per_org, x='Counts', y='Organisms', 
                labels={
                    'Counts':xlabel, 
                    'Organisms':ylabel
                    }, 
                    title=title) 
    # update the background colors 
    fig.update_layout(
        {
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )
    return fig


def plot_change_in_presentation_between_experiment(
    change_in_presentation_array: np.ndarray, index_first: int, index_second:int, 
    plotting_kwargs: Dict[str,str] = {},
    title='Change in protein presentation',
    xlabel="Proteins", ylabel="magnitude of change in protein count") -> plt.Figure:
    """plot the change in protein presentation between two experiment 

    :param change_in_presentation_array: a 3D tensor of shape number of experiments by number of experiment by number of identified proteins. 
    :type change_in_presentation_array: np.ndarray
    :param index_first: [description]
    :type index_first: int
    :param index_second: the index of the first experiment in the tensor
    :type index_second: int
    :param plotting_kwargs: a dict object containing parameters for the sns.scatterplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param title: The title of the figure, defaults to 'Change in protein presentation'
    :type title: str, optional
    :param xlabel: The axis on the x-axis , defaults to "Proteins"
    :type xlabel: str, optional
    :param ylabel: The axis on the y-axis, defaults to "magnitude of change in protein count"
    :type ylabel: str, optional
    :raises ValueError: if the provided tensor is not of rank 3 
    :raises IndexError: if the provided indices are out of bound 
    """
    # check the correct index 
    if len(change_in_presentation_array.shape)!=3:
        raise ValueError(f"The provided tensor must of be of rank 3, current tensor has rank: {len(change_in_presentation_array.shape)}")    
    
    if index_first > change_in_presentation_array.shape[0]:
        raise IndexError(f"The provided index for the first experiment is out of bound. Number of elements along the first axis is: {change_in_presentation_array.shape[0]}")
    
    if index_second > change_in_presentation_array.shape[1]:
        raise IndexError(f"The provided index for the second experiment is out of bound. Number of elements along the second axis is: {change_in_presentation_array.shape[1]}")
    # get the tensor 
    change_tensor=change_in_presentation_array[index_first,index_second,:].reshape(-1)
    # sort the tensor 
    change_tensor.sort()
    change_tensor=change_tensor[::-1]
    # plot the tensor elements 
    fig=plt.figure()
    ax=sns.scatterplot(x=np.arange(len(change_tensor)),y= change_tensor, **plotting_kwargs)
    plt.hlines(y=np.median(change_tensor), xmin=0,xmax=len(change_tensor),color='red')
    # add the axis labels 
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # add the title
    ax.set_title(title)
    # return the results 
    return fig

def plot_experiment_set_counts(counts_table:pd.DataFrame, 
            log_scale: bool =False, plotting_kwargs: Dict[str,str] = {}) -> plt.figure:
    """visualize the number of peptides and number of peptides-per-organism per experiment.

    :param counts_table: a pandas dataframe that contain the count and the organism name.
    :type counts_table: pd.DataFrame
    :param log_scale: Normalize the peptide counts using log 10, defaults to False
    :type log_scale: bool, optional
    :param plotting_kwargs: a dict object containing parameters for the sns.catplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    """
    # Perform log 10 Normalization 
    if log_scale:
        counts_table['Counts'] = np.log10(counts_table['Counts'])
    # create a figure 
    fig=plt.figure()
    # plot the results 
    ax=sns.catplot(data=counts_table,x="Organisms",y="Counts",kind="bar",hue="Names")
    # change the y-label if log_scale is used 
    if log_scale: 
        ax.set_ylabels('Log10 of the Counts')
    # return the results 
    return fig

def plot_peptide_length_per_experiment(counts_table:pd.DataFrame,
                plotting_kwargs: Dict[str,str] = {}) -> plt.figure: 
    """visualize the peptide length distribution among the experiments defined in the set 
    

    :param counts_table: a pandas dataframe that contain the length of each peptide defined in the experiment along with \
    the experiment name 
    :type counts_table: pd.DataFrame
    :param plotting_kwargs: a dict object containing parameters for the sns.catplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    """
    # create a figure
    fig=plt.figure()
    # plot the results 
    ax= sns.catplot(x="Names",y="Peptide_length", data=counts_table, 
                    hue="Names",kind="box", **plotting_kwargs)
    # set the labels 
    ax.set_xlabels("Experiment IDs")
    ax.set_ylabels("Peptide length")
    # return the figure 
    return fig 

def plot_coverage_and_annotation(protein_coverage:Dict[str,np.ndarray],
                                temp_dir: str = ".",
                                figure_size: Tuple[int]=(7,8),
                                figure_dpi:int =600,
                                coverage_track_dict: Dict[str,Dict[str,Union[str,dict]]]={
                                  "coverage_dict":{"color":"grey",
                                                    "width":1.2}},
                                chains_track:bool =True,
                                chains_track_dict: Dict[str,Dict[str,Union[str,dict]]]={
                                  "track_label_dict":{"fontsize":6,"color":"black"},
                                "track_elements_dict":{"color":"magenta","capstyle":"butt"}},
                              domains_track:bool=True, 
                              domain_track_dict:Dict[str,Dict[str,Union[str,dict]]]={
                                  "track_label_dict":{"fontsize":6,"color":"black"},
                                   "track_element_names_dict":{"fontsize":4,"color":"black"},
                                   "track_elements_dict":{"color":"blue","capstyle":"butt"}
                                  },
                                transmembrane_track:bool=True, 
                                transmembrane_track_dict:Dict[str,Dict[str,Union[str,dict]]]={
                                  "track_label_dict":{"fontsize":6,"color":"black"},
                                   "track_element_names_dict":{"fontsize":4,"color":"black"},
                                   "track_elements_dict":{"color":"blue","capstyle":"butt"}
                                  },
                              modifications_track:bool=True,
                              modifications_track_dict:Dict[str,Dict[str,Union[str,dict]]]={
                                  "height_frac":0.5,
                                 "track_label_dict":{"fontsize":6,"color":"black"}},
                              glyco_track:bool=True,
                              glyco_track_dict: Dict[str,Dict[str,Union[str,dict]]]={
                                  "height_frac":0.5,
                              "track_label_dict":{"fontsize":6,"color":"black"}},
                              disulfide_track:bool =True,
                              disulfide_track_dict={
                                  "height_frac":0.5,
                              "track_label_dict":{"fontsize":6,"color":"black"}},
                              sequence_variants_track:bool =True,
                              sequence_variants_track_dict: Dict[str,Dict[str,Union[str,dict]]]={
                                  "height_frac":0.5,
                              "track_label_dict":{"fontsize":6,"color":"black"}},
                              )->plt.Figure:
    """The function plot the annotation track which summarizes all the known information about the protein and its associated peptides.
        
        Parameters
        ----------
        protein_coverage: Dict[str,np.ndarray]
            A dict that contain the uniprot accession of the protein as key and the protein coverage as value 

        temp_dir : str 
            The name of the temp directory to download the protein XML scheme to it. 

        figure_size: list/tuple, optional
            The figure size in inches. The default is 7 by 8 in which is 
            (17.8*20.32) cm. 
        
        figure_dpi: int, optional
            The resolution of the figure in dots-per-inch. The default is 600.
            
        coverage_track : bool, optional
            whether or not to plot the coverage track. The default is True.
        
        coverage_track_dict: dict, optional 
            The parameters of the function ``Annotator.add_coverage_track``.
            it is only used if coverage_track is set to True.
            The default is {"coverage_dict":{"color":"grey","width":1.2}}.
        
        chains_track : bool, optional
            whether or not to plot the chain track. The default is True.
        
        chains_track_dict: dict, optional 
            The parameters of the function ``Annotator.add_stacked_track``.
            it is only used if chains_track is set to True.
            The default is {"track_label_dict":{"fontsize":6,"color":"black"},
            "track_elements_dict":{"color":"magenta","capstyle":"butt"}}.
        
        domains_track : bool, optional
            whether or not to plot the domains track. The default is True.
            
        domains_track_dict: dict, optional 
            The parameters of the function ``Annotator.add_segmented_track``.
            it is only used if domains_track is set to True. 
            The default is {"track_label_dict":{"fontsize":6,"color":"black"},
            "track_element_names_dict":{"fontsize":4,"color":"black"},
            "track_elements_dict":{"color":"blue","capstyle":"butt"}}

        transmembrane_track : bool, optional
            whether or not to plot the transmembrane track. The default is True.
            
        transmembrane_track_dict: dict, optional 
            The parameters of the function ``Annotator.add_segmented_track``.
            it is only used if transmembrane_track is set to True. 
            The default is {"track_label_dict":{"fontsize":6,"color":"black"},
            "track_element_names_dict":{"fontsize":4,"color":"black"},
            "track_elements_dict":{"color":"blue","capstyle":"butt"}}
        
        modifications_track : bool, optional
            whether or not to plot the generic modification track. 
            The default is True.
            
        modifications_track_dict: dict, optional
            The parameters of the function ``Annotator.add_marked_positions_track``.
            it is only used if modifications_track is set to True. 
            The default is {"height_frac":0.5,
            "track_label_dict":{"fontsize":6,"color":"black"}}
            
        glyco_track : bool, optional
            whether or not to plot the generic glycosylation track.
            The default is True.
        
        glyco_track_dict: dict, optional 
            The parameters of the function ``Annotator.add_marked_positions_track``.
            it is only used if glyco_track is set to True. 
            The default is {"height_frac":0.5,
            "track_label_dict":{"fontsize":6,"color":"black"}}
            
        disulfide_track : bool, optional
            whether or not to plot the generic disulfide bond track. 
            The default is True.
        
        disulfide_track_dict: dict, optional
            The parameters of the function ``Annotator.add_marked_positions_track``.
            it is only used if disulfide_track is set to True. 
            The default is {"height_frac":0.5,
            "track_label_dict":{"fontsize":6,"color":"black"}}
            
        sequence_varients_track : bool, optional
            whether or not to plot the sequence varients track. 
            The default is True.
        
        sequence_varients_track_dict: dict, optional
            The parameters of the function ``Annotator.add_marked_positions_track``.
            it is only used if sequence_varients_track is set to True. 
            The default is {"height_frac":0.5,
            "track_label_dict":{"fontsize":6,"color":"black"}}.
        
        splice_varient_track : bool, optional
            whether or not to plot the solice varients track. The default is True.
        
        splice_varients_track_dict:
            The parameters of the function ``Annotator.add_stacked_track``.
            it is only used if chains_track is set to True.
            The default is {"track_label_dict":{"fontsize":6,"color":"black"},
            "track_elements_dict":{"color":"magenta","capstyle":"butt"}}.
        
        Examples
        --------
        In the following series of examples we are going to optimize the plotting of an annotator 
        track for the Syntenin-1 protein, protein accession O00560 which has 298 amino acids 
        
        ## first, let's simulate some coverage data  
        >>> sim_coverage: np.ndarray = np.concatenate([np.zeros(50),np.ones(50),np.zeros(198)])
        >>> panel_trial1: plt.Figure = plot_coverage_and_annotation({'O00560':sim_coverage})
        
        Ass seen from panel_trial1 the figure looks overcrowded an un-balanced. To polish the figure, let's first start 
        by adjusting the size of the axes' axis 

        >>> panel_trial2: plt.Figure = plot_coverage_and_annotation(
            {'O00560':sim_coverage},
            coverage_track_dict={
                "xlabel_dict":{
                                "fontsize":2
                                },
                "ylabel_dict":{
                                "fontsize":2
                                }
                            },
            chains_track_dict={
                "track_label_dict":{
                                "fontsize":2
                                    }
                            },
            domain_track_dict={
                "track_label_dict":{
                                "fontsize":2
                                    }
                            },
     modifications_track_dict={
                "track_label_dict":{
                                "fontsize":2
                                }
                            }
            )
        
        Okay, 2nd trial looks much better than the first trial, but we can polish it further by adjusting the font size
        and the size of domains and chains 

        >>> panel_trial3: plt.Figure = plot_coverage_and_annotation(
            {'O00560':sim_coverage},
            coverage_track_dict={
                "xlabel_dict":{
                                "fontsize":2
                                },
                "ylabel_dict":{
                                "fontsize":2
                                }
                            },
            chains_track_dict={
                "track_label_dict":{
                                "fontsize":2
                                   }, 
                "track_element_names_dict":{
                                        "fontsize":2
                                            }
                            },
            domain_track_dict={
                "track_label_dict":{
                                "fontsize":2
                                   },

                 "track_element_names_dict":{
                                        "fontsize":2
                                            }
                            },
             modifications_track_dict={
                "track_label_dict":{
                                "fontsize":2
                                   
                                   }
                            }
            )
        
        Okay, trial 3 looks much better than the previous tow trials, let's polish the graph further 
        by adjust the size of elments in the track as well as font-size 

        >>> panel_trial4: plt.Figure = plot_coverage_and_annotation(
            {'O00560':sim_coverage},
            coverage_track_dict={
                "xlabel_dict":{
                                "fontsize":3
                                },
                
                "ylabel_dict":{
                                "fontsize":3
                                }, 
                
                "coverage_dict":{
                    "color":"grey",
                    "width":1, 
                },
                "number_ticks":25, 

                            },
            chains_track_dict={
                "track_label_dict":{
                                "fontsize":3
                                   }, 
                "track_element_names_dict":{
                                        "fontsize":3
                                            },
                "track_elements_dict":{
                                    "color":"blue"
                                        }
                            },
            domain_track_dict={
                "track_label_dict":{
                                "fontsize":3
                                   },
                "track_element_names_dict":{
                                        "fontsize":3
                                        },
                "track_elements_dict":{
                                    "color":"lime"
                                      }
                            },
             modifications_track_dict={
                "track_label_dict":{
                                "fontsize":3
                                   },
                "marker_bar_dict": {
                                "color":"black",
                                "linestyles":"solid",
                                "linewidth":0.5,
                                "alpha":0.4
                                }, 
                "marker_dict":{
                                "s":2,
                                "color":"red"
                             }
                            }
            )
        Returns
        -------
        vispanel: matplotlib.figure.Figure
            The resulting panel.
    """
    # get the protein id 
    protin_id: str = list(protein_coverage.keys())[0]
    # adjust the protein shape 
    if len(protein_coverage[protin_id]) == 2: 
        protein_coverage[protin_id]=protein_coverage[protin_id].reshape(-1)
    # get all the information about the protein 
    protein_features=Features(protin_id,temp_dir)
    #create the panel
    panel=Annotator(protein_length=protein_coverage[protin_id].shape[0],
                       figure_size=figure_size, figure_dpi=figure_dpi)
    # add the base track   
    panel.add_coverage_track(protein_coverage[protin_id],
                  coverage_as_base=True,**coverage_track_dict)
    # add the chains
    if chains_track:
        chains=protein_features.get_chains() 
        if chains is not None:
            # update the chain dictionary names from chainId --> Name, more generic for the plotting function
            for chain_name in chains.keys(): 
                chains[chain_name]['Name']=chains[chain_name].pop('chainId')  
            # plot the results 
            panel.add_stacked_track(track_dict=chains,
                                    track_label="Chains",
                                    **chains_track_dict)  
        else:
            print("No chains are associated with this protein") 
    # add the domain track
    if domains_track:
        domains=protein_features.get_domains()
        if domains is not None:
            # add the domain name to the dict to make it as generic as possible --> more generic prinintg 
            for domain in domains.keys(): 
                domains[domain]['Name']=domain
            panel.add_segmented_track(track_dict=domains,
                                      track_label="Domains",
                                      **domain_track_dict)
        else:
            print("No domains are known in this protein")
    # add the transmembrane track
    if transmembrane_track:
        transmembrane_regions=protein_features.get_transmembrane_regions()
        if transmembrane_regions is not None:
            # prepear a dict containing the transmembrane regions, 
            trans_regions=dict()
            counter=1
            for region in transmembrane_regions:
                trans_regions['TM'+str(counter)]={
                    "Name":"TM"+str(counter),
                    "startIdx":region[0],
                    "endIdx":region[1]}
                counter+=1
            panel.add_segmented_track(track_dict=trans_regions,
                                      track_label="Transmembrane",
                                      **transmembrane_track_dict)
        else:
            print("No transmembrane regions are known in this protein")
    # add the modification track
    if modifications_track:
        modifications=get_PTMs_modifications_positions(protein_features)
        if len(modifications) !=0:
                panel.add_marked_positions_track(positions=modifications,
                                                 track_label="Modifications",
                                                **modifications_track_dict)
        else:
            print("No modifications are known in this protein")
    # add the glycozilation  track
    if glyco_track:
        glyco_positions=get_PTMs_glycosylation_positions(protein_features)
        if len(glyco_positions) !=0:
                panel.add_marked_positions_track(positions=glyco_positions,
                                                 track_label="Glycosylation",
                                                            **glyco_track_dict)
        else: 
            print("No glycosylation sites are known in this protein")
    # add disulfide_track
    if disulfide_track:
        sulfide_positions=get_PTMs_disuldfide_bonds(protein_features)
        if len(sulfide_positions) !=0:
            panel.add_marked_positions_track(positions=sulfide_positions,
                                                 track_label="Disulfide Bonds",
                                                            **disulfide_track_dict)
        else:
            print("No disulfide sites are known in this protein")
        
    # add sequence varient track:
    if sequence_variants_track: 
        sequence_variants_positions=get_sequence_variants_positions(protein_features)
        if len(sequence_variants_positions)!=0:
            panel.add_marked_positions_track(
                     positions=sequence_variants_positions,
                     track_label="Sequence Variants",
                     **sequence_variants_track_dict)     
        else:
            print("No sequence variants sites are known in this protein")
    plt.tight_layout() # adjust and scale the figure sizes 
    return panel.get_figure()

def plot_MDS_from_ic_coverage(distance_matrix: pd.DataFrame, 
    plotting_kwargs: Dict[str,str]={'c':'b'},
    title:str ='MDS plot using immunopeptidomic coverage as a distance metric', 
    MDS_params: Dict[str,str]={'random_state':93}, 
    plot_params: Dict[str,str]={}, 
    )->plt.Figure:
    """ plot MDS using immunopeptiomic coverage distance as a precomputed distance metric 

    :param distance_matrix: A pandas dataframe that contain the results of the distance matrix between each pair of experiments 
    :type distance_matrix: pd.DataFrame
    :param plotting_kwargs: a dict object containing parameters for the plt.plot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param title: the title of the MDS plot 
    :type title: str 
    :param MDS_params: a dict of parameters that will be forwarded to the function manifold.MDS, defaults to {'random_state':93}
    :type MDS_params: Dict[str:str]
    """
    # Get the matrix 
    dist=np.array(distance_matrix)
    # Normalize the matrix to values between 0 and 1
    dist/=np.amax(dist)
    # Compute the MDS 
    mds=manifold.MDS(n_components=2,dissimilarity="precomputed",**MDS_params)
    mds.fit(dist)
    # get the co-ordinate of the MDS 
    coords = mds.fit_transform(dist)
    names=distance_matrix.columns
    # plot the results 
    fig=plt.figure()
    plt.scatter(coords[:,0],coords[:,1],**plotting_kwargs)
    for label, x, y in zip(names, coords[:,0], coords[:,1]):
        plt.annotate(label, (x,y), xycoords = 'data')
    plt.xlabel('First dimension')
    plt.ylabel('Second dimension')
    plt.title(title)    
    # return the results 
    return fig 

def plot_num_peptides_per_protein_hist(num_pep_per_protein: pd.DataFrame,
            plotting_kwargs: Dict[str,str]={'c':'b'},
            title:str ='The Distribution of the number of peptides per protein') -> plt.Figure: 
    """plot a dist plot showing the distribution of the number of peptides observed from each protein. 

    :param num_pep_per_protein: a table containing the number of peptides observed from each protein inferred
    :type num_pep_per_protein: pd.DataFrame
    :param plotting_kwargs: a dict object containing parameters for the sns.histplot function, defaults to {}
    :type plotting_kwargs: Dict[str,str], optional
    :param title: the title of the MDS plot 
    :type title: str 
    :return: a dist plot constructed 
    :rtype: plt.Figure
    """
    fig=plt.figure()
    ax=sns.histplot(num_pep_per_protein.iloc[:,-1])
    ax.set_xlabel('Number of peptides per protein')
    ax.set_ylabel('Frequency')
    ax.title('The distribution of number of peptides per protein')  
    return fig 

def plot_MS_spectrum(spectrum: poms.pyopenms_2.MSSpectrum,
                    log_scale: bool=False,
                    color_specs: str ='black', 
                    color_base: str ='grey', 
                    grid_on: bool = True, 
                    spect_param: Dict[str,str]={},
                    base_param: Dict[str, str]={}, 
                    grid_param: Dict[str,str]={}, 
                    xlabel: str='M/Z',
                    ylabel: str= 'Intensity',
                    title: str= None
                    )->Union[plt.Figure]:
    """Plotting an Mass spectrometry spectrum 

    :param spectrum: a spectrum instance to plot it's peaks 
    :type spectrum: poms.pyopenms_2.MSSpectrum
    :param log_scale: a boolean flag of whether or not to normalize the y axis using a log scale, defaults to False
    :type log_scale: bool, optional
    :param color_specs: the color of the peaks in the figure, defaults to 'black'
    :type color_specs: str, optional
    :param color_base: the color of the baseline in the figure, defaults to 'grey'
    :type color_base: str, optional
    :param grid_on: , defaults to True
    :type grid_on: a boolean flag of whether or not to add a grid to the plot, optional
    :param spect_param: [description], defaults to {}
    :type spect_param: Dict[str,str], optional
    :param base_param: [description], defaults to {}
    :type base_param: Dict[str, str], optional
    :param grid_param: [description], defaults to {}
    :type grid_param: Dict[str,str], optional
    :param xlabel: [description], defaults to 'M/Z'
    :type xlabel: str, optional
    :param ylabel: [description], defaults to 'Intensity'
    :type ylabel: str, optional
    :param title: [description], defaults to None
    :type title: str, optional
    :return: [description]
    :rtype: Union[plt.Figure]
    """
    fig=plt.Figure()
    for x,y in zip(*spectrum.get_peaks()):
        if log_scale: 
            plt.vlines(x=x,ymax=np.log10(y),ymin=0,color=color_specs, 
                    **spect_param)
        else: 
            plt.vlines(x=x,ymax=y,ymin=0,color=color_specs,
                    **spect_param)
    plt.hlines(y=0, 
        xmin=min(spectrum.get_peaks()[0]),
        xmax=max(spectrum.get_peaks()[0]),
        color=color_base, **base_param)
    if grid_on: 
        plt.grid(**grid_param)
    ## add the X & Y label 
    plt.xlabel(xlabel)
    if log_scale: 
        plt.ylabel('log10'+ylabel)
    ## Add the title 
    if title is not None: 
        plt.title(title)
    return fig

def plot_chord_diagram_among_set(exp_set,
            level:str='protein', filename='results_fig.svg',
            fig_params:Dict[str,str]={
                'node_cmap':'PuBu',
                'edge_cmap':'Category20',
                'height':700,
                'width':700,
                'title':"Shared Protein Representation among Individuals",
                'edge_alpha':0.5,
                'edge_line_width':1,
                'label_text_color':'blue'
            }):
    """plot a chord diagram showing the overlap among a group of immunopeptiomics experiments 
    
    :param exp_set: an experimental set instance containing
    :type exp_set: An ExperimentSet
    :param level: the level of similiarity, currently, the library support protein and peptide levels 
    :type level: str 
    :param filename: the name of the file to save the results as a SVG figure, defaults to results_fig.svg in the current working directory
    :type filename: str 
    :param fig_params: a dict of parameters to optimize the output of the figure
    :type fig_params: dict 
    """
    ## Compute a dataframe with the experiment name 1 experiment name 2, # number sahred protein
    if level=='protein':
        overlap_table=exp_set.compute_protein_overlap_matrix()
    else:
        overlap_table=exp_set.compute_peptide_overlap_matrix()
    ## Unrol the matrix into the long form 
    overlap_table=overlap_table.stack().rename_index(['Individal_A','Individal_B']).reset_index()
    ## Remove the diagonal elements 
    overlap_table=overlap_table[overlap_table.loc['Individal_A']==overlap_table.loc['Individal_B']]
    ## Get the unique experimental names 
    unique_experimental_name=list(exp_set.get_experiments().keys())
    unique_experimental_ds=hv.Dataset(pd.DataFrame({'Exp_id':unique_experimental_name}))
    ## Generate a chord diagram among the nodes of the Chord  
    chord_diagram=hv.Chord((overlap_table,unique_experimental_ds))
    ## cutamize the output of the figure 
    chord_diagram=chord_diagram.opts(node_color="Exp_id",edge_color="Individal_A",**fig_params)
    ## generate the outout as SVG 
    export_svg(hv.render(chord_diagram),filename=filename)
    return

def plotly_goea_results(res_df:pd.DataFrame, focus:str='CC')->go.Figure:
    """ plot a bubble plot of the generate plotly figure.

    Args:
        res_df (pd.DataFrame): the input table as computed by the GOEngine 
        focus (str, optional): The class of GO term to focus on, can be any of 'CC' for cellular_component,
        MF for molecular functions or 'MF' for molecular functions. Defaults to 'CC'.

    Returns:
        go.Figure: a plotly figure showing the results 
    """
    if focus.upper() not in ['CC','MF','BP']:
        raise ValueError(f"The value of focus: {focus} is not supported, currently, supporting: {','.join(['CC','MF','BP'])}")
    filtered_df=res_df.loc[res_df.NS==focus].copy(deep=True)
    filtered_df['p_fdr_bh']= -1*np.log10([float(num) for num in filtered_df['p_fdr_bh']])
    fig=px.scatter(filtered_df,
            y="name", x="p_fdr_bh",
	        size="study_count",
            labels={'p_fdr_bh':'-log10(p-val)','name':''},
            template='plotly_white')
    return fig




