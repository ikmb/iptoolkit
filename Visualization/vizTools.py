#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: The module contain visualization functions the can be used to plot the results obtained from the 
datastructures API or from the analysis functions defined in the Analysis Module.   
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
from scipy.stats import ttest_ind
from statannot import add_stat_annotation
import plotly.express as px 
from plotly.graph_objects import Figure
import plotly.graph_objects as go 
from plotly import tools 
import math
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
def plot_overlap_heatmap(results_df:pd.DataFrame, plotting_kwargs: PlottingKeywards={})->sns.matrix.ClusterGrid:
    """
    @brief: plot a user provided dataframe as a cluster heatmap using seaborn library
    @param: results_df: a pandas dataframe table that hold the overlapping number.
    """
    fig=sns.clustermap(results_df, **plotting_kwargs)
    return fig
# define the plotly version of the library 
def plotly_overlap_heatmap(results_df: pd.DataFrame)->Figure:
    """
    @brief: plot a user provided dataframe as a heatmap using plotly library
    @param: results_df: a pandas dataframe table that hold the overlapping number.
    """
    return px.imshow(results_df) 

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
                            alpha: float = 0.9, title=" Parired protein representation" ) -> plt.Figure:  
    """ 
    @brief: compare the representation between two experiments using matplotlib library.
    @param: protein_one_repr: a dict object containing the legand of the first protein along with its mapped array
    @param: protein_two_repr: a dict object containing the legand of the second protein along with its mapped aray  
    @param: color_first: the color of representation for the first protein 
    @param: color_second: the color of the second protein 
    @param: alpha: the transparency of the figure. 
    @param: title: the title of the figure. 
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
    """
    @brief: compare the peptide coverage for the same protein under different conditions using the same protein using plotly library.
    @param: protein_one_repr: a dict object containing the legand of the first protein along with its mapped array
    @param: protein_two_repr: a dict object containing the legand of the second protein along with its mapped array  
    @param: title: the title of the figure. 
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

def plotly_multi_traced_coverage_representation(proteins, 
                title: str= "Protein Coverage Across  ") -> Figure : 
    """
    @brief: plot a multi-traced representation for the same protein accross 
    @param: proteins: a dict object containing for each protein the corresponding mapped array. 
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
    fig.update_layout(title=title+" "+str(trace_counter)+" Conditions",
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
    """
    @brief: plot the peptide coverage for a given protein. 
    @param:  mapped_protein: a NumPy array with shape of 1 by protein length or shape protein-length
    @param: col: the color of the coverage respresentation
    @param: prot_name: the default protein name
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

def plotly_protein_coverage(mapped_protein: np.ndarray, prot_name: str = None)->plt.Figure:
    """
    @brief: plot the peptide coverage for a given protein. 
    @param:  mapped_protein: a numpy array with shape of 1 by protein length or shape protein-length
    @param: prot_name: the default protein name
    """
    if len(mapped_protein.shape)==2:
        mapped_protein=mapped_protein.reshape(-1) 
    # translate the array into a dataframe 
    df_res= pd.DataFrame(mapped_protein)
    df_res.columns=['Coverage']
    # plot the figure 
    fig=px.area(x,y="Coverage",
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
                                    background_color: str ='black', low: str ='red', high: str ='blue')->None: 
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
    """
    @brief: visualize a histogram of the eluted peptide length using seaborn library. 
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

def plotly_peptide_length_dist(pep_length: List[int],
                         plotting_kwargs: Dict[str,str]={}, 
                         x_label: str = 'Peptide Length',
                         y_label: str = 'Counts',
                         title: str = 'Peptide Length distribution'):
    """
    @brief: visualize a histogram of the eluted peptide length using plotly library 
    @param: pep_length: a list of integer containing the peptides' lengths
    @param: x_label: the label of the x-axis 
    @param: y_label: the label of the y-axis 
    @param: title: the title of the figure
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

def plotly_num_peptides_per_parent(nums_table: pd.DataFrame,
                        num_prot: int = -1, 
                        plotting_kwargs: Dict[str,str]={}, 
                        x_label: str = 'Number of peptides',
                        y_label: str = 'Protein ID',
                        title: str = 'Number of peptides per protein'):
    """
    @brief: visualize a histogram of the the number of peptides per each inferred protein.  
    @param: nums_table: a pandas dataframe containing number of peptides identified from each protein. 
    @param: num_prot, the number of protein to show relative to the first element, for example, the first 10, 20 etc.
    If the default value of -1 is used then all protein will be plotted, however, this might lead to a crowded figure.
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
            'plot_bgcolor':'rgba(0,0,0,0)',
            'paper_bgcolor':'rgba(0,0,0,0)'
        }
    )
    return fig

def plot_parent_protein_expression_in_tissue(expression_table: pd.DataFrame, 
    ref_expression: pd.DataFrame , tissue_name: str, sampling_num: int = 10,
    plotting_kwargs: Dict[str,str]={'orient':'v'}, def_value: float = -1,
    ylabel: str = 'Normalized Expression') -> plt.Figure:
    """
    @brief: plot the parent protein expression in tissue relative a sampled collection of non-presented data using seaborn library.
    @param: expression_table: The protein expression table which contains the expresion value for each parent protein
    @param: ref_expression: The reference expression of the tissue under investigation. 
    @param: sampling_num: the number of times to sample from the non-prsenter. 
    @param: tissue_name: The name of the tissue 
    @param: def_value: The default value for proteins that could not be mapped to the expression database 
    @param: plotting_kwargs: a dict object containing parameters for the sns.violinplot function.
    @param: ylabel: the label on the y-axis. 
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
    title=f't-test score: {ttest_res.statistic}, P-val: {ttest_res.pvalue }'
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
    """
    @brief: plot the parent protein expression in tissue relative a sampled collection of non-presented data using plotly library.
    @param: expression_table: The protein expression table which contains the expresion value for each parent proteins. 
    @param: ref_expression: The reference expression of the tissue under investigation. 
    @param: sampling_num: the number of times to sample from the non-prsenter. 
    @param: tissue_name: The name of the tissue. 
    @param: def_value: The default value for proteins that could not be mapped to the expression database. 
    @param: ylabel: the label on the y-axis. 
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
    """
    @brief: plot the correlation between the gene expression and the num of peptids per protein using seaborn library 
    @param: exp_count_table: A table that contain the number of peptides and the expresion value for each protein in the database. 
    @param: tissue_name: The name of the tissue 
    @param: def_value: The default value for proteins that could not be mapped to the expression database 
    @param: plotting_kwargs: a dict object containing parameters for the sns.scatter function.
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis.
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

def plotly_gene_expression_vs_num_peptides(exp_count_table: pd.DataFrame, tissue_name: str,
    def_value: float = -1, xlabel: str = 'Number of peptides', ylabel: str = 'Expression value', 
    title: str= 'Peptides per protein Vs. Protein Expression Level')->plt.Figure:
    """
    @brief: plot the correlation between the gene expression and the number of peptids per protein using plotly library. 
    @param: exp_count_table: A table that contain the number of peptides and the expresion value for each protein in the database. 
    @param: tissue_name: The name of the tissue 
    @param: def_value: The default value for proteins that could not be mapped to the expression database 
    @param: plotting_kwargs: a dict object containing parameters for the sns.scatter function.
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of proteins per each sub-cellular compartment
    @param: protein_loc: A table that contain the count of protein from each location.  
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of proteins per each sub-cellular compartment
    @param: protein_loc: A table that contain the count of protein from each location.  
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: title: the title of the figure.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of proteins per each GO Term 
    @param: protein2goTerm: A table that contain the count of proteins from each GO-Term  
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of proteins per each GO Term using plotly library. 
    @param: protein2goTerm: A table that contain the count of proteins from each GO-Term  
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: title: the title of the figure.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of peptides obtained from each compartment using seaborn library. 
    @param: pep2code: A table that contain the count of peptides from each  location 
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis. 
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of peptides obtained from each compartment using plotly library. 
    @param: pep2code: A table that contain the count of peptides from each  location 
    @param: drop_unknown: whether or not to drop protein with unknown location. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: title: the title of the figure.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of peptides obtained per Go-Term  using matplotlib library.
    @param: pep2goTerm: A table that contain the count of peptides from each GO-Term 
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: drop_unknown: whether or not to drop peptide with unknown GO-term. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of peptides obtained per Go-Term  using plotly library.
    @param: pep2goTerm: A table that contain the count of peptides from each GO-Term 
    @param: drop_unknown: whether or not to drop peptide with unknown GO-term. Default is False.  
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis.
    @param: title: the title of the figure.
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
    """
    @brief: plot the number of peptides per each organism inferred from the experiment using seaborn and matlotlib.
    @param: pep_per_org: A table that contain the number of peptides belonging to each organism.
    @param: log_scale: Whether or not to scale the number of peptides using a log scale, default is False.
    @param: plotting_kwargs: a dict object containing parameters for the sns.barplot function.
    @param: ylabel: the label on the y-axis. 
    @param: xlabel: the label on the x-axis. 
    @param: title: the title of the figure.
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
    title: str= 'Number of peptides per organism', 
    paper_bg_color: str ='rgba(0,0,0,0)') -> Figure:
    """
    @brief: plot the number of peptides per each organism inferred from the experiment using plotly library. 
    @param: pep_per_org: A table that contain the count of peptides from each organism.
    @param: log_scale: Whether or not to scale the number of peptide using a log scale, default is False.
    @param: xlabel: the label on the x-axis 
    @param: ylabel: the label on the y-axis 
    @param: title: the title of the figure.
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
    """ 
    @brief: plot the change in protein presentation between two experiment 
    @param: change_in_presentation_array: a 3D tensor of shape number of experiments by 
    number of experiment by number of identified proteins. 
    @param: index_first: the index of the first experiment in the tensor. 
    @param: index_second: the index of the second experiment in the tensor. 
    @param: plotting_kwargs: a dict object containing parameters for the sns.scatterplot function.
    @param: title: The title of the figure
    @param: xlabel: The x-axis label of the figure 
    @param: ylabel: the y-axis label of the figure
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
    """
    @brief: visualize the number of peptides and number of peptides-per-organism per experiment. 
    @param: counts_table: a pandas dataframe that contain the count, organism name and the count.
    @param: plotting_kwargs: a dict object containing parameters for the sns.catplot function.  
    @param: log_scale: Normalize the peptide counts one log 10, default is False 
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
    """
    @brief: visualize the peptide length distribution among the experiments defined in the set 
    @param: counts_table: a pandas dataframe that contain the length of each peptide defined in the experiment along with 
    the experiment name 
    @param: plotting_kwargs: a dict object containing parameters for the sns.catplot function.  
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