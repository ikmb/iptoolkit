#!/usr/bin/env Python 
"""Utility functions that are used through the library 
"""
# load the models 
from Bio import SeqIO
from IPTK.Utils.Types import FastaSet
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib 
import string
import random
import pickle 
import urllib
from typing import List, Dict
# define the functions 
def pad_mapped_proteins(list_array: List[np.ndarray],
                pre_pad:bool =True, padding_char: int =-1)->np.ndarray:
    """ Pad the provided list of array into a 2D tensor of shape number of arrays by maxlength. 

    :param list_array: A list of NumPy arrays where each array is a mapped_protein array, \
    the expected shape of these arrays is 1 by protein length.
    :type list_array: List[np.ndarray]
    :param pre_pad: pre or post padding of shorter array in the list_array. Defaults to True, which mean prepadding
    :type pre_pad: bool, optional
    :param padding_char: The padding char, defaults to -1
    :type padding_char: int, optional
    :return: A 2D tensor of shape number of arrays by maxlength. 
    :rtype: np.ndarray
    """
    # reshape the resulting arrays
    for idx in range(len(list_array)):
        list_array[idx].reshape(1,-1) # reshape te array into shape 1* number of elements 
    # getting the max 
    max_len: int =max([elem.shape[1] for elem in list_array])
    # compute the padding distance 
    paddling_lens: List[int]=[max_len-elem.shape[1] for elem in list_array]
    # making arrays to hold padding_arrays
    padding_arrays: List[np.ndarray]= []
    # generate the padding distance 
    for idx in range(len(list_array)):
        padding_arrays.append(np.array([padding_char]*paddling_lens[idx]).reshape(1,-1))
    # allocate a list to hold the results 
    resulting_arrays: List[np.ndarray] = []
    # fuse the arrays 
    for idx in range(len(list_array)):
        if pre_pad:
            resulting_arrays.append(
                np.concatenate([padding_arrays[idx], list_array[idx]],axis=1)
            )
        else: 
            resulting_arrays.append(
                np.concatenate([list_array[idx], padding_arrays[idx] ],axis=1)
            )
    # concat the results array 
    results_array: np.ndarray = np.concatenate(resulting_arrays,axis=0)
    return results_array

def generate_random_name(name_length: int)->str:
    """
    :param name_length: Generate a random ASCII based string
    :type name_length: int
    :return: [description]
    :rtype: str
    """
    chars=[char for char in string.ascii_uppercase] # get ASCII upper chars 
    chars.extend([str(elem) for elem in list(range(10))]) # extend the vocab 
    return ''.join([random.choice(chars) for _ in range(name_length)]) # generate the name 

def append_to_calling_string(param:str, def_value, cur_val, calling_string:str, is_flag: bool = False)->str:
    """ help function that take a calling string, a parameter, a default value and current value \
    if the parameter does not equal its default value the function append the parameter with its current \ 
    value to the calling string adding a space before the calling_string. 

    :param param: The name of the parameter that will be append to the calling string 
    :type param: str
    :param def_value: The default value for the parameter 
    :type def_value: [type]
    :param cur_val: The current value for the parameter
    :type cur_val: [type]
    :param calling_string: The calling string in which the parameter and the current value might be appended to it 
    :type calling_string: str
    :param is_flag: If the parameter is a control flag, i.e. a boolean switch, it append the parameter to the calling string without associating a value to it , defaults to False
    :type is_flag: bool, optional
    :return: the updated version of the calling string 
    :rtype: str
    """
    if is_flag: 
        if def_value != cur_val:
            calling_string += ' -'+param
        return calling_string
    if def_value != cur_val: 
        calling_string+=" -"+param+" "+str(cur_val)
    return calling_string

def generate_random_protein_mapping(protein_len: int , max_coverage: int) -> np.ndarray:
    """Generate a NumPy array with shape of 1 by protein_len where the elements in the array 
    is a random integer between zero &  max_coverage. 

    :param protein_len: The length of the protein 
    :type protein_len: int
    :param max_coverage: The maximum peptide coverage at each position 
    :type max_coverage: int
    :return: a NumPy array containing a simulated protein coverage 
    :rtype: np.ndarray
    """
    return np.random.randint(low=0, high= max_coverage, size=(protein_len,))

def generate_color_scale(color_ranges: int )-> matplotlib.colors.LinearSegmentedColormap: 
    """generate a color gradient with number of steps equal to color_ranges -1 
    
    :param color_ranges:  the number of colors in the range
    :type color_ranges: int
    :return: A color gradient palette
    :rtype: matplotlib.colors.LinearSegmentedColormap
    """
    return plt.cm.get_cmap('hsv',color_ranges)

def simulate_protein_representation(num_conditions : int , protein_len: int ,
         protein_coverage: int )->Dict[str, np.ndarray]:
    """ Simulate protein peptide coverage under-different conditions
    
    :param num_conditions: The number of condition to simulate 
    :type num_conditions: [type]
    :param protein_len: The length of the protein  
    :type protein_len: [type]
    :param protein_coverage: The maximum protein coverage 
    :type protein_coverage: [type]
    :return: a dict of length num_conditions containing the condition index and a simulated protein array   
    :rtype: Dict[str, np.ndarray]
    """
    sim_res=dict()
    color_gradient=generate_color_scale(num_conditions)
    for idx in range(num_conditions): 
        sim_res['COND_'+str(idx)]={'name': 'COND_'+str(idx),
                                   'mapped_protein':generate_random_protein_mapping(protein_len,protein_coverage),
                                   'color':color_gradient(idx)
                                   }
    return sim_res

def simulate_protein_binary_represention(num_conditions: int, protein_length: int): 
    """
    :param num_conditions: The number of conditions to simulate 
    :type num_conditions: int
    :param protein_length: The Length of the protein  
    :type protein_length: int
    :return: A 2D matrix of shape protein_length by number of conditions, where each element can be either zero or 1.
    :rtype: np.ndarray
    """
    return np.random.randint(low=0, high=2, size=(protein_length,num_conditions)).astype(np.float64)
 
def save_3d_figure(outpath: str, fig2save: plt.Figure) ->None:
    """write a pickled version of the a 3D figure so it can be loaded later for more interactive analysis
    
    :param outpath: The output path of the writer function 
    :type outpath: str
    :param fig2save: The figure to save to the output file
    :type fig2save: plt.Figure
    :raises IOError: In case writing the file failed 
    """
    try: 
        with open(outpath,'wb') as writer_buf:
            pickle.dump(fig2save,writer_buf)
    except Exception as exp:
        raise IOError(f'While writing the figure to the provided path the following error was encountered: {exp}')

def load_3d_figure(file_path: str) ->plt.Figure: 
    """
    :param file_path: Load a pickled 3D figure from the provided path 
    :type file_path: str
    :raises IOError: The path of the pickled figure. 
    :return: a matplotlib figure 
    :rtype: plt.Figure
    """
    try:
        with open(file_path, 'rb') as reader_buf:
            fig=pickle.load(reader_buf)
        return fig 
    except Exception as exp: 
        raise IOError(f'While loading your figure, the following error was encountered: {exp}')   

def build_sequence_table(sequence_dict:Dict[str,str])->pd.DataFrame:
    """construct a sequences database from a sequences dict object  
    
    :param sequence_dict: a dict that contain the protein ids as keys and sequences as values. 
    :type sequence_dict: Dict[str,str]
    :return: pandas dataframe that contain the protein ID and the associated protein sequence 
    :rtype: pd.DataFrame
    """
    return pd.DataFrame(sequence_dict)

def get_idx_peptide_in_sequence_table(sequence_table:pd.DataFrame, peptide:str)->List[str]:
    """check the sequences table if the provided peptide is locate in one of its sequences and returns 
    a list of protein identifiers containing the identifier of the hit proteins.

    :param sequence_table:  pandas dataframe that contain the protein ID and the associated protein sequence 
    :type sequence_table: pd.DataFrame
    :param peptide: The peptide sequence to query the protein with 
    :type peptide: str
    :return: A list of protein identifiers containing the identifier of the hit proteins
    :rtype: List[str]
    """
    if sequence_table.columns != ['Sequences']: 
        sequence_table.columns=['Sequences']
    return sequence_table.loc[sequence_table['Sequences'].str.contains(peptide)].index.tolist()

def check_peptide_made_of_std_20_aa(peptide:str)->str:
    """Check if the peptide is made of the standard 20 amino acids, if this is the case, 
    it return the peptide sequence, otherwise it return an empty string
    
    :param peptide: a peptide sequence to check its composition
    :type peptide: str
    :return: True, if the peptide is made of the standard 20 amino acids, False otherwise. 
    :rtype: str
    """
    amino_acids=['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    for amino_acid in peptide: 
        if amino_acid not in amino_acids: 
            return ''
    return peptide

def get_experiment_summary(ident_table: pd.DataFrame) -> pd.DataFrame:
    """takes as an input an identification table and return a summary table containing the count of unique peptides,\
    unique proteins, maximum peptide length, minmum peptide length, median and mean peptide length 

    :param ident_table: the identification table as returned by one of the parser functions defined in the IO modules 
    :type ident_table: pd.DataFrame
    :return: The summary table 
    :rtype: pd.DataFrame
    """
    # compute the statsitics 
    num_unique_peptides: int = len(set(ident_table.peptide))
    num_unique_proteins: int = len(set(ident_table.protein))
    peptide_len: np.ndarray= np.array(ident_table.end_index) - np.array(ident_table.start_index)
    max_len: int = np.max(peptide_len) 
    min_len: int = np.min(peptide_len)
    medain_len: float = np.median(peptide_len)
    mean_len: float = np.mean(peptide_len)
    # generate the data frame 
    res: pd.DataFrame = pd.DataFrame({
        'num_unique_peptides':[num_unique_peptides], 
        'num_unique_proteins': [num_unique_proteins], 
        'max_len': [max_len], 
        'min_len':[min_len],
        'median_len':[medain_len], 
        'mean_len': [mean_len]
    })
    # return the results 
    return res 

def combine_summary(child_dfs: List[pd.DataFrame], root_df: pd.DataFrame = None) -> pd.DataFrame: 
    """combine multiple summaray dataframes into one dataframe 

    :param child_dfs: a list of summary dataframes to conctinate into one 
    :type child_dfs: List[pd.DataFrame]
    :param root_df: a dataframe to append the child dataframe to its tail, defaults to None
    :type root_df: pd.DataFrame, optional
    :return: a dataframe containing the root and the child dataframes 
    :rtype: pd.DataFrame
    """
    if root_df is None: 
        root_df= pd.DataFrame(columns=['num_unique_peptides','num_unique_proteins','max_len','min_len','median_len','mean_len'])
    # loop over the child dataframes
    for child_df in child_dfs: 
        root_df=pd.concat([root_df, child_df], axis=0)
    # return the results 
    return root_df