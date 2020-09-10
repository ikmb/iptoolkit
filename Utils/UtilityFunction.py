#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: utility functions that are used through the library 
@version: 0.0.1
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
    """
    @brief pad the provided list of array into a 2D tensor of shape
    number of arrays by maxlength. 
    @param: list_array: a list of numpy arrays where each array is a mapped_protein array, 
    the expected shape of these arrays is 1 by protein length.
    @param: pre_pad: pre or post padding of shorter array in the library.Default is pre-padding.
    @param: padding_char: The padding char, Default is -1. 
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
    @brief generate a random ASCII based string
    @param name_length the number of charachter in the name 
    """
    chars=[char for char in string.ascii_uppercase] # get ASCII upper chars 
    chars.extend([str(elem) for elem in list(range(10))]) # extend the vocab 
    return ''.join([random.choice(chars) for _ in range(name_length)]) # generate the name 

def append_to_calling_string(param:str, def_value, cur_val, calling_string:str, is_flag: bool = False)->str:
    """
    @brief: a help function that take a calling string, a parameter, a default value and current value 
    if the parameter does not equal its default value the function append the parameter with its current 
    value to the calling string adding a space before the calling_string 
    @param: param: the name of the parameter that will be append to the calling string 
    @param: def_value: the default value for the parameter 
    @param: cur_val: the current value for the parameter
    @param: calling_string: the calling string in which the parameter and the current value might be appended to it 
    @param: is_flag: if the parameter is a control flag, i.e. a boolean switch, it append the parameter to the calling 
    string without associating a value to it 
    """
    if is_flag: 
        if def_value != cur_val:
            calling_string += ' -'+param
        return calling_string
    if def_value != cur_val: 
        calling_string+=" -"+param+" "+str(cur_val)
    return calling_string

def generate_random_protein_mapping(protein_len: int , max_coverage: int) -> np.ndarray:
    """
    @brief: generate a numpy array with shape of 1 by protein_len where the elements in the array 
    is a random integer between zero &  max_coverage 
    @param: protein_len: the protein length 
    @param: max_coverage: the maximum coverage at each position in the amino acid position  
    """
    return np.random.randint(low=0, high= max_coverage, size=(protein_len,))

def generate_color_scale(color_ranges: int )-> matplotlib.colors.LinearSegmentedColormap: 
    """
    @brief: generate a color gradient with number of steps equal to color_ranges -1 
    @param: color_ranges: the number of colors in the range, 
    """
    return plt.cm.get_cmap('hsv',color_ranges)

def simulate_protein_representation(num_conditions, protein_len, protein_coverage):
    """
    @brief: simulate protein peptide coverage under-different conditions
    @param: num_conditions: The number of condition to simulate 
    @param: protein_len: The length of the protein  
    @param: protein_coverage: The maximum protein coverage 
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
    @brief: Return a 2D matrix of shape protein_length by number of conditions, where each element can be either 
    zero.
    @num_conditions: The number of conditions to simulate 
    @protein_length: The Length of the protein  
    """
    return np.random.randint(low=0, high=2, size=(protein_length,num_conditions)).astype(np.float64)
 
def save_3d_figure(outpath: str, fig2save: plt.Figure) ->None:
    """
    @brief: write a picklized version of the a 3D figure so it can be loaded later for more interactive analysis
    @param: outpath: The output path of the writer function 
    @param: fig2save: The figure to save to the output file
    """    
    try: 
        with open(outpath,'wb') as writer_buf:
            pickle.dump(fig2save,writer_buf)
    except Exception as exp:
        raise IOError(f'While writing the figure to the provided path the following error was encountered: {exp}')

def load_3d_figure(file_path: str) ->plt.Figure: 
    """
    @brief: load a picklized 3D figure from thr provided path 
    @param: file_path: the path of the pickilized figure. 
    """    
    try:
        with open(file_path, 'rb') as reader_buf:
            fig=pickle.load(reader_buf)
        return fig 
    except Exception as exp: 
        raise IOError(f'While loading your figure, the following error was encountered: {exp}')   

def build_sequence_table(sequence_dict:dict)->pd.DataFrame:
    """
    @brief construct a sequecnes database from sequecnes dict object  
    """
    return pd.DataFrame(sequence_dict)

def get_idx_peptide_in_sequence_table(sequence_table:pd.DataFrame, peptide:str):
    """
    @brief check the sequences table if the provided peptide is locate in one of its sequences and returns 
    a list of protein identifiers containing the identifier of the hit proteins. 
    """
    if sequence_table.columns != ['Sequences']: 
        sequence_table.columns=['Sequences']
    return sequence_table.loc[sequence_table['Sequences'].str.contains(peptide)].index.tolist()

def check_peptide_made_of_std_20_aa(peptide:str)->str:
    """
    @brief: check if the peptide is made of the standard 20 amino acids, if this is the case, 
    it return the peptide sequence, otherwise it return an empty string
    @param: peptide: the sequence of the peptide to check 
    """
    amino_acids=['A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    for amino_acid in peptide: 
        if amino_acid not in amino_acids: 
            return ''
    return peptide