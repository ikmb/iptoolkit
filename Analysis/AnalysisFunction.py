#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: CPU impelementation of the analysis functions
@version: 0.0.1
@date: 19.08.2020  
"""
# load the modules: 
import numpy as np 
import pandas as pd 
import subprocess as sp 
from typing import List, Callable
from IPTK.IO import MEMEInterface as memeIF
from IPTK.IO import OutFunctions as out_func
from Bio.PDB import PDBList
from IPTK.DataStructure.Experiment import Experiment
# define some types 
Peptides=List[str]
Proteins=List[str]
## Define the peptide overlap 
def get_binnary_peptide_overlap(exp1:Experiment, exp2:Experiment)->Peptides:
    """
    @brief: compare the peptide overlap between two experimental objects
    @param: exp1: an instance of class Experiment
    @param: exp2: an instance of class Experiment
    @return: a list of peptides that have been identified in both experiment   
    """ 
    peptide_one=exp1.get_peptides()
    peptide_two=exp2.get_peptides()
    return peptide_one.intersection(peptide_two)
     
def get_binnary_protein_overlap(exp1:Experiment, exp2:Experiment)->Proteins:
    """
    @brief: compare the protein overlap between two experimental objects
    @param: exp1: an instance of class Experiment
    @param: exp2: an instance of class Experiment
    @return: a list of peptides that have been identified in both experiment   
    """ 
    protein_one=exp1.get_proteins()
    protein_two=exp2.get_proteins()
    return list(protein_one.intersection(protein_two))

def compute_binary_distance(peptides: List[str],dist_func:Callable)->np.ndarray:
    """
    @brief: compare the distance between every pair of peptides in a collection of peptides. 
    @param: peptides: a collection of peptides sequence
    @param: dist_func: function to compute the distance between each pair of peptides 
    @note: make sure that the dist_function is suitable with the peptides which might have different lengths
    """
    num_peptides=len(peptides)
    distance_matrix=np.zeros(shape=(num_peptides,num_peptides))
    # compute the pair wise distance 
    for raw_idx in range(num_peptides):
        for col_idx in range(num_peptides):
            try:
                distance_matrix[raw_idx,col_idx]=dist_func(peptides[raw_idx],peptides[col_idx])
            except Exception as exp:
                raise RuntimeError(f'While computing the pairwise distance between {peptides[raw_idx]} and {peptides[col_idx]} The following error was observed: {exp}')
    # return the results 
    return distance_matrix

def get_sequence_motif(method_params,peptides:Peptides, method:str='meme', keep_temp:bool=False,
                       temp_dir: str ="./TEMP_DIR", verbose: bool = False, 
                       ):
    """
    @brief: compute the sequences motif from a collection of peptide sequences.
    @details: The function use the subprocess module to make a system call to the function 
    specified using the @param: method. and call the function @@xx@@ defined in @@yy@@ to parse the method
    output to a specific parameters.  
    @param: peptides: a list of string containing the peptide sequences 
    @param: method: the method to call to compute the motif, can be any of MEME,
    NNAlgin, mhc x , case insenstive.   
    @param: keep_temp: whether or not to keep the temp file written to the temp_directory, default is False. 
    @param: temp_dir: the temp directory to write temp-directory to it
    @param: verbose: whether or not to print the output of the motif discovery to the screen.  
    @param: method_params: a dict object that contain method specific parameters 
    """
    if method=='meme':
        # check the meme is installed 
        if not memeIF.is_meme_callable():
            raise FileNotFoundError(memeIF.MEME_IS_NOT_INSTALLED_ERROR)
        # check the temp directory exists
        try:
            os.mkdir(temp_dir)
        except FileExistsError: 
            pass 
        # write the sequences 
        outfile_seq_file=os.join.path(temp_dir,'TEMP_FASTA_SEQ.fasta')
        out_func.write_auto_named_peptide_to_fasta(peptides,outfile_seq_file)
        # call meme to compute the motif(s)
        meme_results_dir=os.path.join(temp_dir,'TEMP_MEME_RES')
        memeIF.call_meme(outfile_seq_file,
                        output_dir=os.path.join(temp_dir,'TEMP_MEME_RES'), 
                        verbose=verbose, **method_params)
        # parse the results 
        motifs=memeIF.parse_meme_results(meme_results_dir)
        return motfis 
    pass

def download_structure_file(pdb_id:str)->None:
    """
    @brief:  downlaod PDB/mmCIF file containgthe pbd_id from PDB using BioPython Module 
    @param: pdb_id: the protein id in protein databank 
    """
    pdb_list=PDBList()
    pdb_list.retrieve_pdb_file(pdb_id)
    return 

