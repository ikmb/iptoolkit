#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: The module contain a collection of analysis function that can be used by the methods of the analysis
module. 
@version: 0.0.1
"""
# load the modules: 
import numpy as np 
import pandas as pd 
import subprocess as sp 
import os
from Bio.PDB import PDBList
from Bio.motifs.meme import Motif 
from typing import List, Callable, Dict, Set
from IPTK.IO import MEMEInterface as memeIF
from IPTK.IO import OutFunctions as out_func
from IPTK.DataStructure.Experiment import Experiment
from IPTK.Utils.UtilityFunction import check_peptide_made_of_std_20_aa
from IPTK.Utils.Mapping import map_from_uniprot_gene
from scipy.stats import pearsonr
# define some types 
Peptides=List[str]
Proteins=List[str]
## Define the peptide overlap 
def get_binnary_peptide_overlap(exp1:Experiment, exp2:Experiment)->Peptides:
    """
    @brief: compare the peptide overlap between two experimental objects.
    @param: exp1: an instance of class Experiment.
    @param: exp2: an instance of class Experiment.
    @return: a list of peptides that have been identified in both experiments.    
    """ 
    peptide_one=exp1.get_peptides()
    peptide_two=exp2.get_peptides()
    return list(peptide_one.intersection(peptide_two))
     
def get_binnary_protein_overlap(exp1:Experiment, exp2:Experiment)->Proteins:
    """
    @brief: compare the protein overlap between two experimental objects.
    @param: exp1: an instance of class Experiment.
    @param: exp2: an instance of class Experiment.
    @return: a list of proteins that have been identified in both experiments.   
    """ 
    protein_one=exp1.get_proteins()
    protein_two=exp2.get_proteins()
    return list(protein_one.intersection(protein_two))

def compute_binary_distance(peptides: List[str],dist_func:Callable)->np.ndarray:
    """
    @brief: compare the distance between every pair of peptides in a collection of peptides. 
    @param: peptides: a collection of peptides sequences.
    @param: dist_func: function to compute the distance between each pair of peptides. 
    @note: make sure that the dist_function is suitable with the peptides which might have different lengths.
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

def get_sequence_motif(peptides:Peptides, keep_temp:bool=False, 
                       temp_dir: str ="./TEMP_DIR", verbose: bool = False, 
                       meme_params:Dict[str,str]={}
                       )->None:
    """
    @brief: compute the sequences motif from a collection of peptide sequences using meme software.
    @param: peptides: a list of string containing the peptide sequences 
    @param: keep_temp: whether or not to keep the temp file written to the temp_directory, default is False. 
    @param: temp_dir: the temp directory to write temp-files to it
    @param: verbose: whether or not to print the output of the motif discovery tool to the stdout.  
    @param: meme_params: a dict object that contain meme controlling parameters.
    @see: IPTK.IO.MEMEInterface for more details 
    """
    # check the meme is installed 
    if not memeIF.is_meme_callable():
        raise FileNotFoundError(f"meme is either not installed or not part of the PATH ==> {os.environ['PATH']}")
    # check the temp directory
    try:
        os.mkdir(temp_dir)
    except FileExistsError: 
        pass 
    # check that amino acids have the same length 
    if len(set([len(pep) for pep in peptides])) !=1:
        raise ValueError('The provided peptides MUST have the same length!')
    # check that only the 20 classical amino acids are in the motif object 
    curated_sequences: List[str] =  [check_peptide_made_of_std_20_aa(pep) for pep in peptides]
    curated_sequences = [elem for elem in curated_sequences if elem !='']
    # write the sequences 
    outfile_seq_file=os.path.join(temp_dir,'TEMP_FASTA_SEQ.fasta')
    out_func.write_auto_named_peptide_to_fasta(peptides,outfile_seq_file)
    # call meme to compute the motif(s)
    meme_results_dir=os.path.join(temp_dir,'TEMP_MEME_RES')
    memeIF.call_meme(outfile_seq_file,
                output_dir=meme_results_dir, 
                verbose=verbose, **meme_params)
    print(f"MEME has finished execution. Results can be found at: {meme_results_dir}")
    
def download_structure_file(pdb_id:str)->None:
    """
    @brief: Download PDB/mmCIF file containing the pbd_id from PDB using BioPython library 
    @param: pdb_id: the protein id in protein databank 
    """
    pdb_list=PDBList()
    pdb_list.retrieve_pdb_file(pdb_id)
    return 

def compute_expression_correlation(exp1:Experiment,exp2:Experiment)->float:
    """
    @brief: compute the correlation in the gene expression between two experiments by constructing a union
    of all the proteins expressed in the first and second experiments, extract the gene expression 
    of these genes and then compute the correlation using scipy stat module. 
    @param: exp1: The first experimental object 
    @param: exp2: The second experimental object 
    """
    # get the expression tables 
    protein_exp1: Set[str] = set(exp1.get_proteins())
    protein_exp2: Set[str] = set(exp2.get_proteins())
    unique_proteins = list(protein_exp1.union(protein_exp2))
    # get the gene id 
    prot2Ense: pd.DataFrame = map_from_uniprot_gene(unique_proteins)
    # allocate lists to hold the results 
    gene_expression_exp1: List[float] = []
    gene_expression_exp2: List[float] = [] 
    # get the expression from experiment one 
    for prot in unique_proteins:
        temp_df: pd.DataFrame = prot2Ense.loc[prot2Ense.iloc[:,0]==prot]
        if temp_df.shape[0] == 1: # we got only one match 
            gene_id: str = temp_df['Gene-ID'].tolist()[0]
            try: 
                gene_expression_exp1.append(exp1.get_tissue().get_expression_profile().get_gene_id_expression(gene_id=gene_id))
            except KeyError:
                gene_expression_exp1.append(-1)
        else:
            temp_gene_expression: List[float] = []
            for gene in temp_df.iloc[:,1].tolist():
                try: 
                    temp_gene_expression.append(
                        exp1.get_tissue().get_expression_profile().get_gene_id_expression(gene_id=gene))
                except KeyError: 
                    temp_gene_expression.append(-1)
            # filter the temp_genes for default value 
            temp_gene_process: List[str] = [elem for elem in temp_gene_expression if elem != -1]
            # add the gene expression as the average if all values have been filtered 
            if len(temp_gene_process) ==0:
                gene_expression_exp1.append(-1)
            else: 
                gene_expression_exp1.append(np.mean(temp_gene_process))
    # ge the expression from exp2: 
    for prot in unique_proteins:
        temp_df: pd.DataFrame = prot2Ense.loc[prot2Ense.iloc[:,0]==prot]
        if temp_df.shape[0] == 1: # we got only one match 
            gene_id: str = temp_df['Gene-ID'].tolist()[0]
            try: 
                gene_expression_exp2.append(exp2.get_tissue().get_expression_profile().get_gene_id_expression(gene_id=gene_id))
            except KeyError:
                gene_expression_exp2.append(-1)
        else:
            temp_gene_expression: List[float] = []
            for gene in temp_df.iloc[:,1].tolist():
                try: 
                    temp_gene_expression.append(
                        exp2.get_tissue().get_expression_profile().get_gene_id_expression(gene_id=gene))
                except KeyError: 
                    temp_gene_expression.append(-1)
            # filter the temp_genes for default value 
            temp_gene_process: List[str] = [elem for elem in temp_gene_expression if elem != -1]
            # add the gene expression as the average 
            if len(temp_gene_process) ==0:
                gene_expression_exp2.append(-1)
            else: 
                gene_expression_exp2.append(np.mean(temp_gene_process))
    # compute construct a dataframe 
    temp_paired_exp_df: pd.DataFrame = pd.DataFrame({
        'exp2': gene_expression_exp1, 
        'exp1': gene_expression_exp2
    })
    # filter the un-mapped from exp1
    temp_paired_exp_df= temp_paired_exp_df.loc[temp_paired_exp_df.iloc[:,0]!=-1,]
    # filter the unmapped from exp2 
    temp_paired_exp_df= temp_paired_exp_df.loc[temp_paired_exp_df.iloc[:,1]!=-1,]
    
    # compute the correlation 
    return pearsonr(temp_paired_exp_df.iloc[:,0],temp_paired_exp_df.iloc[:,1])[0]

def compute_change_in_protein_representation(mapped_prot_cond1: np.ndarray, 
    mapped_prot_cond2: np.ndarray)->float:
    """
    @brief: compute the change in the protein representation between two conditions, by computing 
    the difference in the area under the curve, AUC.
    @param: mapped_prot_cond1: a mapped protein instance containing the protein coverage in the first condition
    @param: mapped_prot_cond2: a mapped protein instance containing the protein coverage in the second condition  
    """
    # un-roll the arrays to 1D arrays.
    if len(mapped_prot_cond1.shape)==2:
        mapped_prot_cond1=mapped_prot_cond1.reshape(-1)
    
    if len(mapped_prot_cond2.shape)==2:
        mapped_prot_cond2=mapped_prot_cond2.reshape(-1)
    
    # assert that the protein have same length 
    if len(mapped_prot_cond1)!=len(mapped_prot_cond2):
        raise ValueError(f'The provided proteins are of different length, found: {len(mapped_prot_cond1)} and {len(mapped_prot_cond1)}') 
    
    # compute the AUC 
    sum_array_one: float = np.sum(mapped_prot_cond1)
    sum_array_two: float = np.sum(mapped_prot_cond2)
    # compute the differences 
    difference: float = np.abs(sum_array_one-sum_array_two)
    # return the results 
    return difference

def compute_difference_in_representation(mapped_prot_cond1: np.ndarray, 
    mapped_prot_cond2: np.ndarray) -> np.ndarray:
    """
    @brief: return the difference in the representation of protein between two proteins.
    @param: mapped_prot_cond1: a mapped protein instance containing the protein coverage in the first condition
    @param: mapped_prot_cond2: a mapped protein instance containing the protein coverage in the second condition  
    """
     # un-roll the arrays to 1D arrays.
    if len(mapped_prot_cond1.shape)==2:
        mapped_prot_cond1=mapped_prot_cond1.reshape(-1)
    
    if len(mapped_prot_cond2.shape)==2:
        mapped_prot_cond2=mapped_prot_cond2.reshape(-1)
    
    # assert that the protein have same length 
    if len(mapped_prot_cond1)!=len(mapped_prot_cond2):
        raise ValueError(f'The provided proteins are of different length, found: {len(mapped_prot_cond1)} and {len(mapped_prot_cond1)}') 

    return mapped_prot_cond1-mapped_prot_cond2
























