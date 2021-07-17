#!/usr/bin/env Python 
"""The module contains a collection of analysis functions that can be used by the methods of
the classes defined in the IPTK.Classes module. 
"""
# load the modules: 
from multiprocessing import Value
import numpy as np 
import pandas as pd 
import subprocess as sp 
import os
from Bio.PDB import PDBList
from numba import jit 
from Bio.motifs.meme import Motif 
from typing import List, Callable, Dict, Set, Union
from IPTK.IO import MEMEInterface as memeIF
from IPTK.IO import OutFunctions as out_func
from IPTK.Classes.Experiment import Experiment
from IPTK.Utils.UtilityFunction import check_peptide_made_of_std_20_aa
from IPTK.Utils.Mapping import map_from_uniprot_gene
from scipy.stats import pearsonr
from IPTK.Classes.Features import Features
# define some types 
Peptides=List[str]
Proteins=List[str]
## Define the peptide overlap 
def get_binnary_peptide_overlap(exp1:Experiment, exp2:Experiment)->Peptides:
    """compare the peptide overlap between two experimental objects.

    :param exp1: an instance of class Experiment.
    :type exp1: Experiment
    :param exp2: an instance of class Experiment.
    :type exp2: Experiment
    :return: a list of peptides that have been identified in both experiments.
    :rtype: Peptides
    """
    peptide_one=exp1.get_peptides()
    peptide_two=exp2.get_peptides()
    return list(peptide_one.intersection(peptide_two))
     
def get_binnary_protein_overlap(exp1:Experiment, exp2:Experiment)->Proteins:
    """compare the protein overlap between two experimental objects.

    :param exp1: an instance of class Experiment.
    :type exp1: Experiment
    :param exp2: an instance of class Experiment.
    :type exp2: Experiment
    :return: a list of proteins that have been identified or inferred in both experiments. 
    :rtype: Proteins
    """
    protein_one=exp1.get_proteins()
    protein_two=exp2.get_proteins()
    return list(protein_one.intersection(protein_two))


@jit(nopython=True,nogil=True,cache=True,parallel=True)
def compute_binary_distance_axf(peptides: List[str], dist_func:Callable)->np.ndarray:
    """Compare the distance between every pair of peptides in a collection of peptides. 
    
    :param peptides: a collection of peptides.
    :type peptides: List[str]
    :param dist_func: a function to compute the distance between each pair of peptides. 
    :type dist_func: Callable, that have been generated using Numba JIT.
    :raises RuntimeError: make sure that the dist_function is suitable with respect to the input peptides. For example, peptides which might have different lengths.
    :return: the distance between each pair of peptides in the provided list of peptides
    :rtype: np.ndarray
    """
    num_peptides=len(peptides)
    distance_matrix=np.zeros(shape=(num_peptides,num_peptides))
    # compute the pair wise distance 
    for raw_idx in range(num_peptides):
        for col_idx in range(num_peptides):
            distance_matrix[raw_idx,col_idx]=dist_func(peptides[raw_idx],peptides[col_idx])
    # return the results 
    return distance_matrix

def compute_binary_distance(peptides: List[str], dist_func:Callable)->np.ndarray:
    """compare the distance between every pair of peptides in a collection of peptides. 
    
    :param peptides: a collection of peptides.
    :type peptides: List[str]
    :param dist_func: a function to compute the distance between each pair of peptides. 
    :type dist_func: Callable
    :raises RuntimeError: make sure that the dist_function is suitable with respect to the input peptides. For example, peptides which might have different lengths.
    :return: the distance between each pair of peptides in the provided list of peptides
    :rtype: np.ndarray
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

def get_sequence_motif(peptides:Peptides, 
                       temp_dir: str ="./TEMP_DIR", verbose: bool = False, 
                       meme_params:Dict[str,str]={}
                       )->None:
    """compute the sequences motif from a collection of peptide sequences using meme software.
    
    :param peptides: a list of string containing the peptide sequences 
    :type peptides: Peptides
    :param temp_dir: the temp directory to write temp-files to it, defaults to "./TEMP_DIR"
    :type temp_dir: str, optional
    :param verbose: whether or not to print the output of the motif discovery tool to the stdout, defaults to False
    :type verbose: bool, optional
    :param meme_params: a dict object that contains meme controlling parameters, defaults to {}
    :type meme_params: Dict[str,str], optional
    :raises FileNotFoundError: incase meme is not installed or could not be found in the path!
    :raises ValueError: Incase the peptides have different length! 
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
    """Download PDB/mmCIF file with a user provided identifer from PDB using BioPython library 

    :param pdb_id: the protein id in protein data bank 
    :type pdb_id: str
    """
    pdb_list=PDBList()
    pdb_list.retrieve_pdb_file(pdb_id)
    return 

def compute_expression_correlation(exp1:Experiment,exp2:Experiment)->float:
    """compute the correlation in the gene expression between two experiments by constructing a union
    of all the proteins expressed in the first and second experiments, extract the gene expression 
    of these genes and then compute the correlation using SciPy stat module. 
    
    :param exp1: The first experimental object 
    :type exp1: Experiment
    :param exp2: The second experimental object 
    :type exp2: Experiment
    :return: the correlation in gene expression of the proteins inferred in the provided pair of experiment
    :rtype: float
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
    """Compute the change in the protein representation between two conditions, by computing 
    the differences in the area under the curve, AUC.

    :param mapped_prot_cond1: a mapped protein instance containing the protein coverage in the first condition
    :type mapped_prot_cond1: np.ndarray
    :param mapped_prot_cond2: a mapped protein instance containing the protein coverage in the second condition  
    :type mapped_prot_cond2: np.ndarray
    :raises ValueError: if the provided pair of proteins is of different length 
    :return: the difference in the area under the coverage curve between the two experiments. 
    :rtype: float
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
    """return the difference in the representation of a protein between two conditions
    by substracting the coverage of the first protein from the second protein.

    :param mapped_prot_cond1: a mapped protein instance containing the protein coverage in the first condition
    :type mapped_prot_cond1: np.ndarray
    :param mapped_prot_cond2: a mapped protein instance containing the protein coverage in the second condition  
    :type mapped_prot_cond2: np.ndarray
    :return: an array that shows the difference in coverage between the two proteins at each amino acid position. 
    :rtype: np.ndarray
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

def get_PTMs_modifications_positions(protein_feature: Features) ->List[int] : 
    r""" return a list of integers containing the position of PTMs in the protein.

    :param protein_feature: A protein feature instance containing all protein features 
    :type protein_feature: Features
    :return: A list that contains the position of the generic modifications in \
            the protein. If no modifications is known the function returns an empty list.
    :rtype: List[int]
    """
    containing=[]
    modification_dict=protein_feature.get_PTMs_modifications()
    if modification_dict==None: 
        return containing
    for key in modification_dict.keys():
            containing.append(modification_dict[key]["startIdx"])
    return containing

 
def get_PTMs_glycosylation_positions(protein_feature: Features) ->List[int]:
    r""" return a list of integers containing glycosylation position in the protein.

    :param protein_feature: a protein feature instance containing all protein features 
    :type protein_feature: Features
    :return: a list that contains the position of the generic glycosylation in \
            the protein. If no glycosylation site(s) are/is known in the protein \
            the function returns an empty list.
    :rtype: List[int]
    """
    glyco_positions=[]
    glyco_dict=protein_feature.get_PTMs_glycosylation()
    if glyco_dict==None: 
        return glyco_positions
    for key in glyco_dict.keys():
        glyco_positions.append(glyco_dict[key]["startIdx"])
    return glyco_positions 


def get_PTMs_disuldfide_bonds(protein_feature)->List[int]:
    """  return a list of integers containing dissulfide bound position in the protein.

    :param protein_feature: a protein feature instance containing all protein features 
    :type protein_feature: Features
    :return: A list that contains the position of the known disulfide bonds in \
            the protein. If no disulfide bond(s) is/are known in the protein \
            the function returns an empty list.
    :rtype: List[int]
    """
    disulfide_positions=[]
    disulfide_dict=protein_feature.get_disulfide_bonds()
    if disulfide_dict==None:
         return disulfide_positions
    for key in disulfide_dict.keys():
        disulfide_positions.append(disulfide_dict[key]["startIdx"])
    return disulfide_positions

def get_sequence_variants_positions(protein_feature)->List[int]:
    """return a list of integers containing known sequence variants positions in the protein.

    :param protein_feature: a protein feature instance containing all protein features 
    :type protein_feature: Features
    :return: A list that contains the position of the known sequence variants in 
            the protein. If no sequence variant/variants that is/are known in the protein,\
            the function returns an empty list.
    :rtype: List[int]
    """
    seqVar_positions=[]
    seqVar_dict=protein_feature.get_sequence_variants()
    if seqVar_dict==None:
        return seqVar_positions
    for key in seqVar_dict.keys():
        seqVar_positions.append(seqVar_dict[key]["startIdx"])
    return seqVar_positions

def get_chain_positions(protein_feature)->List[List[int]]:
    """ return a list of lists containing the start and end position of known chains in the protein 

    :param protein_feature: a protein feature instance containing all protein features 
    :type protein_feature: Features
    :return:  a list of list that contains the position of the known chain/chains \
            in the protein. each list has two elements which are the start and \
            the end position of the position lists. \
            If no chains are known in  the protein, \
            the function returns an empty list. \
    :rtype: List[List[int]]
    """
    chain_positions=[]
    chains_dict=protein_feature.get_chains()
    if chains_dict==None: 
        return chain_positions
    for key in chains_dict.keys():
            chain_positions.append([chains_dict[key]["startIdx"],
                                    chains_dict[key]["endIdx"]])
    return chain_positions


def get_domains_positions(protein_feature)->List[List[int]]:
    """ return a list of lists containing the start and end position of known domains in the protein.  

    :param protein_feature: a protein feature instance containing all protein features 
    :type protein_feature: Features
    :return: a list of list that contains the position of the known domain/domains \
            in the protein. each list has two elements which are the start and \
            the end position of the domain. \
            If there are no domain/domains that is/are known in the protein, \
            the function returns an empty list. \
    :rtype: List[List[int]]
    """
    domains_positions=[]
    domains_dict=protein_feature.get_domains()
    if domains_dict==None: 
        return domains_positions
    for key in domains_dict.keys():
        domains_positions.append([domains_dict[key]["startIdx"],
                                    domains_dict[key]["endIdx"]])
    return domains_positions

def get_splice_variants_positions(protein_feature)->List[List[int]]:
    """ return a list of lists containing the start and end position of known splice variants in the protein.  

    :param protein_feature: a protein feature instance containing all protein features 
    :type protein_feature: Features
    :return:  a list of list that contains the position of the known splice variant/variants \ 
            in the protein. each list has two elements which are the start and \
            the end position of the splice variants. \
            If no splice variant/variants that is/are known in the protein, \
            the function returns an empty list. \
    :rtype: List[List[int]]
    """
    splice_variant_positions=[]
    variants=protein_feature.get_splice_variants()
    if variants==None: 
        return splice_variant_positions
    for key in variants.keys():
        splice_variant_positions.append([variants[key]["startIdx"],
                                    variants[key]["endIdx"]])
    return splice_variant_positions

def compute_ic_distance_protein(protein_id:str, experiment_set,
                                mode="restrictive") -> pd.DataFrame: 
    """compute the immunopeptidomic coverage distance between a group of experiments \
        in the experiment set using one protein defined in protein_id. 

    :param protein_id: the uniprot or the identifier of the protein in Experiment objects
    :type protein_id: str
    :param experiment_set: An experiment set object containing all the experiments 
    :type experiment_set: ExperimentSet
    :param mode: mode of calculations, if restrictive, the protein defined by protein_id MUST\
         be present in every experiment. If is not defined an error will be through. However,\
         incase the mode is "permissive" the absent protein will be treated as an array of zeros.\
         it defaults to "restrictive".
    :type mode: str, optional
    :return: a square distance matrix contain the differences in immunopeptidomics coverage between each pair of experiments\
         in the set of experiments. 
    :rtype: pd.DataFrame
    """
    # extract the mapping array: 
    exps=experiment_set.get_experiments()
    mapped_arrays=dict()
    ## loop over each experiment and extract the experiment 
    for exp_id in exps.keys(): 
        try: 
            mapped_arrays[exp_id]=exps[exp_id].get_mapped_protein(protein_id)
        except KeyError as exp: 
            if mode=="permissive": 
                mapped_arrays[exp_id]=-1
            else: 
                raise KeyError(f"protein id: {protein_id} is not defined in experiment: {exp_id}")
    ## compute the pair wise distance 
    ic_dist: np.ndarray = np.zeros((len(mapped_arrays),len(mapped_arrays)))
    ## get the names of the arrays  
    row_names=col_names=list(mapped_arrays.keys())
    ## loop over the array element to fill 
    for row_idx in range(len(row_names)): 
        for col_idx in range(len(col_names)): 
            # handling the 4 possible cases that can be encountered during the execution 
            if isinstance(mapped_arrays[row_names[row_idx]], int): 
                if isinstance(mapped_arrays[row_names[col_idx]], int):
                    ic_dist[row_idx,col_idx]=0
                else: 
                    ic_dist[row_idx,col_idx]=np.sum(np.abs(mapped_arrays[row_names[col_idx]]))
            else: 
                if isinstance(mapped_arrays[row_names[col_idx]], int):
                    ic_dist[row_idx,col_idx]=np.sum(np.abs(mapped_arrays[row_names[row_idx]]))
                else: 
                    ic_dist[row_idx,col_idx]=np.sum(np.abs(mapped_arrays[row_names[row_idx]]-mapped_arrays[row_names[col_idx]]))
    ## construct a data frame out of the results 
    results_df: pd.DataFrame = pd.DataFrame(mapped_arrays)   
    results_df.columns=col_names
    results_df.index=row_names
    ## return the results 
    return results_df

def compute_ic_distance_experiments(experiment_set,mode="restrictive")->pd.DataFrame:
    """compute the immunopeptidomic coverage distance between a group of experiments \
        using all proteins in the intersection or the union, depending on the mode, of the set\
        of proteins defined in the set. 

    :param experiment_set: An experiment set object containing all the experiments 
    :type experiment_set: ExperimentSet
    :param mode: mode of calculations, if restrictive the proteins defined by protein_id MUST\
        be defined in every experiment, if it is not defined and error will be through. However,\
        incase it is "permissive" the absent protein will be treated as an array of zeros.\
        it defaults to "restrictive".
    :type mode: str, optional
    :return: A square distance matrix containing the avarage difference in the immunopeptidomic coverage\
        between each pair of experiments in the experiment set.  
    :rtype: pd.DataFrame
    """
    ## GETTING the set of proteins depending on the mode
    if mode=="restrictive": 
        protin_set=experiment_set.get_proteins_present_in_all()
    else: 
        protein_set=[]
        for exp_id in experiment_set.get_experiments().keys(): 
            protein_set.extend(experiment_set[exp_id].get_proteins())
    protin_set=list(set(protin_set))
    ## getting the dict of experiments 
    exps=experiment_set.get_experiments()
    ## computing the distance tensor Tensor 
    distance_tensor = np.zeros((len(exps),len(exps), len(protin_set)))
    ## looping over all the elements in the tensor 
    row_names=col_names= list(exps.keys())
    ## FILL the tensor  
    # ----------------
    for row_idx in range(len(row_names)):
        for col_idx in range(len(col_names)): 
            for protein_idx in range(len(protin_set)):
                ## obtain the mapped array of the rows
                try: 
                    mapped_array_er_pp=exps[row_names[row_idx]].get_mapped_protein(protin_set[protein_idx])
                except KeyError: 
                    mapped_array_er_pp=-1
                ## obtain the mapped array of the columns 
                try: 
                    mapped_array_ec_pp= exps[col_names[col_idx]].get_mapped_protein(protin_set[protein_idx]) 
                except KeyError: 
                    mapped_array_er_pp=-1
                ## fill the tensor using the mapped array information 
                if isinstance(mapped_array_er_pp, int): # i.e. set to -1 --> not present   
                    if isinstance(mapped_array_ec_pp, int): # i.e. set to -1 --> not present  
                        distance_tensor[row_idx,col_idx,protein_idx]=0
                    else: 
                        distance_tensor[row_idx,col_idx,protein_idx]=np.sum(np.abs(mapped_array_ec_pp))
                else: 
                    if isinstance(mapped_array_ec_pp, int):
                        distance_tensor[row_idx,col_idx,protein_idx]=np.sum(np.abs(mapped_array_er_pp))
                    else:
                        distance_tensor[row_idx,col_idx,protein_idx]=np.sum(np.abs(mapped_array_er_pp-mapped_array_ec_pp)) 
    ## COMPUTE the mean across the protein set: 
    distance_matrix= np.mean(distance_tensor,-1) 
    ## CONSTRUCT A DF out of the results 
    results_df=pd.DataFrame(distance_matrix)
    ## set the names of the dataframe 
    results_df.columns=col_names
    results_df.index=row_names
    ## return the results 
    return results_df

def compute_jaccard_index(exp1:Experiment,exp2:Experiment, level:str='peptide')->float:
    """Compute Jaccard index between samples two samples 

    Args:
        exp1 (Experiment): The first experimental instance 
        exp2 (Experiment): The first experimental instance 
        level (str): The level of computing the overlap between samples, can be any of peptide or protein 

    Returns:
        float: Jaccard index computed with regard to the to provide level
    """
    if level != 'peptide' and level != 'protein': 
        raise ValueError(f"Level: {level} is not supported, currently only level, peptide and protein are supported")
    if level=='peptide':
        return (len(exp1.get_peptides().intersection(exp2.get_peptides())) / len(exp1.get_peptides().union(exp2.get_peptides())))
    if level=='protein':
        return (len(exp1.get_proteins().intersection(exp2.get_proteins())) / len(exp1.get_proteins().union(exp2.get_proteins())))



