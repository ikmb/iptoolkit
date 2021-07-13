#!/usr/bin/env python 
"""The module contain functions that can be used for developing & testing other functions of the library 
"""
# import the modules 
import pandas as pd
from Bio import SeqIO
import numpy as np
from typing import List , Set, Dict 
import random
from IPTK.Classes.Experiment import Experiment
from IPTK.Classes.Proband import Proband
from IPTK.Classes.HLASet import HLASet 
from IPTK.Classes.Tissue import Tissue
from IPTK.Classes.Database import SeqDB
from IPTK.Utils.UtilityFunction import generate_random_name
# define the functions
def simulate_an_experimental_ident_table_from_fasta(path2load: str, num_pep: int, num_prot: int )->pd.DataFrame:
    """Simulate an IP identification table from a fasta file. Please Note,  if the reminder of num_pep over num_prot does not equal to zero,
    the floor of this ratio will be used to sample peptides from each proteins 

    :param path2load:  The path to load the Fasta files 
    :type path2load: str
    :param num_pep: The number of peptides in the tables 
    :type num_pep: int
    :param num_prot: The number of UNIQUE proteins in the table
    :type num_prot: int
    :raises ValueError: if number of proteins or number of peptide is zero 
    :return: an identification table 
    :rtype: pd.DataFrame
    """
    # load the fasta files
    # check that number of proteins & number of peptides does not equal zero 
    if not num_pep or not num_prot:
        raise ValueError('number of peptides or number of proteins can not be zero')
    fasta_gene=SeqIO.parse(path2load,'fasta')
    proteins_dict=dict()
    # load the desired number of proteins
    idx=0 
    while idx < num_prot:
        temp_seq=next(fasta_gene)
        if len(str(temp_seq.seq)) <50:
            idx-=1
            continue
        proteins_dict[temp_seq.id.split('|')[1]]=str(temp_seq.seq)
        idx+=1
    # get the average number of peptides from proteins 
    num_peptides_=int(np.floor(num_pep/num_prot))
    # generate the results 
    peptides,proteins,start_idx,end_idx=[],[],[],[]
    proteins_names=list(proteins_dict.keys())
    for prot_idx in range(num_prot):
        protine_seq=proteins_dict[proteins_names[prot_idx]]
        for idx in range(num_peptides_):
            pep_length=int(np.random.normal(15,3))
            anchor_point=np.random.randint(low=pep_length,high=len(protine_seq)-pep_length)
            peptides.append(protine_seq[anchor_point:anchor_point+pep_length])
            proteins.append(proteins_names[prot_idx])
            start_idx.append(anchor_point)
            end_idx.append(anchor_point+pep_length)
    # return the results as a data.frame 
    res=pd.DataFrame({
        'peptides':peptides,
        'proteins':proteins,
        'start_index':start_idx,
        'end_index':end_idx
    })
    return res

def simulate_an_expression_table(num_transcripts: int = 100) ->pd.DataFrame:
    """ create a dummy expression table to be used for testing and developing Tissue based classes 

    :param num_transcripts: The number of transcripts that shall be contained in the transcript, defaults to 100
    :type num_transcripts: int, optional
    :raises ValueError: incase number of transcripts is 0
    :return: a stimulated expression label 
    :rtype: pd.DataFrame
    """
    # check that the input is a positive integer 
    if num_transcripts<=0:
        raise ValueError(f'Number of Transcripts must be bigger than zero')
    # simulate the data 
    trans_id, prot_name, trans_value= [],[],[]
    for _ in range(num_transcripts):
        trans_id.append('ENS_'+generate_random_name(12))
        prot_name.append(generate_random_name(12))
        trans_value.append(np.random.normal(25,10))
    # create the data table 
    return pd.DataFrame(
        {
            'Transcript_id':trans_id,
            'Protein_name': prot_name,
            'Expression_value':trans_value
        }
    )

def simulate_random_experiment(alleles: List[str],  path2fasta: str,  tissue_name: str='TEST_TISSUE',
        num_pep: int = 10, num_prot: int = 5, proband_name: str = None )->Experiment:
    """ Simulate a random experiment objects 

    :param alleles: a list of alleles names. 
    :type alleles: List[str]
    :param path2fasta: The path to load the database objects 
    :type path2fasta: str
    :param tissue_name: the name of the tissue, defaults to 'TEST_TISSUE'
    :type tissue_name: str, optional
    :param num_pep: the number of peptides in the table, defaults to 10
    :type num_pep: int, optional
    :param num_prot: number of proteins, defaults to 5
    :type num_prot: int, optional
    :param proband_name: The name of the Proband, defaults to None
    :type proband_name: str, optional
    :return: A simulated experimental object 
    :rtype: Experiment
    """
    if proband_name is None:
        proband_name=generate_random_name(12)
    
    proband: Proband = Proband(name=proband_name)
    hla_set: HLASet = HLASet(alleles)
    ident_table: pd.DataFrame =simulate_an_experimental_ident_table_from_fasta(path2fasta,num_pep,num_pep)
    # to be upgraded to the new version of the Tissue class 
    tissue: Tissue = Tissue(tissue_name,simulate_an_expression_table(num_transcripts=1000),
            simulate_an_expression_table(num_transcripts=100))
    database: SeqDB = SeqDB(path2fasta)
    return Experiment(proband=proband,hla_set=hla_set,tissue=tissue,database=database,
    ident_table=ident_table)

def generate_random_peptide_seq(peptide_length:int, num_peptides: int)->List[str]:
    """generate a list of random peptides for testing and development purposes.
    
    :param peptide_length: The peptide length 
    :type peptide_length: int
    :param num_peptides: The number of peptides in the generate list  
    :type num_peptides: int
    :return: a list of random peptides 
    :rtype: List[str]
    """
    # define amino acids 
    amino_acids=['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
    results=[]
    # generate the sequences
    for _ in range(num_peptides):
        results.append(
            ''.join([random.choice(amino_acids) for _ in range(peptide_length)])
        )
    # return the results 
    return results

def simulate_mapped_array_list(min_len: int = 20, max_len: int = 100, num_elem: int = 100)->List[np.ndarray]:
    """Simulate a list of mapped arrays proteins to be used for development purposes 
    
    :param min_len: the minmum length of the protein, defaults to 20
    :type min_len: int, optional
    :param max_len: the maximum length for the protein, defaults to 100
    :type max_len: int, optional
    :param num_elem: the number of arrays in the protein, defaults to 100
    :type num_elem: int, optional
    :return: A list of simulated NumPy array that represent protein peptide coverage
    :rtype: List[np.ndarray]
    """
    # allocate a list to hold the results 
    results: List[np.ndarray] = []
    # fill the list with arrays: 
    for _ in range(num_elem):
        array_length: int = np.random.randint(low=min_len,high=max_len)
        results.append(np.random.randint(low=0,high=2,size=(1,array_length)))
    # return the results 
    return results