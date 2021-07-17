"""Impelement Numba JIT accelerated functions for scaling up IPTK code
"""
## load the modules 
from numba.np.ufunc import parallel
import numpy as np
from numba import jit
from typing import List, Dict, Set, Union

@jit(nopython=True,nogil=True,parallel=True,cache=True)
def _get_mapped_proteins(start_idx:np.ndarray, end_idx:np.ndarray, protin_len:int)->np.ndarray:
    """Compute mapped protein array 

    Args:
        start_idx:(np.ndarray): the start boundires of child peptide 
        end_idx (np.ndarray): the end boundires of child peptide 
        protin_len (int): the length of the parent protein

    Returns:
        np.ndarray: an np.ndarray of shape(1,protein_len) containing the mapping of peptides to protein
    """
    protein_backbone=np.zeros(shape=(1,protin_len),dtype=np.int32)
    for i,j in zip(start_idx,end_idx): protein_backbone[0,i:j]+=1
    return protein_backbone

@jit(nopython=True,nogil=True,parallel=False,cache=True)
def _get_non_presented_peptides(exc_reg_s_idx: int, exc_reg_e_idx: int, length: int, seq:str)->str:
    """ get a non-presented peptides from the protein backbone 

    Args:
        exc_reg_s_idx (int): the start position of the represented peptide 
        exc_reg_e_idx (int): the end position of the represented peptide
        length (int): the length of the non-represented peptide 
        seq (str): the protin sequence 

    Returns:
        str: a non-presented peptide which is a subsequence of the protein sequence 
    """
    try_flag=True
    while try_flag: 
        anc_point=np.random.randint(low=0,high=len(seq)-length)
        if anc_point < exc_reg_s_idx or anc_point > exc_reg_e_idx: 
            neg_seq=seq[anc_point:anc_point+length] 
            try_flag=False
    return neg_seq