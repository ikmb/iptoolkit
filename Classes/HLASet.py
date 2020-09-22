#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@contact: h.elabd@ikmb.uni-kiel.de
"""
# load the modules 
from __future__ import annotations
import numpy as np 
from IPTK.Classes.HLAMolecules import HLAMolecule
from IPTK.Utils.Types import HLA_Names
from mhcnames import allele_name
from typing import List 
# define the class 
class HLASet: 
    def __init__(self, hlas: HLA_Names, gene_sep: str = ':')-> HLASet:
        """
        @brief: initialize an HLASet class which is a collection of HLA molecules 
        @param: hlas: a list of alleles to be added to the set  
        @param: gene_sep: Incase of HLA-DP and HLA-DQ, what is the gene separator that separate the genes' names, 
        for example, in DQA1*0303:DQB1*0202 the separator is the colon, Default is : 
        """
        # check that the provided list is not empty 
        if len(hlas)==0:
            raise ValueError(f"The provided list is empty!") 
        # add the alleles to the dict of objects 
        self._hlas=dict()
        for hla in hlas: 
            try:
                if 'DQ' in hla or 'DP' in hla:
                    chains=hla.split(gene_sep)
                    if len(chains)>2:
                        raise ValueError(f'splitting the allele name: {hla} with separator: {gene_sep} generated {len(chains)} chains!!, at max two chain are allowed')
                    elif len(chains)==2:
                        self._hlas[hla]=HLAMolecule(chain_one=chains[0], chain_two=chains[1])
                    else: 
                        self._hlas[hla]=HLAMolecule(chain=hla[0])
                else: 
                    self._hlas[hla]=HLAMolecule(chain=hla)
            except Exception as exp: 
                raise RuntimeError(f'While adding the allele: {hla} to the set, the following error : {exp} was encountered.')
        # check that all alleles belong to the same class 
        self._assert_same_class()
        return 

    def get_hla_count(self)->int:
        """
        @brief: get the count of HLA molecules in the set 
        """
        return len(self)
    
    def _assert_same_class(self)->None:
        """
        @brief: assert that all the alleles in the set are of the same class, i.e. HLA-I or HLA-II
        """
        # get the HLA-class 
        # get a radom class from the input HLA 
        default_class=self._hlas[list(self._hlas.keys())[0]].get_class() # get the class of a random molecules 
        # check that all the alleles belong to the same class 
        for allele_name in self._hlas.keys(): 
            if self._hlas[allele_name].get_class()!=default_class: 
                raise ValueError(f"ASSERTION FAILED: the alleles in the set belong to a different classes")
        
    def get_class(self)->int:
        """
        @brief: return the class of the HLA-alleles in the current instance 
        """
        return self._hlas[list(self._hlas.keys())[0]].get_class()
    
    def get_alleles(self)-> List[str]:
        """
        @brief: returns the alleles of the set 
        """
        return list(self._hlas.keys())
    
    def has_allele(self,allele: str)->bool:
        """
        @brief: check whether or not the provided allele existed in the list of alleles
        @param: allele: the name of the allele to check the database for 
        """
        allele_names=list(self._hlas.keys())
        for name in allele_names:
            if allele == name or allele == name.strip('HLA').strip('-') or allele.strip('HLA').strip('-') == name :
                return True
        return False
    
    def has_gene(self,gene_name:str)->bool:
        """
        @brief check if at least one of the alleles in the set belongs to the provided locus
        @param gene_name:  the gene name to search the set against.
        """
        for allele_name in self._hlas.keys(): 
            if gene_name in self._hlas[allele_name].get_gene(): 
                return True 
        return False 
    
    def has_allele_group(self, allele_group:str)->bool:
        """
        @brief check if at least one allele in the set belongs to the provided allele group
        @param allele_name: the allele group to search the set for 
        """
        for allele_name in self._hlas.keys(): 
            if allele_group in self._hlas[allele_name].get_allele_group():
                return True 
        return False 
    
    def has_protein_group(self,protein_group:str)->bool:
        """
        @brief check if at least one allele in the set belongs to the provided protein group
        @param protein_group: the protein group to search the set for
        """
        for allele_name in self._hlas.keys(): 
            if protein_group in self._hlas[allele_name].get_protein_group():
                return True 
        return False 

    # define some magic functions
    def __len__(self)->int:
        """
        @brief: a magic len function that return the number of alleles in the set 
        """
        return len(self._hlas)

    def __str__(self) ->str:
        """
        @brief: a magic function to compute a string form of the class  
        """
        return f"An HLASet containing {len(self)} alleles"
    
    def __repr__(self)->str:
        """
        @brief: a magic function to compute the class representation 
        """
        return str(self)
            
            
            
            
            