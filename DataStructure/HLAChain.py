#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: The implementation of an HLA molecules  
@version: 0.0.1
@date: 19.08.2020  
"""
# load the modules
from __future__ import annotations 
import mhcnames
# define the class 
class HLAChain: 
    def __init__(self, name:str)->HLAChain:
        """
        @breif: create an instance an HLA form a name string 
        @param name: the name of the HLA instance 
        """
        self._name=name
        name=mhcnames.parse_allele_name(name)
        # extract fields from the name
        self._hla_class=self.get_chain_class(name.gene)
        self._gene=name.gene
        self._allele_group=name.allele_family
        self._protein_group=name.allele_code

    def get_class(self)->int:
        """
        @brief: get the HLA class 
        """
        return self._hla_class 
    
    def get_gene(self)->str:
        """
        @brief get the gene name
        """
        return self._gene
    
    def get_allele_group(self)->str:
        """
        @brief get the allele group 
        """
        return self._allele_group
    
    def get_protein_group(self)->str:
        """
        @brief get the protein name
        """
        return self._protein_group
    
    def get_chain_class(self, gene_name:str)->int:
        """
        @brief check the allele name and return 1 if the gene belongs to class one 
        and 2 if it belong to class two 
        @param gene_name the name of the gene 
        """
        if gene_name in ['A','B','C']: 
            return 1
        return 2
         
    def get_name(self)->str:
        """
        @brief: return the chain name
        """
        return self._name

    # define some magic number for the class 
    def __str__(self)->str:
        """
        @brief return a string representation of the class 
        """
        return f"""An HLA chain of class: {self.get_class()} from gene: {self.get_gene()},
                With an allele group of: {self.get_allele_group()} and a protein group of: 
                {self.get_protein_group()}   
                """

    def __repr__(self)->str:
        """
        @brief compute the class string representation
        """
        return str(self)
        
