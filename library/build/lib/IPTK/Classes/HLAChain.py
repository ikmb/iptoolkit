#!/usr/bin/env Python 
"""The implementation of an HLA molecules  
"""
# load the modules
from __future__ import annotations 
import mhcnames
# define the class 
class HLAChain: 
    def __init__(self, name:str)->HLAChain:
        """Create an HLAChain instance. 
        :param name: the allele name 
        :type name: str
        :return: an HLAChain instance 
        :rtype: HLAChain
        """
        self._name:str=name
        name=mhcnames.parse_allele_name(name)
        # extract fields from the name
        self._hla_class=self.get_chain_class(name.gene)
        self._gene=name.gene
        self._allele_group=name.allele_family
        self._protein_group=name.allele_code

    def get_class(self)->int:
        """
        :returns: The HLA class 
        :rtype: int
        """
        return self._hla_class 
    
    def get_gene(self)->str:
        """
        :return: The gene name
        :rtype: str
        """
        return self._gene
    
    def get_allele_group(self)->str:
        """
        :return: The allele group 
        :rtype: str
        """
        return self._allele_group
    
    def get_protein_group(self)->str:
        """
        :return: The protein name
        :rtype: str
        """
        return self._protein_group
    
    def get_chain_class(self, gene_name:str)->int:
        """
        :param gene_name: the name of the gene 
        :type gene_name: str
        :return: 1 if the gene belongs to class one and 2 if it belong to class two 
        :rtype: int
        """
        if gene_name in ['A','B','C']: 
            return 1
        return 2
         
    def get_name(self)->str:
        """
        :return: The chain name
        :rtype: str
        """
        return self._name

    # define some magic number for the class 
    def __str__(self)->str:
        """
        :return: a string representation of the class
        :rtype: str
        """
        return f"""An HLA chain of class: {self.get_class()} from gene: {self.get_gene()},
                With an allele group of: {self.get_allele_group()} and a protein group of: 
                {self.get_protein_group()}.   
                """

    def __repr__(self)->str:
        return str(self)
        
