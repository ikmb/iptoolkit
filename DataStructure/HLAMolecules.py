#!/usr/bin/env python 
"""
@brief: a representation of an HLA molecules 
@author: Hesham ElAbd
@version: 0.0.1
"""
# load the modules 
from __future__ import annotations
from IPTK.DataStructure.HLAChain import HLAChain
from IPTK.Utils.Types import HLA_Chains,Genes,AlleleGroup,ProteinGroup
# define the class 
class HLAMolecule:
    def __init__(self,**hla_chains)->HLAMolecule:
        """
        @brief: an HLA molecule which is represent as a molecules that is composite of at max two HLA chains
        @param: at max two key value pairs where the key is just a place holder and the value is the chain name, see the example below
        @example: 
        >>> test_molecule=HLAMolecule(chain_one=DQA1*03:03,chain_two=DQB1:0202)
        >>> test_molecule.get_names() # get the name 
        DQA1*0303:DQB1*0202
        """
        # check that the dict contain at max two chains 
        if len(hla_chains) > 2:
            raise ValueError(f"The provided HLA-Chains MUST contain at max 2 chains, however, your input contains {len(hla_chains)}")
        if len(hla_chains)==0:
            raise ValueError('To construct an HLA molecules, you need at least one hla chains, however, the number of provided chain is 0')
        # create a dict to store the results
        self._chains=dict() 
        # check that if the dict has more than one chain that they belong to the same gene 
        names=list(hla_chains.values())
        number_of_chains=len(names)
        while number_of_chains+1:
            # parse the first chain 
            try:
                self._chains[names[number_of_chains-1]]=HLAChain(names[number_of_chains-1])
            except Exception as exp: 
                raise RuntimeError(f'while parsing the chain named: {names[number_of_chains-1]} The following error was encountered: {exp}')
            number_of_chains-=1
        # Check the chains have the same class incase more than one chain have been provided 
        if len(names)==2:
            if self._chains[names[0]].get_class() != 2 or  self._chains[names[0]].get_class() !=2: 
                raise ValueError(f'Incase more than one chain is provided, it must be from an HLA-II molecules')
    
    def get_name(self,sep:str = ':')->str:
        """
        @brief: return the name of the allele by concatenating the names of the individual chains using 
        a separator 
        """
        names=list(self._chains.keys())
        names.sort() # sort the list to ensure Alpha->Beta order 
        if len(names)==2:
            return names[0].replace(':','')+sep+names[1].replace(':','')
        return names[0]
    
    def get_class(self)->int: 
        """
        @brief: return the class of the HLA molecules 
        """
        return self._chains[list(self._chains.keys())[0]].get_class()
    
    def get_gene(self)->Genes:
        """
        @brief: return gene/pair of genes coding for the current HLA molecules 
        """
        if len(self._chains)==2:
            names=list(self._chains.keys())
            return [self._chains[names[0]].get_gene(), self._chains[names[1]].get_gene()]
        return [self._chains[list(self._chains.keys())[0]].get_gene()]
    
    def get_allele_group(self)->AlleleGroup:
        """
        @brief: return the allele group for the instance chain/pair of chains 
        """
        if len(self._chains)==2:
            names=list(self._chains.keys())
            return [self._chains[names[0]].get_allele_group(), self._chains[names[1]].get_allele_group()]
        return [self._chains[list(self._chains.keys())[0]].get_allele_group()]
    
    def get_protein_group(self)->ProteinGroup:
        """
        @brief: return the protein group for the instance chain/pair of chains 
        """
        if len(self._chains)==2:
            names=list(self._chains.keys())
            return [self._chains[names[0]].get_protein_group(), self._chains[names[1]].get_protein_group()]
        return [self._chains[list(self._chains.keys())[0]].get_protein_group()]



