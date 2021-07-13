#!/usr/bin/env Python 
""" An abstraction for a collection of HLA alleles
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
        """create an HLASet class which is a collection of HLA molecules 
        :param hlas: a list of alleles to be added to the set  
        :type hlas: HLA_Names
        :param gene_sep: Incase of HLA-DP and HLA-DQ, what is the gene separator that separate the genes' names, 
        for example, in DQA1*0303:DQB1*0202 the separator is the colon, defaults to ':'.
        :type gene_sep: str, optional
        :raises ValueError: incase the list of alleles is empty 
        :raises ValueError: incase of HLA-DP and HLA-DQ with a mismatch allele names.
        :raises RuntimeError: incase adding the alleles failed for any reason.
        :return: an HLASet instance 
        :rtype: HLASet
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
        :return: The count of HLA molecules in the set 
        :rtype: int
        """
        return len(self)
    
    def _assert_same_class(self)->None:
        """
        :raises ValueError: if the alleles in the set are of different classes, for example,  HLA-I or HLA-II
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
        :return: The class of the HLA-alleles in the current instance 
        :rtype: int
        """
        return self._hlas[list(self._hlas.keys())[0]].get_class()
    
    def get_alleles(self)-> List[str]:
        """
        :return: The current alleles in the set
        :rtype: int
        """
        return list(self._hlas.keys())
    
    def has_allele(self,allele: str)->bool:
        """
        :param allele: The name of the alleles to check for its occurrence in the instance. 
        :type allele: str
        :return: True, if the provided allele is in the current instance, False otherwise. 
        :rtype: bool
        """
        allele_names=list(self._hlas.keys())
        for name in allele_names:
            if allele == name or allele == name.strip('HLA').strip('-') or allele.strip('HLA').strip('-') == name :
                return True
        return False
    
    def has_gene(self,gene_name:str)->bool:
        """
        :param gene_name: the gene name to search the set against.
        :type gene_name: str
        :return: True, if at least one of the alleles in the set belongs to the provided gene. False otherwise
        :rtype: bool
        """
        for allele_name in self._hlas.keys(): 
            if gene_name in self._hlas[allele_name].get_gene(): 
                return True 
        return False 
    
    def has_allele_group(self, allele_group:str)->bool:
        """
        :param allele_group: The allele group to search the set for 
        :type allele_group: str
        :return: True, if at least one allele in the set belongs to the provided allele group, False otherwise. 
        :rtype: bool
        """
        for allele_name in self._hlas.keys(): 
            if allele_group in self._hlas[allele_name].get_allele_group():
                return True 
        return False 
    
    def has_protein_group(self,protein_group:str)->bool:
        """
        :param protein_group: The protein group to search the set for
        :type protein_group: 
        :return: True, if at least one allele in the set belongs to the provided protein group
        :rtype: bool
        """
        for allele_name in self._hlas.keys(): 
            if protein_group in self._hlas[allele_name].get_protein_group():
                return True 
        return False 

    # define some magic functions
    def __len__(self)->int:
        """
        :return: The number of alleles in the set 
        :rtype: int
        """
        return len(self._hlas)

    def __str__(self) ->str:
        """
        :return: A string form of the class  
        :rtype: str
        """
        return f"An HLASet containing {len(self)} alleles."
    
    def __repr__(self)->str:
        return str(self)
    
    def get_names(self)->List[str]: 
        """Return a list of all HLA allele names defined in the set

        Returns:
            List[str]: [description]
        """
        return [molecule.get_names() for molecule in self._hlas.values()]
            
            
            
            