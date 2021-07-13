#!/usr/bin/env python 
"""a representation of an HLA molecules 
"""
# load the modules 
from __future__ import annotations
from IPTK.Classes.HLAChain import HLAChain
from IPTK.Utils.Types import HLA_Chains,Genes,AlleleGroup,ProteinGroup
# define the class 
class HLAMolecule:
    def __init__(self,**hla_chains)->HLAMolecule:
        """ Create an HLA molecule which is represented as a molecule that is composite of at max two HLA chains
        :raises ValueError: If number of chains is bigger than 2
        :raises ValueError: If the number of chains is 0
        :raises RuntimeError: Captsure any exception than might be encountered while creating the chains.
        :raises ValueError: If the provided chain belong to different classes, for example class one and class two
        :return: an HLAMolecule instance 
        :rtype: HLAMolecule
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
        :param sep: The name of the allele by concatenating the names of the individual chains using \
        a separator, defaults to ':'
        :type sep: str, optional
        :return: [description]
        :rtype: str
        """
        names=list(self._chains.keys())
        names.sort() # sort the list to ensure Alpha->Beta order 
        if len(names)==2:
            return names[0].replace(':','')+sep+names[1].replace(':','')
        return names[0]
    
    def get_class(self)->int: 
        """
        :return: The class of the HLA molecules 
        :rtype: int
        """
        return self._chains[list(self._chains.keys())[0]].get_class()
    
    def get_gene(self)->Genes:
        """
        :return: return gene/pair of genes coding for the current HLA molecules 
        :rtype: Genes
        """
        if len(self._chains)==2:
            names=list(self._chains.keys())
            return [self._chains[names[0]].get_gene(), self._chains[names[1]].get_gene()]
        return [self._chains[list(self._chains.keys())[0]].get_gene()]
    
    def get_allele_group(self)->AlleleGroup:
        """
        :return: The allele group for the instance chain/pair of chains 
        :rtype: AlleleGroup
        """
        if len(self._chains)==2:
            names=list(self._chains.keys())
            return [self._chains[names[0]].get_allele_group(), self._chains[names[1]].get_allele_group()]
        return [self._chains[list(self._chains.keys())[0]].get_allele_group()]
    
    def get_protein_group(self)->ProteinGroup:
        """
        :return: The protein group for the instance chain/pair of chains 
        :rtype: ProteinGroup
        """
        if len(self._chains)==2:
            names=list(self._chains.keys())
            return [self._chains[names[0]].get_protein_group(), self._chains[names[1]].get_protein_group()]
        return [self._chains[list(self._chains.keys())[0]].get_protein_group()]



