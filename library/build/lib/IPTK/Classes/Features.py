#!/usr/bin/env python 
"""Parses the XML scheme of a uniprot protein and provides a python API  
for quering and accessing the results 
"""
# load the modules 
from __future__ import annotations
import numpy as np 
from Bio import SeqIO
import urllib
import os 
from typing import List, Tuple, Dict, Union
class Features:
    r"""The class provides a template for the features associated with a protein. \
    The following features are associated with the protein \ 
    #signal peptide: dict \
        The range of the signal peptides, if the protein has no signal, for example, a globular \
        cytologic protein. None is used as a default, placeholder \
        value. \ 
    #chains:dict \ 
        the chains making up the mature protein, the protein should at least have one chain. \ 
    #domain: dict\ 
        the known domains in the protein, if no domain is defined, None is used. \ 
    #modification sites: nested dict \ 
        that contains information about the PTM sites, glycosylation site and disulfide bonds. \ 
    #sequence variances: dict \ 
        which contains information about the sequence variants of a protein structure. \ 
    #split variance: dict \ 
        which contain known splice variants \ 
    ** Notes: Although disulfide bond is not a PTMs, it is being treated as a \ 
    one here to simplify the workflow. \ 
    """
    def __init__(self,uniprot_id:str, temp_dir: str = None)->Features:
        """Download the features associated with a protein from uniprot and then parse the results using SeqIO to extract the features
        
        :param uniprot_id: Uniprot id to download its XML scheme from Uniprot 
        :type uniprot_id: str 
        :param temp_dir: a temporary directory to download the XML scheme to it, if not provided files are download to the current working directory, defaults to None.
        :type temp_dir: str
        """
        # download the files 
        if temp_dir is None: 
            temp_dir="." 
        # define the file path 
        file_path: str = os.path.join(temp_dir,f"{uniprot_id}.xml")
        try: 
            urllib.request.urlretrieve(f"https://www.uniprot.org/uniprot/{uniprot_id}.xml", file_path) 
        except urllib.error.HTTPError: 
            raise IOError("Downloading the provided fail, Check your provided uniprot id")
        # read the sequence object 
        record: SeqIO.SeqRecord = SeqIO.read(file_path,"uniprot-xml")
        # fill the provided record 
        self.extracted_features=dict()
        # parse the features in the record.feature object
        for feature in record.features:
            # extract sequences signal peptide
            if feature.qualifiers["type"]=="signal peptide":
                # Try to extract the start and end position for the chain
                # if this could no be extracted None is used as a default value
                try:
                    start_index=int(feature.location.start)
                except(TypeError):
                    start_index=None
                try:
                    end_index=int(feature.location.end)
                except(TypeError):
                    end_index=None
                self.extracted_features["SignalPeptide"]={
                    "startIdx":start_index,
                    "endIdx":end_index
                    }
            # extract the chain information
            elif feature.qualifiers["type"]=="chain":
                if "Chains" in self.extracted_features.keys():
                    chainIdx=len(self.extracted_features["Chains"]) # get the chain index
                    chainName="chain_number_"+str(chainIdx)
                    # Try to extract the start and end position for the chain
                    # if this colud no be extracted None is used as a default value
                    try:
                        start_index=int(feature.location.start)
                    except(TypeError):
                        start_index=None
                    
                    try:
                        end_index=int(feature.location.end)
                    except(TypeError):
                        end_index=None
                    self.extracted_features["Chains"][chainName]={
                        "chainId":feature.id,
                        "startIdx":start_index,
                        "endIdx":end_index
                        }
                else:
                    
                    try:
                        start_index=int(feature.location.start)
                    except(TypeError):
                        start_index=None
                    
                    try:
                        end_index=int(feature.location.end)
                    except(TypeError):
                        end_index=None
                        
                    self.extracted_features["Chains"]={
                        "chain_number_0":{
                        "chainId":feature.id,
                        "startIdx":start_index,
                        "endIdx":end_index
                            }}
            # extract the domain information
            elif feature.qualifiers["type"]=="domain":
                if "Domains" in self.extracted_features.keys():
                    self.extracted_features["Domains"][
                    feature.qualifiers["description"]]={
                        "startIdx":int(feature.location.start),
                        "endIdx":int(feature.location.end)}
                else: 
                    self.extracted_features["Domains"]={
                        feature.qualifiers["description"]:{
                        "startIdx":int(feature.location.start),
                        "endIdx":int(feature.location.end)
                            }
                        }
            # extract transmembrane features:
            elif feature.qualifiers["type"]=='transmembrane region':
                if "transmembrane_region" in self.extracted_features.keys():
                    self.extracted_features["transmembrane_region"].append((int(feature.location.start),int(feature.location.end)))
                else:
                    self.extracted_features["transmembrane_region"]=[(int(feature.location.start),int(feature.location.end))]
            # extract the topological information:
            elif feature.qualifiers["type"]=="modified residue":
                if "PTMs" in self.extracted_features.keys():
                    if "Modifications" in self.extracted_features["PTMs"].keys():
                        modificationIdx=len(self.extracted_features["PTMs"][
                            "Modifications"])
                        modificationName="SeqModification_num_"+str(modificationIdx)
                        self.extracted_features["PTMs"]["Modifications"][
                            modificationName]={
                                "Name":feature.qualifiers["description"],
                                "startIdx":int(feature.location.start),
                                "endIdx":int(feature.location.end)
                                }
                    else:
                        self.extracted_features["PTMs"]["Modifications"]={}
                        self.extracted_features["PTMs"]["Modifications"][
                            "SeqModification_num_0"]={
                            "Name":feature.qualifiers["description"],
                            "startIdx":int(feature.location.start),
                            "endIdx":int(feature.location.end)
                            }
                else: 
                    self.extracted_features["PTMs"]={}
                    self.extracted_features["PTMs"]["Modifications"]={}
                    self.extracted_features["PTMs"]["Modifications"][
                        "SeqModification_num_0"]={
                            "Name":feature.qualifiers["description"],
                            "startIdx":int(feature.location.start),
                            "endIdx":int(feature.location.end)
                            }
            # extract and add glycosylation sites to the features class
            elif feature.qualifiers["type"]=="glycosylation site":
                if "PTMs" in self.extracted_features.keys():
                    if "GlycoSite" in self.extracted_features["PTMs"].keys():
                        glycositeIdx=len(self.extracted_features["PTMs"]["GlycoSite"])
                        glyco_site_name="Glyco_Site_number_"+str(glycositeIdx)
                        self.extracted_features["PTMs"]["GlycoSite"][
                            glyco_site_name]={
                                "Name":feature.qualifiers["description"],
                                "startIdx":int(feature.location.start),
                                "endIdx":int(feature.location.end)
                                }
                    else:
                        self.extracted_features["PTMs"]["GlycoSite"]={}
                        self.extracted_features["PTMs"]["GlycoSite"][
                            "Glyco_Site_number_0"]={
                                "Name":feature.qualifiers["description"],
                                "startIdx":int(feature.location.start),
                                "endIdx":int(feature.location.end)
                                }
                else:
                    self.extracted_features["PTMs"]={}
                    self.extracted_features["PTMs"]["GlycoSite"]={}
                    self.extracted_features["PTMs"]["GlycoSite"][
                            "Glyco_Site_number_0"]={
                                "Name":feature.qualifiers["description"],
                                "startIdx":int(feature.location.start),
                                "endIdx":int(feature.location.end)
                                }
            # extract and add disulfide site to the feature class
            elif feature.qualifiers["type"]=="disulfide bond":
                if "PTMs" in self.extracted_features.keys():
                    if "DisulfideBond" in self.extracted_features["PTMs"].keys():
                        disulfideBondIdx=len(self.extracted_features["PTMs"]["DisulfideBond"])
                        disulfide_site_name="disulfide_site_number_"+str(disulfideBondIdx)
                        # try to extract the start and the end position.
                        try:
                            start_index=int(feature.location.start)
                        except(TypeError):
                            start_index=None
                        try:
                            end_index=int(feature.location.end)
                        except(TypeError):
                            end_index=None
                        self.extracted_features["PTMs"]["DisulfideBond"][
                            disulfide_site_name]={
                                "Name":feature.qualifiers["type"],
                                "startIdx":start_index,
                                "endIdx":end_index
                                }
                    else:
                        try:
                            start_index=int(feature.location.start)
                        except(TypeError):
                            start_index=None
                        try:
                            end_index=int(feature.location.end)
                        except(TypeError):
                            end_index=None
                        self.extracted_features["PTMs"]["DisulfideBond"]={}
                        self.extracted_features["PTMs"]["DisulfideBond"][
                            "disulfide_site_number_0"]={
                                "Name":feature.qualifiers["type"],
                                "startIdx":start_index,
                                "endIdx":end_index
                                }
                else:
                    self.extracted_features["PTMs"]={}
                    self.extracted_features["PTMs"]["DisulfideBond"]={}
                    self.extracted_features["PTMs"]["DisulfideBond"][
                            "disulfide_site_number_0"]={
                                "Name":feature.qualifiers["type"],
                                "startIdx":int(feature.location.start),
                                "endIdx":int(feature.location.end)
                                }
            # extract sequence variant from the protein sequence:
            elif feature.qualifiers["type"]=="sequence variant":
                if "SeqVar" in self.extracted_features.keys():
                    varient_index=len( self.extracted_features["SeqVar"])
                    varient_name="Sequence_varient_number_"+str(varient_index)
                    # check if the feature has a description entrey
                    if "description" in feature.qualifiers.keys():
                        snp_id=feature.qualifiers["description"]
                    else:
                        snp_id=None
                    # check if the feature has a original entrey
                    if "original" in feature.qualifiers.keys():
                        original_amino_acid=feature.qualifiers["original"]
                    else:
                        original_amino_acid=None
                    # check if the feature has a varient entrey
                    if "variation" in feature.qualifiers.keys():
                        varient_amino_acid=feature.qualifiers["variation"]
                    else:
                        varient_amino_acid=None
                    # fill the entries
                    self.extracted_features["SeqVar"][varient_name]={
                        "VarientId":feature.qualifiers["id"],
                        "SNP_Id":snp_id,
                        "original":original_amino_acid,
                        "variation":varient_amino_acid,
                        "startIdx":int(feature.location.start),
                        "endIdx":int(feature.location.end)
                        }
                else:
                    self.extracted_features["SeqVar"]={}
                    if "description" in feature.qualifiers.keys():
                        snp_id=feature.qualifiers["description"]
                    else:
                        snp_id=None
                    # check if the feature has a description entrey
                    if "description" in feature.qualifiers.keys():
                        snp_id=feature.qualifiers["description"]
                    else:
                        snp_id=None
                    # check if the feature has a original entrey
                    if "original" in feature.qualifiers.keys():
                        original_amino_acid=feature.qualifiers["original"]
                    else:
                        original_amino_acid=None
                    # check if the feature has a varient entrey
                    if "variation" in feature.qualifiers.keys():
                        varient_amino_acid=feature.qualifiers["variation"]
                    else:
                        varient_amino_acid=None
                    self.extracted_features["SeqVar"][
                        "Sequence_varient_number_0"]={
                        "VarientId":feature.qualifiers["id"],
                        "SNP_Id":snp_id,
                        "original":original_amino_acid,
                        "variation":varient_amino_acid,
                        "startIdx":int(feature.location.start),
                        "endIdx":int(feature.location.end)
                        }
            # extract splice vaients from the protein sequences: 
            elif feature.qualifiers["type"]=="splice variant":
                if "SpliceVar" in self.extracted_features.keys():
                    SpliceVarIdx=len(self.extracted_features["SpliceVar"])
                    spliceVarient_name="splice_varient_number_"+str(SpliceVarIdx)
                    self.extracted_features["SpliceVar"][spliceVarient_name]={
                        "Name":feature.qualifiers["id"],
                        "Isoform":feature.qualifiers["description"],
                        "startIdx":int(feature.location.start),
                        "endIdx":int(feature.location.end)
                        }
                else: 
                    self.extracted_features["SpliceVar"]={}
                    self.extracted_features["SpliceVar"][
                        "splice_varient_number_0"]={
                        "Name":feature.qualifiers["id"],
                        "Isoform":feature.qualifiers["description"],
                        "startIdx":int(feature.location.start),
                        "endIdx":int(feature.location.end)
                        }
        # fill in the empty object with None:
        if "SignalPeptide" not in self.extracted_features.keys():
            self.extracted_features["SignalPeptide"]=None
        if "Chains" not in self.extracted_features.keys():
            self.extracted_features["Chains"]=None
        if "Domains" not in self.extracted_features.keys():
            self.extracted_features["Domains"]=None
        if "PTMs" not in self.extracted_features.keys():
            self.extracted_features["PTMs"]=None
        else:
            if "Modifications" not in self.extracted_features["PTMs"].keys():
                self.extracted_features["PTMs"]["Modifications"]=None 
            if "GlycoSite" not in self.extracted_features["PTMs"].keys():
                self.extracted_features["PTMs"]["GlycoSite"]=None 
            if "DisulfideBond" not in self.extracted_features["PTMs"].keys():
                self.extracted_features["PTMs"]["DisulfideBond"]=None 
        if "SeqVar" not in self.extracted_features.keys():
            self.extracted_features["SeqVar"]=None
        if "SpliceVar" not in self.extracted_features.keys():
            self.extracted_features["SpliceVar"]=None
        if 'transmembrane_region' not in self.extracted_features.keys():
            self.extracted_features["transmembrane_region"]=None
        return 
    # accessor methods:
    def has_signal_peptide(self)->bool:
        """
        :return: True if the protein has a signal peptide and False other wise.
        :rtype: bool
        """
        if self.extracted_features["SignalPeptide"]==None:
            return False
        return True
        
    def get_signal_peptide_index(self)->Tuple[int,int]:
        """
        :return:  The Index of the signal-peptide in the protein, if not signal peptide is defined, it returns None
        :rtype: Tuple[int,int]
        """
        if self.extracted_features["SignalPeptide"]==None:
            return None,None
        startIdx=self.extracted_features["SignalPeptide"]["startIdx"]
        endIdx=self.extracted_features["SignalPeptide"]["endIdx"]
        return startIdx,endIdx
    
    def has_chains(self)->bool:
        """
        :return: True if the protein has/have chain/chains as feature and False otherwise.
        :rtype: [type]
        """
        if self.extracted_features["Chains"]==None:
            return False
        return True
    
    def get_number_chains(self) -> int:
        """
        :return: The number of chains in the protein. if no chain is defined it returns zero.
        :rtype: int
        """
        if not self.has_chains():
            return 0
        return len(self.extracted_features["Chains"])
    
    def get_chains(self)->Dict[Dict[str,Union[str,int]]]:
        """
        :return: A dictionary that contains the chains of the protein, if no chain is defined it return None
        :rtype: Dict[Dict[str,Union[str,int]]]
        """
        return self.extracted_features["Chains"]
    
    def has_transmembrane_domains(self)->bool:
        """
        Returns:
            bool: True if the protein has transmembrane region and false otherwise
        """
        return self.extracted_features["transmembrane_region"]!=None

    def get_transmembrane_regions(self)->List[Tuple[int,int]]:
        """return a list containing the boundaries of transmembrane regions in the protein 

        Returns:
            List[Tuple[int,int]]: a list containing the boundaries of transmembrane regions in the protein 
        """
        return self.extracted_features["transmembrane_region"]
    
    def get_num_transmembrane_regions(self)->int:
        """Return the number of transmembrane regions on the protein 

        Returns:
            int: Return the number of transmembrane regions on the protein
        """
        print()
        if self.extracted_features["transmembrane_region"]!=None:
            return len(self.get_transmembrane_regions())
        return 0

    def has_domains(self)->bool: 
        """
        :return: True if the protein has a defined domain/domains, otherwise it return False.
        :rtype: bool 
        """
        return  self.extracted_features["Domains"]!=None
    
    def get_number_domains(self)->int:
        """
        :return: The number of domains a protein has, if no domain is defined it returns zero.
        :rtype: int
        """
        if self.extracted_features["Domains"] ==None:
            return 0
        return len(self.extracted_features["Domains"])
    
    def get_domains(self)->Dict[str, Dict[str, int]]:
        """
        :return:  The domains defined in the protein sequence, if no domain is defined it returns None.
        :rtype: Dict[str, Dict[str, int]]
        """
        return self.extracted_features["Domains"]
    
    def has_PTMs(self)->bool:
        """
        :return:True if the protein has a PTMs and False other wise
        :rtype: bool
        """
        if self.extracted_features["PTMs"] ==None:
            return False
        return True

    def has_disulfide_bond(self)->bool:
        """
        :return: True is the protein has disulfide and False other wise
        :rtype: bool
        """
        if self.has_PTMs():
            if self.extracted_features["PTMs"]["DisulfideBond"]==None:
                return False
            else: 
                return True
        return False
    
    def has_glycosylation_site(self)->bool:
        """
        :return: True if the protein has a glycosylation site and False otherwise.
        :rtype: [type]
        """
        if self.has_PTMs():
            if self.extracted_features["PTMs"]["GlycoSite"] == None:
                return False
            else:
                return True
        return False
    
    def has_site_modifications(self)->bool:
        """
        :return: True if the protein has a modification site and False otherwise
        :rtype: bool
        """
        if self.has_PTMs():
            if self.extracted_features["PTMs"]["Modifications"] == None:
                return False
            else:
                return True
        return False
    
    def get_PTMs(self)->Dict[str,Dict[str,Dict[str,Union[str,int]]]]:
        """
        :return:  a nested dictionary that contains the PTMs found within the protein \
        the PTMs are classified into three main categories:

            1- Modifications: which is the generic case and contain information \
            about any sequence modification beside disulfide bonds and glycosylation.
            
            2- glycosylation: contains information about glycosylation sites
            
            3- DisulfideBond: contains information about disulfide bond

        :rtype: Dict[str,Dict[str,Dict[str,Union[str,int]]]]
        """
        return self.extracted_features["PTMs"]
    
    def get_PTMs_modifications(self)->Dict[str,Dict[str,Union[str,int]]]:
        """
        :return:  The generic modifications found on the protein. If the protein has no PTM, the function returns None.
        :rtype: Dict[str,Dict[str,Union[str,int]]]
        """
        if self.extracted_features["PTMs"] is None: return None
        if "Modifications" in self.extracted_features["PTMs"].keys():
            return self.extracted_features["PTMs"]["Modifications"]
        return None
    
    def get_PTMs_glycosylation(self)->Dict[str,Dict[str,Union[str,int]]]:
        """
        :return: The glycosylation sites found on the protein. If the protein has no glycosylation sites, the function returns None.
        :rtype: [type]
        """
        if self.extracted_features["PTMs"] is None: return None 
        return self.extracted_features["PTMs"]["GlycoSite"]
    
    def get_disulfide_bonds(self)->Dict[str,Dict[str,Union[str,int]]]:
        """
        :return: The disulfide sites found on the protein. If the protein has no disulfide sites, the function returns None
        :rtype: [type]
        """
        if self.extracted_features["PTMs"] is None: return None
        return self.extracted_features["PTMs"]["DisulfideBond"]
    
    def get_number_PTMs(self)->int:
        """
        :return:  The number of PTMs the sequence has, this include di-sulfide bonds. See Note1 for more details. \
        If the protein has no PTMs the function returns zero
        :rtype: int
        """
        if not self.has_PTMs():
            return 0
        else:
             number_of_PTMs=self.get_number_modifications()+self.get_number_glycosylation_sites()+self.get_number_disulfide_bonds()
             return number_of_PTMs
    
    def get_number_modifications(self)->int:
        """
        :return:  Returns the total number of generic modifications found on the protein. \
        if no modification is found it return 0
        :rtype: int
        """
        if not self.has_site_modifications():
            return 0
        else:
            return len(self.get_PTMs_modifications())
    
    def get_number_glycosylation_sites(self)->int:
        """
        :return: The number of glycosylation_sites the protein has, if the protein has no glycosylation sites, the function returns zero
        :rtype: int
        """
        if not self.has_glycosylation_site():
            return 0
        else:
            return len(self.get_PTMs_glycosylation())
    
    def get_number_disulfide_bonds(self)->int:
        """
        :return: The number of disulfide bonds the protein has, if the protein has no disulfide bonds, the function return zero.
        :rtype: int
        """
        if not self.has_disulfide_bond():
            return 0
        else:
            return len(self.get_disulfide_bonds())
        
    def has_sequence_variants(self) ->bool:
        """
        :return:  True if the protein has a sequence variants, and False otherwise.
        :rtype: bool
        """
        if self.extracted_features["SeqVar"] == None:
            return False
        else: 
            return True
    
    def get_sequence_variants(self) -> Dict[str,Dict[str,Union[str,int]]]:
        """
        :return: A dict object that contains all sequence variants within a protein, if the protein has no sequence variants the function returns None.
        :rtype: Dict[str,Dict[str,Union[str,int]]]
        """
        return self.extracted_features["SeqVar"]
    
    def get_number_sequence_variants(self)->int:
        """
        :return: The number of sequence variants the protein has, if the protein has no sequence varient, the function returns 0.
        :rtype: int
        """
        if not self.has_sequence_variants():
            return 0
        else:
            return len(self.get_sequence_variants())
    
    def has_splice_variants(self)->bool:
        """
        :return: True if the sequence has a splice variants and False otherwise.
        :rtype: bool
        """
        if self.extracted_features["SpliceVar"]==None:
            return False
        else:
            return True
    
    def get_splice_variants(self)->Dict[str,Dict[str,Union[str,int]]]:
        """
        :return: A dict object that contains the splice variants. If the protein has no splice variants the function returns None.
        :rtype: Dict[str,Dict[str,Union[str,int]]]
        """
        return self.extracted_features["SpliceVar"]
    
    def get_number_splice_variants(self)->int:
        """
        :return: The number of slice variants in the protein, if the protein has no splice variants, the function returns zero.
        :rtype: int
        """
        if not self.has_splice_variants():
            return 0
        else:
            return len(self.get_splice_variants())
    
    def summary(self)->Dict[str,Union[str,int]]:
        """
        :return: The function return a dict object that summarizes the features of the protein.
        :rtype: Dict[str,Union[str,int]]
        """
        summary=dict()
        summary["has_signal_peptide"]=self.has_signal_peptide()
        summary["number_of_chains"]=self.get_number_chains()
        summary["number_of_domains"]=self.get_number_domains()
        summary["number_of_PTMs"]=self.get_number_PTMs()
        summary["number_of_modifications"]=self.get_number_modifications()
        summary["number_of_glycosylation_sites"]=self.get_number_glycosylation_sites()
        summary["number_of_disulfide_bonds"]=self.get_number_disulfide_bonds()
        summary["number_of_sequence_varients"]=self.get_number_sequence_variants()
        summary["number_of_splice_varients"]=self.get_number_splice_variants()
        summary["number_of_transmembrane_regions"]=self.get_num_transmembrane_regions()
        return summary
    
    def __str__(self)->str:
        """
        The string representation of the class
        Returns
        -------
        the string representation of the class
        """
        summary=self.summary()
        string_rep=""" A protein feature instance with: {}  chains, {} transmembrane regions,  {} domains.  {} PTMs, {} sequence variants and {} splice variants""".format(
            summary["number_of_chains"],summary["number_of_transmembrane_regions"],summary["number_of_domains"],
            summary["number_of_PTMs"],summary["number_of_sequence_varients"],
            summary["number_of_splice_varients"],
            )
        return string_rep

    def __repr__(self)->str: 
        """
        :return: A formated print statement for the class 
        :rtype: str
        """
        return str(self)
    