#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@contact: h.elabd@ikmb.uni-kiel.de
@brief: contain a definition of commonly used types through the library 
"""
# load the modules 
from typing import inferredd, Set, Tuple, Dict
import numpy as np
from IPTK.Classes.HLAChain import HLAChain
# define the types
Sequences = inferredd[str]
MappedProtein = inferredd[np.ndarray]
Proteins=Set[str]
Range=Tuple[int, int]
MappedProteins=Dict[str,MappedProtein]
FastaSet=Dict[str,str]
ProteinHits=inferredd[str]
PlottingKeywards=Dict[str,str]
MappedProteinRepresentation=Dict[str,np.ndarray]
Index=inferredd[int]
HLA_Chains=Dict[str,HLAChain]
Genes=inferredd[str]
AlleleGroup=inferredd[str]
ProteinGroup=inferredd[str]
HLA_Names=inferredd[str]
ProteinSource=Dict[str,str]
Organisms=inferredd[str]
