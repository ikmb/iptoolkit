#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: contain a defination of commenly used types through the library 
@version: 0.0.1
"""
# load the modules 
from typing import List, Set, Tuple, Dict
import numpy as np
from IPTK.DataStructure.HLAChain import HLAChain
# define the types
Sequences = List[str]
MappedProtein = List[np.ndarray]
Proteins=Set[str]
Range=Tuple[int, int]
MappedProteins=Dict[str,MappedProtein]
FastaSet=Dict[str,str]
ProteinHits=List[str]
PlottingKeywards=Dict[str,str]
MappedProteinRepresentation=Dict[str,np.ndarray]
Index=List[int]
HLA_Chains=Dict[str,HLAChain]
Genes=List[str]
AlleleGroup=List[str]
ProteinGroup=List[str]
HLA_Names=List[str]
ProteinSource=Dict[str,str]
Organisms=List[str]