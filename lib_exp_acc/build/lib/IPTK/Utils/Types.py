#!/usr/bin/env python 
"""Contain a definition of commonly used types through the library 
"""
# load the modules 
from typing import List, Set, Tuple, Dict
import numpy as np
from IPTK.Classes.HLAChain import HLAChain
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
