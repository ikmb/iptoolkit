#!/usr/bin/env python 
"""An abstraction for an MzML file, the class relay heavily on the amazing library\
OpenMS and its python binding PyOpenMS
"""
## LOAD THE MODULES
#------------------
from __future__ import annotations
import pyopenms as poms 
import numpy as np 
import pandas as pd 
from typing import List, Dict  
## CLASS DEFINITION
#------------------
class MzMLExperiment:
    def __init__(self,path2file: str)->MzMLExperiment:
        """class constructor, constructs an MzML object from an MzML file
        
        :param path2file: the path to the mzML file 
        :type path2file: str
        :raises IOError: In case loading the file failed 
        :return: a new object constructed from the provided file 
        :rtype: MzMLExperiment
        """
        self.exp=poms.MSExperiment()
        self._spectra_tree=None
        try:
            poms.MzMLFile().load(path2file,self.exp)
        except RuntimeError as exp: 
            raise IOError(f'While loading the input file: {path2file} the following error was encountered: {exp}')
    
    def get_num_MS1_spectra(self)->int:
        """ 
        :return: the number of MS1 spectra in the current experiment  
        :rtype: int
        """
        return len([spectra for spectra in self.exp.getSpectra() if spectra.getMSLevel()==1]) 

    def get_num_MS2_spectra(self)->int:
        """ 
        :return: the number of MS2 spectra in the current experiment  
        :rtype: int
        """ 
        return len([spectra for spectra in self.exp.getSpectra() if spectra.getMSLevel()==2]) 
        

    def get_MS1_spectra(self)->List[poms.pyopenms_2.MSSpectrum]:
        """
        :return: a list of the MS1 spectra in the mzML file 
        :rtype: List[poms.pyopenms_2.MSSpectrum]
        """ 
        return [spectra for spectra in self.exp.getSpectra() if spectra.getMSLevel()==1]

    def get_MS2_spectra(self):
        """
        :return: a list of the MS2 spectra in the mzML file 
        :rtype: List[poms.pyopenms_2.MSSpectrum]
        """
        return [spectra for spectra in self.exp.getSpectra() if spectra.getMSLevel()==1]

    def __len__(self)->int:
        """a magic function to return the total number of spectra in the provided experiment  

        :return: the number of generated spectra in the file 
        :rtype: int
        """
        return self.exp.getNrSpectra()
    
    def __getitem__(self, idx:int)->poms.pyopenms_2.MSSpectrum:
        """access the ith spectrum in the file 

        :param idx: the index of the spectrum 
        :type idx: int
        :raises IndexError: Incase an the provided index is negative of bigger than the number of the spectra in the file 
        :return: the corresponding spectrum 
        :rtype: poms.pyopenms_2.MSSpectrum
        """
        try:
            return self.exp.getSpectrum(idx)
        except Exception as exp: 
            raise IndexError(f'While getting the {idx}th spectrum the following error was encountered: {exp}')
    
    def get_exp(self)->poms.pyopenms_4.MSExperiment:
        """
        :return: return the instance experiment objects 
        :rtype: poms.pyopenms_4.MSExperiment
        """
        return self.exp

    def get_MS1_RT(self)->np.ndarray:
        """return the Retention time of all MS1 spectra in the file 
        :return: a NumPy array containing the RT for all MS1 spectra in the provided file 
        :rtype: np.ndarray
        """
        return np.array([spectra.getRT() for spectra in self.get_MS1_spectra()])

    def get_MS2_RT(self)->np.ndarray:
        """return the Retention time of all MS2 spectra in the file 
        :return: a NumPy array containing the RT for all MS2 spectra in the provided file 
        :rtype: np.ndarray
        """
        return np.array([spectra.getRT() for spectra in self.get_MS2_spectra()])

   
        










