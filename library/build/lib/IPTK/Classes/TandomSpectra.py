#!/usr/bin/env 
""" TandomSpectra is a simple class to link MS1 and MS2 spectra 
"""
## LOAD THE MODULES
#------------------
from __future__ import annotations
import pyopenms as poms 
## DEFINe THE CLASS
#------------------
class TandomSpectra:
    def __init__(self,parent:poms.pyopenms_4.MSExperiment,
        children:List[poms.pyopenms_4.MSExperiment])->TandomSpectra:
        """class constructor, construct a new instance of class TandomSpectra

        :param parent: an MS1 spectrum that acts as the precursor for the class  
        :type parent: poms.pyopenms_4.MSExperiment
        :param children: MS2 spectra obtained from the precursor MS1 spectrum 
        :type children: a list of MS2 spectrums 
        :return: a new instance 
        :rtype: TandomSpectra
        """
        if parent.getMSLevel()!=1:
            raise ValueError(f'Parent Spectrum must have an MS-level of 1, however, current parent has: {parent.getMSLevel()}')
        else: 
            self._parent=parent
        for child in children:
            if child.getMSLevel() != 2:
                raise ValueError(f'All the children spectra must be of level 2, however, spectra: {child.getNativeID()} has a level of: {child.getMSLevel()}')
        self._children=children
        return 
    
    def __len__(self)->int:
        """a magic function to return the number of children in the class 

        :return: [description]
        :rtype: int
        """
        return len(self._children)
    
    def __repr__(self)->str:
        return str(self)

    def __str__(self)->str:
        return f'a TandomSpectra with {len(self)} children'
        