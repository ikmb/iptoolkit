#!/usr/bin/env python 
""" The class provide some convenient function for computing summary stats and extracting information from\
     a collection of mzML files  
""" 
## IMPORT THE MODULES
#--------------------
from __future__ import annotations
import os
import pyopenms as poms 
from IPTK.Classes.MzMLExperiment import MzMLExperiment
import pandas as pd
from typing import List, Callable
from tqdm import tqdm 
## CLASS DEFINTION
#-----------------
class MzMLExperimentSet:
    def __init__(self, path2files:str)->MzMLExperimentSet:
        """ create a new instance of mzMLExperimentSet from all mzML files present in a provided directory 

        :param path2files: The path to a directory containing all the files that need to be loaded
        :type path2files: str
        :raises ValueError: incase the provided path does not exists 
        :raises RuntimeError: incase the creating an MzMLExperiment instance from one of the files failed
        :return: an MzMLExperimentSet instance 
        :rtype: MzMLExperimentSet
        """
        if not os.path.exists(path2files):
            raise ValueError(f'The provided path: {path2files} sequences does not exists!!')
        mzML_files=[p_file for p_file in os.listdir(path2files) if '.mzML' in p_file]
        self._files={}
        for mzML_file in tqdm(mzML_files): 
            try:
                self._files[mzML_file.split('.')[0]]=MzMLExperiment(os.path.join(path2files,mzML_file))
            except Exception as exp:
                raise RuntimeError(f'While loading: {os.path.join(path2files,mzML_file)} the following error was encountered: {str(exp)}')
        return
    
    def get_total_number_of_spectra_per_file(self)->pd.DataFrame:
        """return the number of spectra in each MS file stored in the set of MzML experiment  

        :return: a table with two columns, the first, is the experiment name and the second is the column name 
        :rtype: pd.DataFrame
        """
        # get the count of MS spectra in each file
        counts_dict=dict()
        for p_file in tqdm(self._files.keys()):
            counts_dict[p_file]=len(self._files[p_file])
        # construct a pandas dataframe from the dict
        table=pd.DataFrame().from_dict(counts_dict,orient='index')
        table['e']=table.index
        table.columns=['Count','Experiment Name']
        table.reset_index(drop=True,inplace=True)
        table.reindex(columns=['Experiment Name','Count'],inplace=True)
        return table
    
    def get_number_of_MS1_per_file(self)->pd.DataFrame:
        """return the number of MS2 spectra in each MS file stored in the set of MzML experiment 

        :return: a table containing the count from each file  
        :rtype: pd.DataFrame
        """
        counts_dict=dict()
        for p_file in tqdm(self._files.keys()):
            counts_dict[p_file]=self._files[p_file].get_num_MS2_spectra()
        # construct a pandas dataframe from the dict
        table=pd.DataFrame().from_dict(counts_dict,orient='index')
        table['e']=table.index
        table.columns=['Count','Experiment Name']
        table.reset_index(drop=True,inplace=True)
        table.reindex(columns=['Experiment Name','Count'],inplace=True)
        return table

    def get_number_of_MS2_per_file(self)->pd.DataFrame:
        """return the number of MS2 spectra in each MS file stored in the set of MzML experiment 

        :return: a table containing the count from each file  
        :rtype: pd.DataFrame
        """
        counts_dict=dict()
        for p_file in tqdm(self._files.keys()):
            counts_dict[p_file]=self._files[p_file].get_num_MS2_spectra()
        # construct a pandas dataframe from the dict
        table=pd.DataFrame().from_dict(counts_dict,orient='index')
        table['e']=table.index
        table.columns=['Count','Experiment Name']
        table.reset_index(drop=True,inplace=True)
        table.reindex(columns=['Experiment Name','Count'],inplace=True)
        return table

