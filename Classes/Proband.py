#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@contact: h.elabd@ikmb.uni-kiel.de
@brief: A description for an IP proband 
"""
# Load the modules 
from __future__ import annotations
from IPTK.Utils.UtilityFunction import generate_random_name
# generate the class 
class Proband: 
    def __init__(self,**info)->Proband:
        """
        @brief initialize a Proband using information provided in the initializer 
        @param: info: an arbitrary number of key-value pair containing information about the proband
        """
        self._meta_info=dict()
        for key, value in info.items(): 
            self._meta_info[key]=value 
        keys=[temp_key.lower() for temp_key,_ in info.items()]
        if 'name' not in keys: 
            self._meta_info['name']=generate_random_name(12)
    # add meta info 
    def update_info(self, **info)->None:
        """
        @brief: add new or update existing info about the patient 
        @param: info an arbitatry number of key-value pair to be added to added to the instance meta-info dict
        """
        for key, value in info.items(): 
            self._meta_info[key]=value
            
    def get_name(self)->str:
        """
        @brief: return the name of the patient
        """
        return self._meta_info['name']
    
    def get_meta_data(self)->dict:
        """
        @brief: return a dict that contain all the meta-data about the patient 
        """
        return self._meta_info