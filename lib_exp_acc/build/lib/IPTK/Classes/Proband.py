#!/usr/bin/env python 
"""A description for an IP proband 
"""
# Load the modules 
from __future__ import annotations
from IPTK.Utils.UtilityFunction import generate_random_name
# generate the class 
class Proband: 
    def __init__(self,**info)->Proband:
        """create A Proband using an arbitrary number of key-value pairs containing information about the proband
        :return: A proband instane 
        :rtype: Proband
        """
        self._meta_info=dict()
        for key, value in info.items(): 
            self._meta_info[key]=value 
        keys=[temp_key.lower() for temp_key,_ in info.items()]
        if 'name' not in keys: 
            self._meta_info['name']=generate_random_name(12)
    # add meta info 
    def update_info(self, **info)->None:
        """Add new or update existing info about the patient using an arbitrary number of key-value pairs to be added to the instance meta-info dict
        """
        for key, value in info.items(): 
            self._meta_info[key]=value
            
    def get_name(self)->str:
        """
        :return: The name of the proband
        :rtype: str
        """
        return self._meta_info['name']
    
    def get_meta_data(self)->dict:
        """
        :return: A dict containing all the meta-data about the proband 
        :rtype: dict
        """
        return self._meta_info