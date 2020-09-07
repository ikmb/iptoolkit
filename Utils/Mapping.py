#!/usr/bin/env python 
"""
@brief: a submodule that contain function to map different database keys 
@author: Hesham ElAbd
@version: 0.0.1
"""
# load the modules: 
import urllib
from typing import List, Dict
import pandas as pd 
# define the mapping function 
def map_from_uniprot_pdb(uniprots: List[str])-> pd.DataFrame:
    """
    @brief: map from uniprot id to protein data bank identifiers
    @param: uniprot_id: a list of uniprot IDs 
    """
    url: str ='https://www.uniprot.org/uploadlists/'
    # define the query parameters 
    q_params: Dict[str, str]={
        'from': 'ACC+ID', 
        'to': 'PDB_ID',
        'format': 'tab',
        'query': ' '.join(uniprots)
    }
    data: bytes =urllib.parse.urlencode(q_params).encode('utf-8')
    request: urllib.request.Request = urllib.request.Request(url,data)
    # read the request
    with urllib.request.urlopen(request) as input_file: 
        results: str =input_file.read().decode('utf-8')
    # parse the resulting strings 
    mapped_pairs: List[str] = results.split('\n')
    # allocate to lists to hold the results 
    unitpot_ids: List[str] = []
    pdb_ids: List[str] = []
    # parse the results 
    for pair in mapped_pairs:
        temp_lists: List[str] = pair.split('\t')
        if len(temp_lists) ==2:
            unitpot_ids.append(temp_lists[0])
            pdb_ids.append(temp_lists[1])
    # combine the data into a dataframe 
    results: pd.DataFrame = pd.DataFrame({
        'Uniprot-ID':unitpot_ids,
        'PDB':pdb_ids
    })
    # return the results 
    return results

def map_from_uniprot_gene(uniprots: List[str])->pd.DataFrame: 
    """
    @brief: map from uniprot id to ensemble gene ids
    @param: uniprot_id: a list of uniprot IDs 
    """
    url: str ='https://www.uniprot.org/uploadlists/'
    # define the query parameters 
    q_params: Dict[str, str]={
        'from': 'ACC+ID', 
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': ' '.join(uniprots)
    }
    data: bytes =urllib.parse.urlencode(q_params).encode('utf-8')
    request: urllib.request.Request = urllib.request.Request(url,data)
    # read the request
    with urllib.request.urlopen(request) as input_file: 
        results: str =input_file.read().decode('utf-8')
    # parse the resulting strings 
    mapped_pairs: List[str] = results.split('\n')
    # allocate to lists to hold the results 
    unitpot_ids: List[str] = []
    ensemble_ids: List[str] = []
    # parse the results 
    for pair in mapped_pairs:
        temp_lists: List[str] = pair.split('\t')
        if len(temp_lists) ==2:
            unitpot_ids.append(temp_lists[0])
            ensemble_ids.append(temp_lists[1])
    # combine the data into a dataframe 
    results: pd.DataFrame = pd.DataFrame({
        'Uniprot-ID':unitpot_ids,
        'PDB':ensemble_ids
    })
    # return the results 
    return results