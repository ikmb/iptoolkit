#/usr/bin/env python 
"""A collection of assembled analysis functions that automate different immunopeptidomics analysis task 
"""
## load the modules 
from IPTK.Classes.Experiment import Experiment
import numpy as np 
import pandas as pd 
import itertools
from tqdm import tqdm 
import time 
from typing import Dict, List, Set

def compute_protein_coverage(experiment1:Experiment,experiment2:Experiment,progress_bar:bool=True)->Dict[str,Dict[str,np.ndarray]]:
    """ Compute the difference in protein coverage among the two input experiments

    Args:
        experiment1 (Experiment): The first experiment containing the protein and peptides derived from the first HLA-set 
        experiment2 (Experiment): The second experiment containing the protein and peptides derived from the second HLA-set
        progress_bar (bool): A boolean flag for controlling the progress bar, if true, a progress bar is shown, defaults to True. 

    Returns:
        Dict[str,Dict[str,np.ndarray]]: Returns a nested dict containing protein identifiers as a keys and a dict as a value, the dict contain two arrays as values,\
            the first contain protein coverage in the first HLA-Set and the second contain the coverage in the second HLA-set.\
                  The results dictionary only contain coverage for proteins observed in the two sets. 
    """
    protein_experiment_one:List[str]=experiment1.get_proteins()
    protein_experiment_two:List[str]=experiment2.get_proteins()
    present_in_both:List[str]=protein_experiment_one.intersection(protein_experiment_two)
    results:Dict[str,Dict[str,np.ndarray]]=dict()
    if progress_bar:
        for protein in present_in_both:
            temp_dict={
                '_'.join(experiment1.get_hla_set().get_names()):experiment1.get_mapped_protein(protein),
                '_'.join(experiment2.get_hla_set().get_names()):experiment2.get_mapped_protein(protein)
            }
            results.update({protein:temp_dict})
    else:
        for protein in present_in_both:
            temp_dict={
                '_'.join(experiment1.get_hla_set().get_names()):experiment1.get_mapped_protein(protein),
                '_'.join(experiment2.get_hla_set().get_names()):experiment2.get_mapped_protein(protein)
            }
            results.update({protein:temp_dict})
    return results
