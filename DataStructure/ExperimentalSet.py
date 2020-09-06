#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: an Experimental set which is a collection of experiments.  ∂
@version: 0.0.1
"""
# load the modules 
from __future__ import annotations
import numpy as np 
import pandas as pd
from IPTK.DataStructure.Experiment import Experiment 
from IPTK.Analysis.AnalysisFunction import get_binnary_peptide_overlap, get_binnary_protein_overlap 
from IPTK.DataStructure.Peptide import Peptide
from typing import Dict, List 
## define some types 
Experiments=Dict[str,Experiment]
Names=List[str]
Counts=Dict[str,int ]
Peptides=List[Peptide]
Proteins=List[str]
# define the class 
class ExperimentSet: 
    def __init__(self,**exp_id_pair)->ExperimentSet:
        """
        @brief: construct an experimentset instance from an arbitrary number of experiment-identifier pairs
        @param exp_id_pair: an arbitrary number of experiments identifier pairs  
        """
        # an experimental dict 
        self._exps=dict()
        # if the number of the provided experiments is zero, we create an an empty instance that will be filled later 
        if len(exp_id_pair)==0: return 
        # get the key and value from the exp-id-pair 
        for ident, exp in exp_id_pair.items():
            if isinstance(exp,Experiment): 
                self._exps[ident]=exp
            else:
                raise ValueError(f"""The constructor expected the value in your input pair to be an instance of class Experiment,
                                     however, your key: {ident} has a type of {type(exp)}""")
    
    def add_experiment(self, **exp_id_pair)->None:
        """
        @brief: add an arbitrary number of experiments to the set 
        @param: exp_id_pair: an arbitrary number of experiments identifier pairs 
        """
        for ident, exp in exp_id_pair.items():
            if isinstance(exp,Experiment): 
                self._exps[ident]=exp
            else:
                raise ValueError(f"""The constructor expected the value in your input pair to be an instance of class Experiment,
                                     however, your key: {ident} has a type of {type(exp)}""")

    def get_num_experiments_in_the_set(self)->int:
        """
        @brief: return the total number of experiments in the set 
        """
        return len(self._exps)
    
    def get_experiments(self)->Dict:
        """
        @brief: return a dict with all the experiment stored in the instance 
        """
        return self._exps
    
    def get_experimental_names(self)->Names:
        """
        @brief: return a list with all the experiments associated names or identifiers  
        """
        return list(self._exps.keys())
    
    def __getitem__(self, name)->Experiment:
        """
        @brief: a magic function for accessing the experiments stored in the set 
        """
        if name in self.get_experimental_names():
            return self._exps[name]
        else:
            raise ValueError(f"Name: {name} is not defined in the current instance")
    
    def __len__(self):
        """
        @brief: a magic function for computing the length of the set, which is the number of experiments 
        stored inside the experiment
        """
        return len(self._exps)
    
    def __str__(self):
        """
        @brief: a magic function for computing a string representation of the class
        """
        return f'an experimental set with {len(self)} Experiments in it'
        
    
    def get_tissue_counts(self)->Counts:
        """
        @brief: return the number of experiment obtained from each tissue in the ExperimentalSet.
        """
        tissues_counts=dict()
        # loop over all the elements in the set 
        for exp in self._exps.keys(): 
            tissue_=self._exps[exp].get_tissue_name()
            if tissue_ in tissues_counts.keys():
                tissues_counts[tissue_]+=1
            else: 
                tissues_counts[tissue_]=1
        # return the results 
        return tissues_counts
    
    def get_allele_count(self)->Counts:
        """
        @brief: return the number of experiment obtained from each allele in the ExperimentalSet.
        """
        allele_counts=dict()
        # loop over all the elements in the set 
        for exp in self._exps.keys(): 
            alleles=self._exps[exp].get_hla_allele()
            for allele in alleles:
                if allele in allele_counts.keys():
                    allele_counts[allele]+=1
                else: 
                    allele_counts[allele]=1
        # return the results 
        return allele_counts
    
    def get_proband_count(self)->Counts:
        """
        @brief: return the number of experiments obtained from each proband in the ExperimentalSet.
        """
        proband_count=dict()
        # loop over all the elements in the set 
        for exp in self._exps.keys(): 
            proband_=self._exps[exp].get_proband_name()
            if proband_ in proband_count.keys():
                proband_count[proband_]+=1
            else: 
                proband_count[proband_]=1
        # return the gene_counts 
        return proband_count
    
    def group_by_tissue(self)->Dict[str,ExperimentSet]:
        """
        @brief: return a map between each tissue and an Experimentalset object representing all experiments 
        belonging to this tissue. 
        """ 
        # define the set of tissues 
        tissues2exps: Dict[str,ExperimentSet]=dict()
        tissue_counter: Dict[str, int]=dict()
        # Initialize the counters 
        tissues: List[str] =list(set([self._exps[exp].get_tissue_name() for exp in self._exps.keys()]))
        for tissue in tissues:
            tissue_counter[tissue]=0
        # loop over all the elements in the set 
        for exp in self._exps.keys(): 
            tissue_=self._exps[exp].get_tissue_name()
            if tissue_ in tissues2exps.keys():
                temp_pair: Dict[str,Experiment]={
                    tissue_+str(tissue_counter[tissue_]):self._exps[exp]
                }
                tissues2exps[tissue_].add_experiment(**temp_pair)
                tissue_counter[tissue_]+=1
            else:
                temp_pair: Dict[str,Experiment]={
                    tissue_+str(tissue_counter[tissue_]):self._exps[exp]
                } 
                tissues2exps[tissue_]=ExperimentSet(**temp_pair)
                tissue_counter[tissue_]+=1
        # return results 
        return tissues2exps
    
    def group_by_proband(self)->Dict[str,ExperimentSet]:
        """
        @brief: return a map between each proband and an Experimentalset object represent all the experiments objects 
        belonging to this proband. 
        """
        # define the set of tissues 
        proband2exps: Dict[str, ExperimentSet]=dict()
        proband_counter: Dict[str,int]=dict() 
        # define the probands 
        probands: List[str] = list(set([self._exps[exp].get_proband_name() for exp in self._exps.keys()]))
        for proband in probands: 
            proband_counter[proband]=0
        # loop over all the elements in the set 
        for exp in self._exps.keys(): 
            proband_=self._exps[exp].get_proband_name()
            if proband_ in proband2exps.keys():
                temp_pair: Dict[str,Experiment]={
                    proband_+str(proband_counter[proband_]):self._exps[exp]
                }
                proband2exps[proband_].add_experiment(**temp_pair)
                proband_counter[proband_]+=1
            else:
                temp_pair: Dict[str,Experiment]={
                    proband_+str(proband_counter[proband_]):self._exps[exp]
                }
                proband2exps[proband_]=ExperimentSet(**temp_pair)
                proband_counter[proband_]+=1
        # return results 
        return proband2exps

    
    def get_unique_peptides(self)->Peptides:
        """
        @brief: compute the set of unique peptides in the experimentalSet 
        @return: A list of all the protein that overlap over the experimentalSet
        """
        res=[]
        for exp_name in self.get_experimental_names(): 
            res.extend(self[exp_name].get_peptides())
        return list(set(res))
    
    def get_unique_proteins(self)->Proteins:
        """
        @brief: compute the set of unique proteins in the experimentalset
        @return: a list of all proteins that overlap over the experimentSet
        """
        res=[]
        for exp in self.get_experimental_names():
            res.extend(self[exp].get_proteins())
        return list(set(res))
    
    def is_peptide_present_in_all(self, peptide:str)->bool:
        """
        @brief: return whether or not a peptide is present in all experiments inside the instance or not 
        @param: peptide: the peptide sequence to search its occurrences in every experiment contained in the set
        """
        # check whether the peptide is defined in the set of unique peptides or not  
        if peptide not in self.get_unique_peptides():
            return False
        # check if the peptide is a member of ALL EXPERIMENTs in the set 
        for exp_name in self.get_experimental_names(): 
            if not self[exp_name].is_member(peptide):
                return False 
        # if the peptide is a member of every experiment in the set, we return True
        return True
    
    def is_protein_present_in_all(self, protein:str)->bool:
        """
        @brief: return whether or not a peptide is present in all experiments inside the set or not 
        @param: protein: the name of the protein to search its occurrences in every experimental in the set 
        """
        if protein not in self.get_unique_proteins():
            return False
        #check that the protein is a member of every experiment 
        for exp_name in self.get_experimental_names():
            if not self[exp_name].is_a_parent_protein(protein):
                return False
        # return true as it is defined in all the experiment in the dataset. 
        return True
    
    
    def get_peptides_present_in_all(self)->Peptides:
        """
        @brief: return the peptides that are observed in every experiments in the set.  
        """
        all_peptides=self.get_unique_peptides()
        results=[]
        for pep in all_peptides: 
            if self.is_peptide_present_in_all(pep):
                results.append(pep)
        return results 
    
    def get_protein_present_in_all(self)->Proteins:
        """
        @brief: return the proteins that are inferred in all experiments of the set 
        """
        all_proteins=self.get_unique_proteins()
        results=[]
        for protein in all_proteins: 
            if self.is_protein_present_in_all(protein):
                results.append(protein)
        return results 
        
    def compute_peptide_overlap_matrix(self)->np.ndarray:
        """
        @brief: return a 2D matrix containing the number of peptide overlapping between each pair 
        of experiments inside the current instance collection of experiments.  
        """
        # allocate the results array 
        results_array=np.zeros(shape=(len(self), len(self))) 
        experiment_names=self.get_experimental_names()
        for raw_idx in range(len(experiment_names)):
            for col_idx in range(len(experiment_names)):
                results_array[raw_idx,col_idx]=len(
                    get_binnary_peptide_overlap(self[experiment_names[raw_idx]],
                                                self[experiment_names[col_idx]])) # compute the peptide overlap    
        # construct a dataframe from the results 
        results_df=pd.DataFrame(results_array)      
        # add the columns and index to the results df 
        results_df.columns=experiment_names
        results_df.index=experiment_names
        # return the results 
        return results_df 

    def compute_protein_overlap_matrix(self)->np.ndarray:
        """
        @brief: return a 2D matrix containing the number of proteins overlapping between each pair 
        of experiments inside the current instance collection of experiment.  
        """
        # allocate the results array 
        results_array=np.zeros(shape=(len(self), len(self))) 
        experiment_names=self.get_experimental_names()
        for raw_idx in range(len(experiment_names)):
            for col_idx in range(len(experiment_names)):
                results_array[raw_idx,col_idx]=len(
                    get_binnary_protein_overlap(self[experiment_names[raw_idx]],
                                                self[experiment_names[col_idx]])) # compute the peptide overlap    
        # construct a dataframe from the results 
        results_df=pd.DataFrame(results_array)      
        # add the columns and index to the results df 
        results_df.columns=experiment_names
        results_df.index=experiment_names
        # return the results 
        return results_df 
    
    def compute_peptide_representation_count(self)->Counts:
        """
        @brief: compute the number of times a peptide was observed accross all experiments in the set 
        """
        results=dict()
        unique_peptides=self.get_unique_peptides()
        # fill the dictionary and initialize the counts to zeros 
        for pep in unique_peptides:
            results[pep]=0
        # loop over all the experiment to count the peptides 
        for peptide in unique_peptides:
            for exp_name in self.get_experimental_names():
                if self[exp_name].is_member(peptide):
                    results[pep]+=1
        # return the results after filling it 
        return results 
    
    def compute_protein_representation_count(self)->Counts:
        """
        @brief: compute the number of times a protein was observed accross all the experiment in the set
        """
        results=dict()
        unique_proteins=self.get_unique_proteins()
        # fill the dictionary and initialize the counts to zeros 
        for prot in unique_proteins:
            results[prot]=0
        # loop over all the experiments to count the peptides 
        for prot in unique_proteins:
            for exp_name in self.get_experimental_names():
                if self[exp_name].is_a_parent_protein(prot):
                    results[prot]+=1
        # return the results after filling it 
        return results  
    
    def compute_protein_coverage_over_the_set(self)->Counts:
        """
        @brief: compute the mapped representation for each protein in the set
        """
        results=dict()
        unique_proteins=self.get_unique_proteins()
        # get a consent representation 
        for prot in unique_proteins:
            for exp_name in self.get_experimental_names():
                if self[exp_name].is_a_parent_protein(prot):
                    if prot in results.keys(): 
                        results[prot]+=self[exp_name].get_mapped_protein(prot)
                    else:
                        results[prot]=self[exp_name].get_mapped_protein(prot)
        # return the results after filling the array 
        return results
    