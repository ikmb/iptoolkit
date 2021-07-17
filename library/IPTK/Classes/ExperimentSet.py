#!/usr/bin/env python 
"""An ExperimentSet which is a collection of experiments.
The class provides an API for integrating and comparing different experimental instances.   
"""
# load the modules 
from __future__ import annotations
import time
from tqdm import tqdm 
import numpy as np 
import pandas as pd
from IPTK.Classes.Experiment import Experiment 
from IPTK.Classes.Peptide import Peptide
from IPTK.Analysis.AnalysisFunction import (get_binnary_peptide_overlap, get_binnary_protein_overlap, 
    compute_jaccard_index, compute_change_in_protein_representation,compute_expression_correlation)
from typing import Dict, List
## define some types 
Experiments=Dict[str,Experiment]
Names=List[str]
Counts=Dict[str,int ]
Peptides=List[Peptide]
Proteins=List[str]
# define the class 
class ExperimentSet: 
    """an API for integrating and comparing different experimental instances
    """
    def __init__(self,**exp_id_pair)->ExperimentSet:
        """Create a new instance using an arbitrary number of experiments. 
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
        """adds an arbitrary number of experiments to the set 
        """
        for ident, exp in exp_id_pair.items():
            if isinstance(exp,Experiment): 
                self._exps[ident]=exp
            else:
                raise ValueError(f"""The constructor expected the value in your input pair to be an instance of class Experiment,
                                     however, your key: {ident} has a type of {type(exp)}""")

    def get_num_experiments_in_the_set(self)->int:
        """
        :return: The number of experiments currently in the set 
        :rtype: int
        """
        return len(self._exps)
    
    def get_experiments(self)->Dict[Experiment]:
        """
        :return: returns a dict with all the experiments stored in the instance as value and ids as keys. 
        :rtype: Dict[Experiment]
        """
        return self._exps
    
    def get_experiment(self, exp_name:str)->Experiment:
        """ returns the experiment pointed to by the provided experimental name

        :param exp_name: the name of the experiment 
        :type exp_name: str
        :raises KeyError: if the provided experimental name is not in the dataset. 
        :return: the experiment corresponds to the provided name 
        :rtype: Experiment
        """
        # define the experiment name 
        if exp_name not in self._exps.keys():
            raise KeyError(f"The provided experimental name: {exp_name} is not defined in the current experiment") 
        # return the experiment 
        return self._exps[exp_name]

    def get_experimental_names(self)->Names:
        """
        :return: A list with all the identifiers of the experiments in the set  
        :rtype: Names
        """
        return list(self._exps.keys())
    
    def get_unique_orgs(self)->List[str]:
        """
        :return: A list of the unique organisms in the set
        :rtype: List[str]
        """
        unique_orgs: List[str] = []
        for name in self.get_experimental_names():
             unique_orgs.extend(self._exps[name].get_orgs())
        # compute the unique organisms in the list of organisms 
        unique_orgs=list(set(unique_orgs))
        return unique_orgs

    def get_total_peptide_per_org_count(self) ->pd.DataFrame:
        """
        :return: The total count of peptides per organism accross the all experiments in the set. 
        :rtype: pd.DataFrame
        """
        # first, get the unique organisms in the set.  
        unique_orgs: List[str] = self.get_unique_orgs()
        # create a counter and initialize it to zero to hold the results 
        org_counter=dict()
        for org in unique_orgs:
            org_counter[org]=0
        # update the counts 
        for name in tqdm(self.get_experimental_names()):
            for _, row in self._exps[name].get_peptides_per_organism().iterrows():
                 org_counter[row['Organisms']]+=row['Counts']
       	# make the data compatible with data frames 
        for org in org_counter.keys():
            org_counter[org]=[org_counter[org]]
		# create a dataframe
        res: pd.DataFrame = pd.DataFrame(org_counter).T
		# add the index as an extra-columns 
        res['Organisms'] = res.index.tolist()
        res.columns=['Counts','Organisms']
        # reformat the dataframe 
        res.reset_index(drop=True,inplace=True)
        res=res.reindex(columns=['Organisms','Counts'])
        # sort the results 
        res=res.sort_values(by='Counts',ascending=False)
        return res

    def compare_org_count_among_exps(self, org:str, abs_count: bool =False) ->pd.DataFrame: 
        """
        :param org: The name of the organism to query the database for it. 
        :type org: str
        :param abs_count: The absolute count, defaults to False
        :type abs_count: bool, optional
        :return: The count of the peptides that belong to a specific organism in the database.
        :rtype: pd.DataFrame
        """
        # allocate and array to hold the results 
        res: np.ndarray = np.zeros((len(self),len(self)))
        # loop over all the experiments in the set 
        # initialize the counters 
        experimental_name: List[str] = self.get_experimental_names()
        for row_idx in tqdm(range(len(experimental_name))):
            # get the counts per column 
            org_row: pd.DataFrame = self._exps[experimental_name[row_idx]].get_peptides_per_organism()
            org_row_count: int = org_row.loc[org_row.iloc[:,0]==org]['Counts'] 
            if org_row_count.empty:
                org_row_count=0
            else:
                org_row_count=org_row_count.tolist()[0]
            # get the counts per row 
            for col_idx in range(len(experimental_name)):
                # get the row-column
                org_col: pd.DataFrame = self._exps[experimental_name[col_idx]].get_peptides_per_organism()
                org_col_count: int = org_col.loc[org_col.iloc[:,0]==org]['Counts'] 
                if org_col_count.empty: 
                    org_col_count=0
                else: 
                    org_col_count=org_col_count.tolist()[0]
                # add the count to the
                res[row_idx,col_idx]= org_row_count-org_col_count
        # return a data frame of the results 
        res=pd.DataFrame(res)
        # add the name of columns and index 
        res.columns=experimental_name
        res.index=experimental_name
        # return the results
        return res 
    
    def drop_peptides_belong_to_org(self, org_name: str) -> None:
        """drops all the peptides that belong to the provided organisms from all experiments in the set.


        :param org_name: The name of the organism to drop.
        :type org_name: str
        """
        for name in self.get_experimental_names():
            self[name].drop_peptide_belong_to_org(org_name)
        return 

    def __getitem__(self, name:str)->Experiment:
        """ A magic function for accessing the experiments stored in the set

        :param name: The experiment name or id 
        :type name: str
        :raises KeyError: if the provided name is not defined in the current instance. 
        :return: the experiment with the corresponding name
        :rtype: Experiment
        """
        if name in self.get_experimental_names():
            return self._exps[name]
        else:
            raise KeyError(f"Name: {name} is not defined in the current instance")
    
    def __len__(self)->int:
        """  A magic function for computing the length of the set, which is the number of experiments 
        stored inside the instance 
        
        :return: The number of experiments stored inside the instance 
        :rtype: int 
        """
        return len(self._exps)
    
    def __str__(self)->str:
        """ A string representation of the class

        :return: A string representation of the class
        :rtype: str
        """
        return f'an experimental set with {len(self)} Experiments in it.'
        
    
    def get_tissue_counts(self)->Counts:
        """
        :return: The number of experiments obtained from each tissue in the current instance
        :rtype: Counts
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
        :return: The number of experiments obtained from each allele in the instance.
        :rtype: Counts
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
        :return: The number of experiments obtained from each proband in the ExperimentSet.
        :rtype: Counts
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
        :return:  A map between each tissue and an ExperimentSet object representing all experiments belonging to that tissue. 
        :rtype: Dict[str,ExperimentSet]
        """
        # define the set of tissues 
        tissues2exps: Dict[str,ExperimentSet]=dict()
        tissue_counter: Dict[str, int]=dict()
        # Initialize the counters 
        tissues: List[str] =list(set([self._exps[exp].get_tissue_name() for exp in self._exps.keys()]))
        for tissue in tissues:
            tissue_counter[tissue]=0
        # loop over all the elements in the set 
        for exp in tqdm(self._exps.keys()): 
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
        :return: a map between each proband and an Experimentset object representing all experiments objects \
        belonging to each unique proband in the set. 
        :rtype: Dict[str,ExperimentSet]
        """
        # define the set of tissues 
        proband2exps: Dict[str, ExperimentSet]=dict()
        proband_counter: Dict[str,int]=dict() 
        # define the probands 
        probands: List[str] = list(set([self._exps[exp].get_proband_name() for exp in self._exps.keys()]))
        for proband in probands: 
            proband_counter[proband]=0
        # loop over all the elements in the set 
        for exp in tqdm(self._exps.keys()): 
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
        :return: The set of unique peptides in the experimentSet 
        :rtype: Peptides
        """
        res:List[str]=[]
        for exp_name in self.get_experimental_names(): 
            res.extend(self[exp_name].get_peptides())
        return list(set(res))
    
    def get_unique_proteins(self)->Proteins:
        """
        :return: The set of unique proteins in the experimentset
        :rtype: Proteins
        """
        res:List[str]=[]
        for exp in self.get_experimental_names():
            res.extend(self[exp].get_proteins())
        return list(set(res))
    
    def is_peptide_present_in_all(self, peptide:str)->bool:
        """
        :param peptide: The peptide sequence to search its occurrences in every experiment contained in the set
        :type peptide: str
        :return: True if peptide is present in all experiments inside the instance, False otherwise 
        :rtype: bool
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
        :param protein: the protein id to search its occurrences in every experimental in the set 
        :type protein: str
        :return: True if peptide is present in all experiments inside the instance, False otherwise 
        :rtype: bool
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
        :return: The peptides that are observed in every experiment in the set.
        :rtype: Peptides
        """
        all_peptides=self.get_unique_peptides()
        results=[]
        for pep in all_peptides: 
            if self.is_peptide_present_in_all(pep):
                results.append(pep)
        return results 
    
    def get_proteins_present_in_all(self)->Proteins:
        """
        :return: The proteins that are inferred in all experiments of the set. 
        :rtype: Proteins
        """
        all_proteins=self.get_unique_proteins()
        results: List[str]=[]
        for protein in all_proteins: 
            if self.is_protein_present_in_all(protein):
                results.append(protein)
        return results 
        
    def compute_peptide_overlap_matrix(self, method:str='count')->np.ndarray:
        """

        :return: A 2D matrix containing the number of peptide overlapping between each pair of experiments inside the current instance collection of experiment.   
        :rtype: np.ndarray
        :param method: The method of computing the overlap methods, can be count or jaccard, incase of count, exact matches between samples are used, i.e. number of overlaps,\
             meanwhile, incase of jaccard, jaccard Index between each pair of sample is used.  
        :type protein: str
        """
        # validate input method: 
        #----------------------
        if method!='count' and method!='jaccard':
            raise ValueError(f"Method: {method} is not supported!, only count and jaccard are currently supported.")
        # allocate the results array 
        results_array=np.zeros(shape=(len(self), len(self))) 
        experiment_names=self.get_experimental_names()
        for raw_idx in tqdm(range(len(experiment_names))):
            for col_idx in range(len(experiment_names)):
                if method=='count':
                    results_array[raw_idx,col_idx]=len(get_binnary_peptide_overlap(self[experiment_names[raw_idx]],self[experiment_names[col_idx]])) # compute the peptide overlap    
                elif method=='jaccard':
                    results_array[raw_idx,col_idx]=compute_jaccard_index(self[experiment_names[raw_idx]],self[experiment_names[col_idx]],level='peptide')
        # construct a dataframe from the results 
        results_df=pd.DataFrame(results_array)      
        # add the columns and index to the results df 
        results_df.columns=experiment_names
        results_df.index=experiment_names
        # return the results 
        return results_df 

    def compute_protein_overlap_matrix(self,method:str='count')->pd.DataFrame:
        """
        :return: returns a 2D matrix containing the number of proteins overlapping between each pair of experiments inside the current instance collection of experiment.  
        :rtype: np.ndarray
        :param method: The method of computing the overlap methods, can be count or jaccard, incase of count, exact matches between samples are used, i.e. number of overlaps,\
             meanwhile, incase of jaccard, jaccard Index between each pair of sample is used.  
        :type protein: str
        """
        # validate input method: 
        #----------------------
        if method!='count' and method!='jaccard':
            raise ValueError(f"Method: {method} is not supported!, only count and jaccard are currently supported.")
        # allocate the results array 
        results_array=np.zeros(shape=(len(self), len(self))) 
        experiment_names=self.get_experimental_names()
        for raw_idx in tqdm(range(len(experiment_names))):
            for col_idx in range(len(experiment_names)):
                if method=='count':
                    results_array[raw_idx,col_idx]=len(get_binnary_protein_overlap(self[experiment_names[raw_idx]],self[experiment_names[col_idx]])) # compute the peptide overlap    
                elif method=='jaccard':
                    results_array[raw_idx,col_idx]=compute_jaccard_index(self[experiment_names[raw_idx]],self[experiment_names[col_idx]],level='protein')
        # construct a dataframe from the results 
        results_df=pd.DataFrame(results_array)      
        # add the columns and index to the results df 
        results_df.columns=experiment_names
        results_df.index=experiment_names
        # return the results 
        return results_df 
    
    def compute_peptide_representation_count(self)->Counts: 
        """
        :return: The number of times a peptide was observed accross experiments in the set 
        :rtype: Counts
        """
        # define the results object 
        results: Dict[str, int]=dict()
        unique_peptides: List[str]=self.get_unique_peptides()
        # fill the dictionary and initialize the counts to zeros 
        for peptide in unique_peptides:
            results[peptide]=0
        # loop over all the experiment to count the peptides 
        for peptide in tqdm(unique_peptides):
            for exp_name in self.get_experimental_names():
                if self[exp_name].is_member(peptide):
                    results[peptide]+=1
        # return the results after filling it 
        return results 
    
    def compute_protein_representation_count(self)->Counts: 
        """
        :return: The number of times a protein was observed accross the experiment in the set
        :rtype: Counts
        """
        results=dict()
        unique_proteins=self.get_unique_proteins()
        # fill the dictionary and initialize the counts to zeros 
        for prot in unique_proteins:
            results[prot]=0
        # loop over all the experiments to count the peptides 
        for prot in tqdm(unique_proteins):
            for exp_name in self.get_experimental_names():
                if self[exp_name].is_a_parent_protein(prot):
                    results[prot]+=1
        # return the results after filling it 
        return results  
    
    def compute_protein_coverage_over_the_set(self)->Dict[str, np.ndarray]:
        """
        :return: The mapped representation for each protein accross the entire set
        :rtype: Dict[str, np.ndarray]
        """
        results:Dict[str, np.ndarray]=dict()
        unique_proteins=self.get_unique_proteins()
        # get a consent representation 
        for prot in tqdm(unique_proteins):
            for exp_name in self.get_experimental_names():
                if self[exp_name].is_a_parent_protein(prot):
                    if prot in results.keys(): 
                        results[prot]+=self[exp_name].get_mapped_protein(prot)
                    else:
                        results[prot]=self[exp_name].get_mapped_protein(prot)
        # return the results after filling the array 
        return results


    def compute_correlation_in_experssion(self)->pd.DataFrame:
        """computes the correlation in parent protein gene-expression across all the experiments
        in the set. See the function **compute_binary_correlation** in the analysis module 
        for information about the computational logic.

        :return: returns a 2D matrix containing the coorelation in gene expression between each pair of experiments inside the current instance collection of experiments.  
        :rtype: pd.DataFrame
        """
        # get the experimental index
        exps_ids: List[str] = self.get_experimental_names()  
        # allocate an array to hold the results
        res_array: np.ndarray = np.zeros((len(exps_ids),len(exps_ids)))
        # fill the array with the expression correlation 
        for row_idx in tqdm(range(len(exps_ids))): 
            for col_idx in range(len(exps_ids)):
                res_array[row_idx,col_idx] = compute_expression_correlation(
                    self.get_experiment(exps_ids[row_idx]),
                    self.get_experiment(exps_ids[col_idx])
                )
        # create a dataframe 
        res: pd.DataFrame =pd.DataFrame(res_array)
        # add the col and row names 
        res.columns=exps_ids
        res.index=exps_ids
        # return the results 
        return res 

    def compute_change_in_protein_representation(self)->np.ndarray:
        """Computes the change in protein representation among the proteins that are presented/ detected in all of the \
        instance's experiments.
        
        :returns: a 3D tensor, T, with shape of (num-experiments, num-experiments, num-proteins), \
        where T[i,j,k] is a the difference between experiment i & j with respect to the kth protein \
        :rtype: np.ndarray
        """
        # get the number of experiments and proteins 
        present_in_all: List[str] = self.get_proteins_present_in_all()
        num_exps: int = self.get_num_experiments_in_the_set()
        # allocate an array to hold the results 
        results_array: np.ndarray = np.zeros(shape=(num_exps, num_exps,len(present_in_all)))
        ## fill the array 
        # create some counters 
        col_counter: int = 0
        row_counter: int= 0
        for prod_idx in tqdm(range(len(present_in_all))): 
            for exp_col in self.get_experiments().keys():
                for exp_row in self.get_experiments().keys():
                    results_array[row_counter,col_counter,prod_idx]=compute_change_in_protein_representation(
                            self.get_experiment(exp_row).get_mapped_protein(present_in_all[prod_idx]), 
                            self.get_experiment(exp_col).get_mapped_protein(present_in_all[prod_idx])
                        )
                    # increase the row counters 
                    row_counter+=1
                row_counter=0
                # increase the columns counter
                col_counter+=1 
            col_counter=0
        # return the results
        return results_array

    def compute_average_distance_between_exps(self)->pd.DataFrame: 
        """computes the average distance between experiments by taking the average over the z-axis
        of the 3D tensor computed by the function compute_change_in_protein_representation.
        
        :return:  A 2D tensor with shape of (num-experiments, num-experiments)
        :rtype: pd.DataFrame
        """
        diff_protein_overlap: np.ndarray = self.compute_change_in_protein_representation()
        # average overlap is: 
        average_scores: np.ndarray = np.sum(diff_protein_overlap,axis=-1)
        # construct a dataframe of the results 
        results_df: pd.DataFrame = pd.DataFrame(average_scores)
        # add the colnames 
        results_df.columns=list(self.get_experiments().keys())
        results_df.index=list(self.get_experiments().keys())
        # return the results 
        return results_df
    
    def compare_peptide_counts(self)->pd.DataFrame:
        """
        :return: A table that contain the total number of peptides and per-organism peptide counts \
        among all experiments in the set
        :rtype: pd.DataFrame
        """
        # allocate a list to hold the results 
        total_count: List[str]= []
        experiment_name: List[str]= []
        # fill the lists 
        for name in self.get_experimental_names():
            experiment_name.append(name)
            total_count.append(len(self[name].get_peptides()))
        # create a dataframe to store the results 
        res: pd.DataFrame = pd.DataFrame({
            'Organisms':['Total']*len(total_count), 
            'Counts':total_count,
            'Names':experiment_name
        })
        # get the count per organism in each experiment 
        for exp_name in self.get_experimental_names():
            # get the temp results 
            temp_res: pd.DataFrame = self[exp_name].get_peptides_per_organism()
            # add the experiment name as an extra column 
            temp_res['Names'] = [exp_name] * temp_res.shape[0]
            # concatenate the results
            res= pd.concat([res,temp_res], axis=0)
        # return the results 
        return res
    
    def compute_peptide_length_table(self)->pd.DataFrame:
        """
        :return: A table that contain the length of each peptide in the experiment
        :rtype: pd.DataFrame
        """
        # allocate a data frame to hold the results 
        res: pd.DataFrame = pd.DataFrame(columns=['Peptide_length','Names'])
        # loop over all the experiment
        for exp_name in self.get_experimental_names():
            temp_res: pd.DataFrame = pd.DataFrame({
                'Peptide_length': self[exp_name].get_peptides_length(),
                'Names': [exp_name]*len(self[exp_name].get_peptides_length())
            }) 
            res= pd.concat([res, temp_res],axis=0)
        # return the results 
        return res

    def get_num_peptide_per_experiment(self)->pd.DataFrame:
        """return a table containing the number of peptides in every experiment in the current set.

        :return: A table with two columns, the first is the experiment name and the second is the peptide count per experiment.
        :rtype: pd.DataFrame
        """
        results:pd.DataFrame= pd.DataFrame({ exp_name:[len(exp.get_peptides())] for exp_name, exp in self._exps.items()},
                                index=['experiment','num_peptide']).T
        results.sort_values(by='num_peptide', inplace=True, ascending=False)
        return results

    def get_num_proteins_per_experiment(self)->pd.DataFrame:
        """return a table containing the number of proteins in every experiment in the current set.

        :return: A table with two columns, the first is the experiment name and the second is the protein count per experiment.
        :rtype: pd.DataFrame
        """
        results:pd.DataFrame= pd.DataFrame({ exp_name:[len(exp.get_proteins())] for exp_name, exp in self._exps.items()},
                                index=['experiment','num_proteins']).T
        results.sort_values(by='num_peptide', inplace=True, ascending=False)
        return results
    
   














