#!/usr/bin/env python 
"""This module provides an abstraction for an IP experiment. 
"""
# load the modules 
from __future__ import annotations
import time
import pandas as pd
import numpy as np 
from IPTK.Classes.Peptide import Peptide
from IPTK.Classes.Protein import Protein
from IPTK.Classes.Proband import Proband
from IPTK.Classes.HLASet import HLASet
from IPTK.Classes.Tissue import Tissue
from IPTK.Classes.Database import SeqDB, OrganismDB
from IPTK.Utils.Mapping import map_from_uniprot_gene
from IPTK.Utils.Types import Sequences, MappedProtein, MappedProteins,ProteinSource
from tqdm import tqdm 
from typing import List, Dict, Set, Tuple
# define the analysis types 
Peptides=List[Peptide]
Proteins=List[Protein]
# define the experiment class 
class Experiment: 
	"""A representation of an immunopeptidomics experiment. 
	"""
	def __init__(self, proband:Proband, hla_set:HLASet, tissue:Tissue,  database:SeqDB,
				ident_table:pd.DataFrame,progress_bar:bool=True)->Experiment: 
		"""Constructs an Experiment instance.

		:param proband: A proband instance containing the proband, name& other meta-data. 
		:type proband: Proband
		:param hla_set: an HLASet instance containing the set of alleles from which the peptide was eluted.
		:type hla_set: HLASet
		:param tissue: an instance of type tissue containing expression values and protein location from the corresponding tissue.
		:type tissue: Tissue
		:param database: a sequence database to exact the sequences of the identified proteins.
		:type database: SeqDB
		:param ident_table: The identification table which contain the peptides inferred from analyzing raw mass spec files.  
		:type ident_table: pd.DataFrame
		:raises ValueError: incase the provided input does not match its proposed type. 
		:return: an experiment instance 
		:rtype: Experiment
		"""
		# check the type of the inputs 
		if not isinstance(proband, Proband):
			raise ValueError(f'proband must be a Proband instance, however, your input is of type: {type(proband)}')
		
		if not isinstance(hla_set,HLASet):
			raise ValueError(f'hla_set must be an HLASet instance, however, your input is of type: {type(hla_set)}')

		if not isinstance(database, SeqDB):
			raise ValueError(f'database must be a DataBase instance, however, your input is of type: {type(database)}')

		if not isinstance(tissue, Tissue):
			raise ValueError(f'tissue must be a Tissue instance, however, your input is of type: {type(SeqDB)}')

		if not isinstance(ident_table, pd.DataFrame):
			raise ValueError(f'ident_table must be a pd.DataFrame instance, however, your input is of type: {type(ident_table)}')	
		
		if ident_table.shape[0]==0:
			raise ValueError('The identification tables can not be an empty table')

		if ident_table.shape[1] != 4:
			raise ValueError('The identification table must have 4 columns see IO modules for more details')

		# add the variables 
		self._proband: Proband =proband
		self._tissue: Tissue =tissue
		self._database: SeqDB =database
		self._hla_set: HLASet =hla_set
		# store the identified peptides as a dict object with peptide sequences as keys and Peptide instance as tables 
		self._peptides: Dict[str,Peptide] =dict()
		# store the identified proteins as set of protein ids 
		self._proteins: List[str]=list(set(ident_table.iloc[:,1].tolist()))
		# fill the instance dict 
		print(f"Creating an experimental object, ... started at: {time.ctime()}")
		for idx in tqdm(range(ident_table.shape[0])): 
			# if the peptide is not in the peptide dictionary we add to it. 
			if ident_table.iloc[idx,0] not in self._peptides.keys():
				self._peptides[ident_table.iloc[idx,0]]= Peptide(ident_table.iloc[idx,0])
			# link the parent protein to its child peptide. 
			self._peptides[ident_table.iloc[idx,0]].add_parent_protein(
				Protein(ident_table.iloc[idx,1], self._database.get_seq(ident_table.iloc[idx,1])), 
				ident_table.iloc[idx,2],ident_table.iloc[idx,3])
		# the instance sequence has been filled 
	
	def get_peptides_length(self)->List[int]: 
		"""return a list containing the length of each unique peptide in the database.

		:return: peptides' lengths 
		:rtype: List[int]
		"""
		peptide_length=[]
		for peptide in self._peptides.keys(): 
			peptide_length.append(len(self._peptides[peptide]))
		return peptide_length

	def annotate_proteins(self, organisms_db: OrganismDB)->None:
		"""Extract the parent organisms of each protein in the experiment from an organism database instance. 

		:param organisms_db: an OrgansimDB instance that will be used to annotate the proteins \
		identified in the experiment.
		:type organisms_db: OrganismDB
		"""
		print(f"annotating proteins with source organism information ..., started at: {time.ctime()}")
		# get the organism
		proteins: List[str]= self.get_proteins()
		# allocate a list to store the parent proteins source organisms
		orgs: Dict[str,str] = dict()
		#  fill the list of elements in the list 
		for protein in tqdm(proteins):
			orgs[protein]=organisms_db.get_org(protein) # This might rise a KeyError incase the protein is not in the database
		# add the organism info 
		self.add_org_info(orgs)
		return 	

	def get_orgs(self)->List[str]:
		"""return a list containing the UNIQUE organisms identified in the current experiment 

		:return: list of all UNIQUE organisms inferred from the inferred proteins.   
		:rtype: List[str]
		"""
		print(f"Getting source organism information ..., started at: {time.ctime()}")
		# allocate a list to hold the results
		unique_results: List[str]=[]
		# get the organisms
		for peptide in tqdm(self.get_peptides()):
			unique_results.extend(self._peptides[peptide].get_parents_org())
		# create a set the contain the unique proteins
		unique_results = list(set(unique_results))
		# return the results 
		return unique_results

	def get_peptides_per_organism(self)->pd.DataFrame:
		"""return a pandas dataframe that contain the count of peptides belonging to each organism in
		the database

		:return: a table with two columns, namely, Organisms and Counts
		:rtype: pd.DataFrame
		"""
		print(f"Computing peptides per organism table ..., starting at: {time.ctime()}")
		# allocate a list to hold the organims per proteins
		peptides_per_organims : Dict[str, int] = dict()
		# get the number of organims 
		organisms: List[str] = self.get_orgs()
		# initialize the counter with zeros 
		for org in organisms:
			 peptides_per_organims[org] = [0]
		# obtain the data from the peptides 
		for pep in tqdm(self.get_peptides()): 
			for org in self._peptides[pep].get_parents_org():
				peptides_per_organims[org][0]+=1 # increment the counter 
		# create a dataframe
		res: pd.DataFrame = pd.DataFrame(peptides_per_organims).T
		# add the index as an extra-columns 
		res['Organisms'] = res.index.tolist()
		res.columns=['Counts','Organisms']
		# reformat the dataframe 
		res.reset_index(drop=True,inplace=True)
		# swap the columns
		res=res.reindex(columns=['Organisms','Counts'])
		# sort the results 
		res=res.sort_values(by='Counts',ascending=False)
		return res
	
	def drop_peptide_belong_to_org(self, org:str)->None: 
		"""Drop the all the peptides that belong to a user provided organism. \
		Note that, this function will IRREVERSIBLY remove the peptide from the experimental object. 

		:param org: the organims name
		:type org: str
		"""
		print(f"Removing peptides belonging to: {org} ..., starting at: {time.ctime()}")
		# remove the peptides belonging to the provided organisms
		for pep in tqdm(self.get_peptides()):
			if org in self._peptides[pep].get_parents_org():
				self._peptides.pop(pep) # remove the peptide from the database 
		# update the list of proteins in the database by dropping proteins with no associated peptide
		print(f"Updating protein list..., starting at: {time.ctime()}") 
		new_proteins: List[str] = []
		for pep in tqdm(self.get_peptides()):
			new_proteins.extend(self._peptides[pep].get_parent_proteins())
		# update the list of parent proteins
		self._proteins=list(set(new_proteins))
		return 
	
	def get_experiment_reference_tissue_expression(self) ->pd.DataFrame:
		"""return the reference gene expression for the current tissue 
	
		:return: A table that contain the expression value for ALL the genes in the instance Tissue 
		:rtype: pd.DataFrame
		"""
		return self._tissue.get_expression_profile().get_table()

	def add_org_info(self, prot2org: ProteinSource)->None:
		"""annotated the inferred proteins with their source organism

		:param prot2org: a dict that contain the protein id as keys and its source organism as values \
		and add this info to each protein inferred in the current experiment.
		:type prot2org: ProteinSource
		:raises RuntimeWarning: If the provided dictionary does cover all proteins in the experimental object.
		"""
		if len(prot2org) < len(self._proteins):
			raise RuntimeWarning(f'The provided dictionary does cover all proteins in the experimental object')
		print(f"Adding organism information ..., starting at: {time.ctime()}")
		# add the info to each parent protein
		for prot in tqdm(prot2org.keys()): 
			for pep in self._peptides.keys():
				if self._peptides[pep].is_child_of(prot):
					self._peptides[pep].add_org_2_parent(prot, prot2org[prot])
		return
	
	def get_flanked_peptides(self,flank_length: int) -> Sequences:
		"""returns a list of sequences containing the peptides identified in the experiment padded with
		the flanking regions from all the parents of each peptide.

		:param flank_length: the length of the flanking region 
		:type flank_length: int
		:return: a list of the peptides + the flanking region. 
		:rtype: Sequences
		"""
		results: List[str]=[]
		for pep_idx in tqdm(self._peptides.keys()): 
			results.extend(self._peptides[pep_idx].get_flanked_peptide(flank_length))
		return results
	
	def get_negative_example(self, fold: int= 2)->Sequences: 
		"""generate negative examples, i.e., non-bounding peptides from the proteins identified in the current experiment.  

		:param fold: the number of negative example to generate relative to the number of unique identified peptides, defaults to 2
		:type fold: int, optional
		:return: list of non-presented peptides from all inferred proteins. 
		:rtype: Sequences
		"""
		results: List[str]=[]
		for _ in tqdm(range(fold)): 
			for pep_idx in self._peptides.keys(): 
				results.extend(self._peptides[pep_idx].get_non_presented_peptides(len(pep_idx)))
		return results
 
	def get_binarized_results(self)->MappedProtein: 
		"""Return a list of NumPy arrays where each array represents a child peptide, parent protein mapped pair.
		Please note that, The function treat each peptide-protein pair individually, that is if two peptides originating from the same protein, 
		it treat them independently and the same protein will be represented twice with the two different peptides. Incase an integrative mapping is needed,
		the function @get_integrated_binarized_results@ shall be used.
		
		:return: a list of NumPy arrays containing the mapping between each peptide protein pair. 
		:rtype: MappedProtein
		"""
		results=[]
		for pep_idx in tqdm(self._peptides.keys()): 
			results.extend(self._peptides[pep_idx].map_to_parent_protein())
		return results 
	
	def get_peptide(self, pep_seq: str) -> Peptide:
		"""return a peptide instance corresponding to the user provided peptide sequence.
	
		:param pep_seq: the peptide sequence 
		:type pep_seq: str
		:raises KeyError: if the peptide sequence has not been inferred from the current database. 
		:return: the peptide instance with the coresponding sequence
		:rtype: Peptide
		"""
		if pep_seq not in tqdm(self._peptides.keys()): 
			raise KeyError(f"your provided petide sequence: {pep_seq} is not defined in the current experiment.") 
		return self._peptides[pep_seq]
	
	def get_mapped_protein(self, pro_id: str)->np.ndarray: 
		"""return an NumPy array of shape 1 x protein length where each number in the array represents 
		the total number of peptides identified in the experiment that have originated from the said position 
		in the protein.
		
		:param pro_id: the protein id
		:type pro_id: str
		:raises KeyError: if the provided protein id was inferred from the current experiment
		:return: a NumPy array that contain the coverage of the protein.
		:rtype: np.ndarray
		"""
		# check if the protein have been identified in the experiment 
		if pro_id not in self._proteins:
			raise KeyError(f"The provided id: {pro_id} is not in the set of identified proteins")
		# store the index of the start and end position of the child peptides   
		start_idxs, end_idxs= [], []
		# obtain all the peptides and a protein instance          
		for pep in self._peptides.keys():
			if self._peptides[pep].is_child_of(pro_id):
				protein=self._peptides[pep].get_parent(pro_id)
				start,end=self._peptides[pep].get_pos_in_parent(pro_id)
				start_idxs.append(start)
				end_idxs.append(end)
		# return the results 
		return protein.get_peptides_map(start_idxs,end_idxs)
		
	def get_expression_of_parent_proteins(self, non_mapped_dval: float = -1)->pd.DataFrame:
		"""return a table containing the expression of the proteins inferred in the current experiment from the current tissue.
		This method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  

		:param non_mapped_dval: A default value to be added incase the parent protein is not define in the expression database, defaults to -1
		:type non_mapped_dval: float, optional
		:return: a table that contain the expression of the protein inferred in the database 
		:rtype: pd.DataFrame
		"""
		proteins: List[str] = list(self.get_proteins())
		print(f"Mapping Uniprot accession to ENSEMBLE IDs ..., starting at: {time.ctime()}")
		map2Ensemble: pd.DataFrame = map_from_uniprot_gene(proteins)
		# allocate a list to hold the expression values 
		expression: List[float] = []
		print(f"Computing the expression of parent proteins ..., starting at: {time.ctime()}")
		for prot in tqdm(proteins):
			# get  all transcripts that map to the protein -> port
			temp_df: pd.DataFrame = map2Ensemble.loc[map2Ensemble.iloc[:,0]==prot]
			# if one-to-one mapping is returned 
			if temp_df.shape[0] ==1:
				# we try to extract the expression value of the protein
				try:
					expression.append(self._tissue.get_expression_profile().get_gene_id_expression(temp_df.iloc[0,1]))
				except KeyError:
					expression.append(non_mapped_dval)
			else: # we have more than one mapping 
				temp_ens_ids: List[str] =  temp_df.iloc[:,1].tolist()
				temp_res_raw: List[int] = []
				for ens_id in temp_ens_ids:
					try: 
						temp_res_raw.append(self._tissue.get_expression_profile().get_gene_id_expression(ens_id))
					except KeyError:
						temp_res_raw.append(non_mapped_dval)
				# filter out default value 
				temp_res_pross: List[int] = [elem for elem in temp_res_raw if elem != non_mapped_dval]
				# if the list is empty, all the transcript can not be mapped 
				if len(temp_res_pross)==0:
					expression.append(non_mapped_dval)
				else: 
					expression.append(np.mean(temp_res_pross))	 
		# construct the dataframe 
		results: pd.DataFrame= pd.DataFrame({'Proteins':proteins, 'Expression':expression})
		return results
	
	def get_main_sub_cellular_location_of_parent_proteins(self, not_mapped_val: str = 'UNK')->pd.DataFrame:
		"""retrun the main cellular location for the identified proteins.
		This method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  


		:param not_mapped_val: The default value to return incase the location of a protein can not be extracted, defaults to 'UNK'
		:type not_mapped_val: str, optional
		:return: A table that contain the main cellular compartment for each protein in the current instance.
		:rtype: pd.DataFrame
		"""
		proteins: List[str] = list(self.get_proteins())
		print(f"Mapping Uniprot accession to ENSEMBLE IDs ..., starting at: {time.ctime()}")
		map2Ensemble: pd.DataFrame = map_from_uniprot_gene(proteins)
		# allocate a list to hold the main location
		main_locations: List[str] = []
		print(f"Getting the Subcellular location of parent proteins ..., starting at: {time.ctime()}")
		for prot in tqdm(proteins):
			# we get a pandas dataframe that contain all the ensemble ids belonging to this protein.  
			temp_df: pd.DataFrame = map2Ensemble.loc[map2Ensemble.iloc[:,0]==prot]
			if temp_df.shape[0]==1: 
				try:
					main_locations.append(';'.join(self._tissue.get_subCellular_locations().get_main_location(temp_df.iloc[0,1])))
				except KeyError: 
					main_locations.append(not_mapped_val)
			else: 
				temp_ens_ids: List[str] = temp_df.iloc[:,1].tolist()
				temp_res_raw: List[str] = []
				for ens_id in temp_ens_ids:
					try: 
						temp_res_raw.append(';'.join(self._tissue.get_subCellular_locations().get_main_location(ens_id)))
					except KeyError: 
						temp_res_raw.append(not_mapped_val)
				# filter out default value 
				temp_res_pross: List[int] = [elem for elem in temp_res_raw if elem != not_mapped_val]
				# if the list is empty, all the proteins can not be mapped 
				if len(temp_res_pross)==0:
					main_locations.append(not_mapped_val)
				else: 
					# get a set of the unique location from different mapping 
					temp_unique_poss: Set[str] = set()
					for elem in temp_res_pross: 
						for loc in elem.split(';'):
							temp_unique_poss.add(loc)
					# append the results into one string and add it to the database elements 		
					main_locations.append(';'.join(temp_unique_poss))
		# construct the dataframe 
		results: pd.DataFrame= pd.DataFrame({'Proteins':proteins, 'Main_locations':main_locations})
		return results

	
	def get_go_location_id_parent_proteins(self, not_mapped_val: str = 'UNK')->pd.DataFrame:
		"""retrun the gene ontology,GO, location terms for all the identified proteins. 

		:param not_mapped_val: The default value to return incase the GO term of the protein can not be extracted, defaults to 'UNK'
		:type not_mapped_val: str, optional
		:return: A table that contain the GO-location term for each protein in the current instance.
		:rtype: pd.DataFrame
		"""
		proteins: List[str] = list(self.get_proteins())
		print(f"Mapping Uniprot accession to ENSEMBLE IDs ..., starting at: {time.ctime()}")
		map2Ensemble: pd.DataFrame = map_from_uniprot_gene(proteins)
		print(f"Getting the GO Subcellular compartment of parent proteins ..., starting at: {time.ctime()}")
		#allocate a list to hold the go terms
		go_terms: List[str] = []
		for prot in tqdm(proteins):
			# we get a pandas dataframe that contain all the ensemble ids belonging to this protein.  
			temp_df: pd.DataFrame = map2Ensemble.loc[map2Ensemble.iloc[:,0]==prot]
			if temp_df.shape[0]==1: 
				try:
					go_terms.append(';'.join(self._tissue.get_subCellular_locations().get_go_names(temp_df.iloc[0,1])))
				except KeyError: 
					go_terms.append(not_mapped_val)
			else: 
				temp_ens_ids: List[str] = temp_df.iloc[:,1].tolist()
				temp_res_raw: List[str] = []
				for ens_id in temp_ens_ids:
					try: 
						temp_res_raw.append(';'.join(self._tissue.get_subCellular_locations().get_go_names(ens_id)))
					except KeyError: 
						temp_res_raw.append(not_mapped_val)
				# filter out default value 
				temp_res_pross: List[int] = [elem for elem in temp_res_raw if elem != not_mapped_val]
				# if the list is empty, all the proteins can not be mapped 
				if len(temp_res_pross)==0:
					go_terms.append(not_mapped_val)
				else: 
					# get a set of the unique location from different mapping 
					temp_unique_poss: Set[str] = set()
					for elem in temp_res_pross: 
						for loc in elem.split(';'):
							temp_unique_poss.add(loc)
					# append the results into one string and add it to the database elements 		
					go_terms.append(';'.join(temp_unique_poss))
		# construct the dataframe 
		results: pd.DataFrame= pd.DataFrame({'Proteins':proteins, 'GO_Terms':go_terms})
		return results

	def get_num_peptide_expression_table(self)->pd.DataFrame:
		"""Get a table that contain the id of all parent proteins, number of peptide per-proteins and the expression value 
		   of these parent transcripts. Please note, this method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  

		:return: the number of peptides per protein table
		:rtype: pd.DataFrame
		"""
		# get the number of tables per peptides 
		num_peptides_per_protein: pd.DataFrame = self.get_peptides_per_protein()
		expression_level: pd.DataFrame = self.get_expression_of_parent_proteins()
		# merge the tables 
		results: pd.DataFrame =pd.merge(num_peptides_per_protein,expression_level)
		# return the results 
		return results

	def get_number_of_proteins_per_compartment(self) -> pd.DataFrame: 
		"""returns the number of proteins from each compartment 
		
		:return: A table that has two columns, namely, Compartment and Counts. 
		:rtype: pd.DataFrame
		"""
		# get the main locations: 
		parent_protein_locations: pd.DataFrame = self.get_main_sub_cellular_location_of_parent_proteins()
		# obtain the locations as a list 
		locations: List[str] = []
		for loc in parent_protein_locations.iloc[:,1].tolist(): 
			locations.extend(loc.split(';'))
		# construct a dictionary to hold the counts 
		unique_compartments: Set[str] = set(locations)
		compartment_counts: Dict[str,int] = dict()
		# initialize the counts 
		for comp in unique_compartments: 
			compartment_counts[comp]=[0]
		# update the counter 
		print(f"counting proteins per compartment ..., starting at: {time.ctime()}")
		for loc in tqdm(locations): 
			compartment_counts[loc][0]+=1
		# construct a data frame from the results 
		res: pd.DataFrame = pd.DataFrame(compartment_counts).T
		# add the index as an extra-columns 
		res['Compartment'] = res.index.tolist()
		res.columns=['Counts','Compartment']
		# reformat the dataframe 
		res.reset_index(drop=True,inplace=True)
		# swap the columns
		res=res.reindex(columns=['Compartment','Counts'])
		# sort the results 
		res=res.sort_values(by='Counts',ascending=False)
		return res

	def get_number_of_proteins_per_go_term(self) -> pd.DataFrame: 
		"""returns the number of proteins from each GO-Term 
		
		:return: A table that has two columns, namely, GO-Terms and Counts. 
		:rtype: pd.DataFrame
		"""
		# get the go-terms
		parent_protein_go_term : pd.DataFrame = self.get_go_location_id_parent_proteins()
		# obtain the locations as a list 
		go_terms: List[str] = []
		for terms in parent_protein_go_term.iloc[:,1].tolist(): 
			go_terms.extend(terms.split(';'))
		# construct a dictionary to hold the counts 
		unique_terms: Set[str] = set(go_terms)
		terms_counts: Dict[str,int] = dict()
		# initialize the counts 
		for term in unique_terms: 
			terms_counts[term]=0
		# update the counter 
		print(f"counting proteins per GO term ..., starting at: {time.ctime()}")
		for term in go_terms: 
			terms_counts[term]+=1
		# prepare the dict to be compatible with dataframe 
		for key in  terms_counts.keys():
			terms_counts[key]= [terms_counts[key]]
		# construct a data frame from the results 
		res: pd.DataFrame= pd.DataFrame(terms_counts).T
		# add the terms as a column
		res['Terms']= res.index.tolist()
		res.columns=['Counts','GO-Terms']
		# reformat the dataframe 
		res.reset_index(drop=True,inplace=True)
		# swap the columns
		res=res.reindex(columns=['GO-Terms','Counts'])
		# sort the results 
		res=res.sort_values(by='Counts',ascending=False)
		return res

	def get_num_peptide_per_location(self)->pd.DataFrame:
		"""retruns the number of peptides obtained from proteins localized to different sub-cellular compartments  
		
		:return: A table that has two columns, namely, Compartment and Counts.
		:rtype: pd.DataFrame
		"""
		# get unique compartments 
		unique_locations: List[str]= self.get_number_of_proteins_per_compartment().iloc[:,0].tolist()
		# initialize the counter to zeros
		pep_per_loc: Dict[str,int]=dict()
		for loc in unique_locations:
			pep_per_loc[loc]=0
		# get the location of each parent protein 
		parent_protein_locs: pd.DataFrame = self.get_main_sub_cellular_location_of_parent_proteins()
		peptide_count_parents: pd.DataFrame = self.get_peptides_per_protein()
		# loop over the parent proteins and update the countert 
		print(f"Computing results table, ... starting at: {time.ctime()}")
		for idx in tqdm(range(parent_protein_locs.shape[0])):
			# get the number of peptides belonging to this protein  
			num_peptides: int = peptide_count_parents.loc[peptide_count_parents.iloc[:,0]==parent_protein_locs.iloc[idx,0]]['Number_of_Peptides'].tolist()[0]
			# get the locations 
			locations: List[str] = parent_protein_locs.iloc[idx,1].split(';')
			# add the locations to the list 
			for loc in locations:
				pep_per_loc[loc]+=num_peptides
		# make the dict compatible with dataframes 
		for key in pep_per_loc.keys():
			pep_per_loc[key]=[pep_per_loc[key]]
		# construct a data frame from the results 
		res: pd.DataFrame = pd.DataFrame(pep_per_loc).T
		# add the index as a columns 
		res['Compartment']=res.index.tolist()
		res.columns=['Counts','Compartment']
		# reformat the dataframe 
		res.reset_index(drop=True,inplace=True)
		# swap the columns
		res=res.reindex(columns=['Compartment','Counts'])
		# sort the results 
		res=res.sort_values(by='Counts',ascending=False)
		# return the results 
		return res 

	def get_num_peptide_per_go_term(self)->pd.DataFrame:
		"""retruns the number of peptides per each GO-Term 
		:return: A table that has two columns, namely, GO-Terms and Counts. 
		:rtype: pd.DataFrame
		"""
		# get GO terms 
		unique_go_terms: List[str]= self.get_number_of_proteins_per_go_term().iloc[:,0].tolist()
		# initialize the counter to zeros
		pep_per_term: Dict[str,int]=dict()  
		for loc in unique_go_terms:
			pep_per_term[loc]=0
		# get the Go-Terms of each parent protein 
		parent_protein_go_terms: pd.DataFrame = self.get_go_location_id_parent_proteins()
		peptide_count_parents: pd.DataFrame = self.get_peptides_per_protein()
		# loop over the parent proteins and update the countert 
		print(f"Computing results table, ... starting at: {time.ctime()}")
		for idx in tqdm(range(parent_protein_go_terms.shape[0])):
			# get the number of peptides belonging to this protein  
			num_peptides: int = peptide_count_parents.loc[peptide_count_parents.iloc[:,0]==parent_protein_go_terms.iloc[idx,0]]['Number_of_Peptides'].tolist()[0]
			# get the locations 
			go_terms: List[str] = parent_protein_go_terms.iloc[idx,1].split(';')
			# add the locations to the list 
			for term in go_terms:
				pep_per_term[term]+=num_peptides
		# make the dict compatible with dataframes 
		for key in pep_per_term.keys():
			pep_per_term[key]=[pep_per_term[key]]
		# construct a data frame from the results 
		res: pd.DataFrame = pd.DataFrame(pep_per_term).T
		# add the index as a columns 
		res['Compartment']=res.index.tolist()
		res.columns=['Counts','GO-Terms']
		# reformat the dataframe 
		res.reset_index(drop=True,inplace=True)
		# swap the columns
		res=res.reindex(columns=['GO-Terms','Counts'])
		# sort the results 
		res=res.sort_values(by='Counts',ascending=False)
		# return the results 
		return res 

	def get_mapped_proteins(self)->MappedProteins: 
		"""returns a dictionary of all the proteins identified in the current experiment with all inferred
		peptides mapped to them. 

		:return: a dictionary that contain the mapped proteins for all the proteins in the current instance. 
		:rtype: MappedProteins
		"""
		results=dict()
		for prot in tqdm(self._proteins):
			results[prot]=self.get_mapped_protein(prot)
		return results 
	
	def get_mono_parent_peptides(self)->Peptides:
		"""
		returns a list of peptides that have only one parent protein
		
		:return: list of peptide instance 
		:rtype: Peptides
		"""
		results: List[Peptide]=[]
		for pep in tqdm(self._peptides.keys()):
			if self._peptides[pep].get_number_of_parents()==1:
				results.append(self._peptides[pep])
		return results
	
	def get_poly_parental_peptides(self)->Peptides:	
		"""returns a list of peptides that have more than one parent protein
		:return: [list of peptide instance 
		:rtype: Peptides
		"""
		results: List[Peptide]=[]
		for pep in tqdm(self._peptides.keys()):
			if self._peptides[pep].get_number_of_parents()>1:
				results.append(self._peptides[pep])
		return results
	
	def get_peptide_number_parent(self, ascending: bool = False)-> pd.DataFrame:
		"""returns a pandas dataframe with the peptide sequence in the first columns and the 
		number of parent proteins in the second column. 
		
		:param ascending: ascending sort the peptide by their number of parent proteins, defaults to False
		:type ascending: bool, optional
		:return: the number of parents for each peptide 
		:rtype: pd.DataFrame
		"""
		peptides=list(self._peptides.keys())
		num_parents=[]
		for second in tqdm(peptides): 
			num_parents.append(self._peptides[second].get_number_of_parents())
		# combine the results into a dataframe 
		res=pd.DataFrame({'Peptides':peptides,'Number_of_parents':num_parents})
		# sort the results 
		res.sort_values(axis=0, by='Number_of_parents', ascending=ascending, inplace=True)
		# return the results 
		return res 
	
	def get_number_of_children(self,pro_id: str)->int:
		"""returns the number of children, i.e. number of peptides belonging to a parent protein

		:param pro_id: the id of the parent protein 
		:type pro_id: str
		:return: the number of peptides
		:rtype: int
		"""
		children=0
		for pep in self._peptides.keys(): 
			if self._peptides[pep].is_child_of(pro_id):
				children+=1
		return children
	
	def get_peptides_per_protein(self, ascending: bool = False)->pd.DataFrame:
		"""returns a pandas dataframe that contain the number of peptides belonging to each protein 
		inferred in the experiment
		
		:param ascending: ascending sort the proteins by their number of parents each child peptide has, defaults to False
		:type ascending: bool, optional
		:return: A table with the following columns, Proteins and Number_of_Peptides
		:rtype: pd.DataFrame
		"""
		proteins=list(self._proteins)
		number_peptides=[]
		for prot in tqdm(proteins): 
			number_peptides.append(self.get_number_of_children(prot))
		# construct a pandas dataframe with the proteins and number of peptides 
		res=pd.DataFrame({'Proteins':proteins,'Number_of_Peptides':number_peptides})
		# sort the resulting table 
		res.sort_values(axis=0, by='Number_of_Peptides',ascending=ascending, inplace=True)
		# return the results 
		return res 

	def get_n_terminal_flanked_seqs(self, flank_length: int)->Peptides:
		"""returns the n-terminal flanking sequences 

		:param flank_length: the length of the flanking region upstream of the N-terminal of the peptide. 
		:type flank_length: int
		:return: a list of sequences containing the N-terminal flanking sequence for each peptide in the instance.  
		:rtype: Peptides
		"""
		results=[]
		for pep_id in tqdm(self._peptides.keys()): 
			results.extend(self._peptides[pep_id].get_n_terminal_flank_seq(flank_length))
		return results 

	def get_c_terminal_flanked_seqs(self, flank_length: int)-> Peptides:
		"""returns the c-terminal flanking sequences

		:param flank_length: the length of the peptide downstream of the C-terminal of the peptide  
		:type flank_length: int
		:return: a list of sequences containing the N-terminal flanking sequence for each peptide in the instance.  
		:rtype: Peptides
		"""
		results=[]
		for pep_id in tqdm(self._peptides.keys()): 
			results.extend(self._peptides[pep_id].get_c_terminal_flank_seq(flank_length))
		return results 

	def get_tissue_name(self)->str:
		"""
		:returns: the name of the tissue 
		:rtype: str
		"""
		return self._tissue.get_name()
		
	def get_proband_name(self)->str:
		"""
		:return: the proband's name
		:rtype: str
		"""
		return self._proband.get_name()
	
	def get_hla_class(self)->int:
		"""
		:return: the HLA class
		:rtype: int
		"""
		return self._hla_set.get_class()
	
	def get_hla_allele(self)-> List[str]:
		"""
		:return: the set of HLA alleles from which the instance peptides have been eluted 
		:rtype: List[str]
		"""
		return self._hla_set.get_alleles()
	
	def get_hla_set(self)->List[str]: 
		"""return the instance hla set 

		Returns:
			HLASet: the instance hla set 
		"""
		return self._hla_set

	def has_hla_allele(self, individual: str)->bool:
		"""returns whether or not the experiment contain an eluted peptides from the provided alleles 
		
		:param individual: the name of the allele as a string
		:type individual: str
		:return: True if the allele is a member of the instance HLASet and False otherwise.  
		:rtype: bool
		"""
		return self._hla_set.has_allele(individual)
	
	def has_gene(self, locus: str)->bool:
		"""returns whether or not the experiment contains peptides eluted from an HLA-allele belonging to the provided locus or not
		
		:param locus: the locus of the allele to query the hla_set against
		:type locus: str
		:return: True if the locus has a member that is a member of the instance HLASet and False otherwise
		:rtype: bool
		"""
		return self._hla_set.has_gene(locus)
	
	def has_allele_group(self, gene_group:str)->bool:
		"""returns whether or not the experiment contains peptides eluted from an HLA-allele belonging to the provided allele group or not
		
		:param gene_group: the gene group to query the hla_set against 
		:type gene_group: str
		:return: True if the gene group has a member that is a member of the instance HLASet and False otherwise
		:rtype: bool
		"""
		return self._hla_set.has_allele_group(gene_group)
	
	def has_protein_group(self, protein_group: str)->bool:
		""" returns whether or not the experiment contains peptides eluted from an HLA-allele belonging to the provided protein group or not

		:param protein_group: The protein group to query the hla_set against 
		:type protein_group: str
		:return: True if the locus has a member that is a member of the instance HLASet and False otherwise
		:rtype: bool
		"""
		return self._hla_set.has_protein_group(protein_group)
	
	def get_peptides(self)->Peptides:
		"""
		:return: a set of all the peptides stored in the experimental object
		:rtype: Peptides
		"""
		return set(self._peptides.keys())

	def get_proteins(self)->Proteins:
		"""
		:return: a set of all the proteins in the experimental object
		:rtype: Proteins
		"""
		return set(self._proteins)
	
	def is_member(self,peptide:str)->bool:
		"""
		:param peptide: check if the peptide is a member of the instance peptides or not.
		:type peptide: str
		:return: True if the peptide has been identified in the current instance, False otherwise. 
		:rtype: bool
		"""
		if peptide in self._peptides.keys():
			return True
		return False
	
	def is_a_parent_protein(self, protein:str)->bool:
		"""
		:param protein: check if the protein is a member of the instance proteins or not.
		:type peptide: str
		:return: True if the protein has been identified in the current instance, False otherwise. 
		:rtype: bool
		"""
		if protein in self._proteins:
			return True
		return False 
	
	def get_tissue(self)->Tissue: 
		"""
		:return: the tissue of the current experiment.
		:rtype: Tissue
		"""
		return self._tissue

	def	__len__(self)->int: 
		""" A magic function for the len function, return the number of unique peptides in the experiment.

		:return: the number of unique peptides in the database
		:rtype: int
		"""
		return len(self._peptides)
	
	def __str__(self)->str:
		"""
		:return: a string representation for the tissue 
		:rtype: str
		"""
		return f"""an experiment from proband: {self._proband.get_name()}, tissue: {self._tissue.get_name()},
				   with an HLA Class: {self._hla_set.get_class()}. The instances contains 
				   {len(self)} peptides identified from {len(self._proteins)} proteins."""
	
	
	def __repr__(self)->str:
		return str(self)
	


				













		