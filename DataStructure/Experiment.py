#!/usr/bin/env python 
"""
@brief: an Experiment class handler 
@author: Hesham ElAbd
@version: 0.0.1
"""
# load the modules 
from __future__ import annotations
import pandas as pd
import numpy as np 
from IPTK.DataStructure.Peptide import Peptide
from IPTK.DataStructure.Protein import Protein
from IPTK.DataStructure.Proband import Proband
from IPTK.DataStructure.HLASet import HLASet
from IPTK.DataStructure.Tissue import Tissue
from IPTK.DataStructure.Database import SeqDB
from IPTK.Utils.Mapping import map_from_uniprot_gene
from IPTK.Utils.Types import Sequences, MappedProtein, MappedProteins,ProteinSource
from typing import List, Dict, Set, Tuple 
# define the analysis types 
Peptides=List[Peptide]
Proteins=List[Protein]
# define the experiment class 
class Experiment: 
	"""
	@brief A representation of an immunopeptidomic experiment. 
	"""
	def __init__(self, proband:Proband, hla_set:HLASet, tissue:Tissue,  database:SeqDB,
				ident_table:pd.DataFrame)->Experiment: 
		"""
		@brief: Construct an experiment instance. 
		@param: Proband_name:a proband instance that contain the proband, name& other  meta-data . 
		@param: Tissue: an instance of type tissue that store expression values for the corresponding tissue, @see tissue for more details
		@param: database: database to exact the sequence of the identified proteins. 
		@param: ident_table: The identification table @see IO.InFunctions for more details. 
		@see: IPTK.IO.InFunctions.load_identification_table for more details. 
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
		for idx in range(ident_table.shape[0]): 
			# if the peptide is not in the peptide dictionary we add to it. 
			if ident_table.iloc[idx,0] not in self._peptides.keys():
				self._peptides[ident_table.iloc[idx,0]]= Peptide(ident_table.iloc[idx,0])
			# link the parent protein to its child peptide. 
			self._peptides[ident_table.iloc[idx,0]].add_parent_protein(
				Protein(ident_table.iloc[idx,1], self._database.get_seq(ident_table.iloc[idx,1])), 
				ident_table.iloc[idx,2],ident_table.iloc[idx,3])
		# the instance sequence has been filled 
	
	def get_peptides_length(self)->List[int]: 
		"""
		@brief: return a list containing the length of each unique peptide in the database 
		"""
		return [len(pep) for pep in self.get_peptides()]


	def add_org_info(self, prot2org: ProteinSource)->None:
		"""
		@brief: annotated the inferred proteins with their source organism
		@param: a dict that contain the protein id as keys and its source organism as values 
		and add this info to each protein inferred in the current experiment.  
		"""
		if len(prot2org) < len(self._proteins):
			raise RuntimeWarning(f'The provided dictionary does cover all proteins in the experimental object')
		# add the info to each parent protein
		for prot in prot2org.keys(): 
			for pep in self._peptides.keys():
				if self._peptides[pep].is_child_of(prot):
					self._peptides[pep].add_org_2_parent(prot, prot2org[prot])
		return
	
	def get_flanked_peptides(self,flank_length: int) -> Sequences:
		"""
		@brief: returns a list of sequences containing the peptides identified in the experiment padded with
		the flanking regions from all the parents of each peptide.
		@param: flank_length: the length of the flanking region 
		""" 
		results: List[str]=[]
		for pep_idx in self._peptides.keys(): 
			results.extend(self._peptides[pep_idx].get_flanked_peptide(flank_length))
		return results
	
	def get_negative_example(self, fold: int= 2)->Sequences: 
		"""
		@brief: generate negative examples, i.e., non-bounding peptides from the proteins identified in the current experiment.  
		@param: fold: the number of negative example to generate relative to the number of identified peptides. Default is 2 
		"""
		results: List[str]=[]
		for _ in range(fold): 
			for pep_idx in self._peptides.keys(): 
				results.extend(self._peptides[pep_idx].get_non_presented_peptides(len(pep_idx)))
		return results
 
	def get_binarized_results(self)->MappedProtein: 
		"""
		@brief: Return a list of numpy array where each array represents a child peptide, parent protein mapped pair.
		@note The function treat each peptide-protein pair individually, that is if two peptides originating from the same protein, 
		it treat them independently and the same protein will be represented twice with the two different peptides. Incase an integrative mapping is need,
		the function @get_integrated_binarized_results@ shall be used.
		"""
		results=[]
		for pep_idx in self._peptides.keys(): 
			results.extend(self._peptides[pep_idx].map_to_parent_protein())
		return results 
	
	def get_peptide(self, pep_seq: str) -> Peptide:
		"""
		@brief: return a peptide instance corresponding to the user provided peptide sequence. 
		@param: pep_seq: the peptide sequence 
		"""
		if pep_seq not in self._peptides.keys(): 
			raise KeyError(f"your provided petide sequence: {pep_seq} is not defined in the current experiment.") 
		return self._peptides[pep_seq]
	
	def get_mapped_protein(self, pro_id: str)->np.ndarray: 
		"""
		@brief return an numpy array of shape 1 x protein length where each number in the array represents 
		the total number of peptides identified in the experiment that have originated from the said position 
		in the protein.
		@param pro_id: the protein id
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
		"""
		@brief: return a table containing the expression of the protein identified in the current experiment in the provided tissue.
		@param: non_mapped_dval: A default value to be added incase the parent protein is not define in the expression database. Default is -1
		@note: This method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  
		"""
		proteins: List[str] = list(self.get_proteins())
		map2Ensemble: pd.DataFrame = map_from_uniprot_gene(proteins)
		# allocate a list to hold the expression values 
		expression: List[float] = []
		for prot in proteins:
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
		"""
		@brief: retrun a the main cellular location for the identified proteins.
		@param: not_mapped_val: the default value to return incase the location of the protein can not be extracted. 
		@note: This method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  
		"""
		proteins: List[str] = list(self.get_proteins())
		map2Ensemble: pd.DataFrame = map_from_uniprot_gene(proteins)
		# allocate a list to hold the main location
		main_locations: List[str] = []
		for prot in proteins:
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
		"""
		@brief: retrun a the gene ontology,GO, location terms for the identified proteins. 
		@param: not_mapped_val: the default value to return incase the GO term of the protein can not be extracted. 
		@note: This method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  
		"""
		proteins: List[str] = list(self.get_proteins())
		map2Ensemble: pd.DataFrame = map_from_uniprot_gene(proteins)
		#allocate a list to hold the go terms
		go_terms: List[str] = []
		for prot in proteins:
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
		"""
		@brief: Get a table that contain the id of all parent proteins, number of peptide per-proteins and the expression value 
		of these parent transcripts. 
		@note: This method need internet connection as it need to access uniprot mapping API to map uniprot IDs to gene IDs.  
		"""
		# get the number of tables per peptides 
		num_peptides_per_protein: pd.DataFrame = self.get_peptides_per_protein()
		expression_level: pd.DataFrame = self.get_expression_of_parent_proteins()
		# merge the tables 
		results: pd.DataFrame =pd.merge(num_peptides_per_protein,expression_level)
		# return the results 
		return results

	def get_number_of_proteins_per_compartment(self) -> pd.DataFrame: 
		"""
		@brief: get number of proteins from each compartment 
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
			compartment_counts[comp]=0
		# update the counter 
		for loc in locations: 
			compartment_counts[loc]+=1
		# prepare the dict to be compatible with dataframe 
		for key in  compartment_counts.keys():
			compartment_counts[key]= [compartment_counts[key]]
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
		"""
		@brief: get the number of proteins per each GO term 
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
		"""
		@brief: retrun the number of peptides obtained from proteins localized to different sub-cellular compartments  
		""" 
		# get unique compartments 
		unique_locations: List[str]= self.get_number_of_proteins_per_go_term().iloc[:,0].tolist()
		# initialize the counter to zeros
		pep_per_loc: Dict[str,int]=dict()
		for loc in unique_locations:
			pep_per_loc[loc]=0
		# get the location of each parent protein 
		parent_protein_locs: pd.DataFrame = self.get_main_sub_cellular_location_of_parent_proteins()
		peptide_count_parents: pd.DataFrame = self.get_peptides_per_protein()
		# loop over the parent proteins and update the countert 
		for idx in range(parent_protein_locs.shape[0]):
			# get the number of peptides belonging to this protein  
			num_peptides: int = peptide_count_parents.loc[peptide_count_parents.iloc[:,0]==parent_protein_locs.iloc[idx,0],1]
			# get the locations 
			locations: List[str] = parent_protein_locs.iloc[idx,1].split(';')
			# add the locations to the list 
			for loc in locations:
				pep_per_loc[loc]+=num_peptides
		# construct a data frame from the results 
		res: pd.DataFrame = pd.DataFrame(pep_per_loc).T
		# add the index as a columns 
		res['Compartment']=res.index.tolist()
		# return the results 
		return res 

	def get_num_peptide_per_go_term(self)->pd.DataFrame:
		"""
		@brief: retrun the number of GO-Terms obtained from proteins localized to different sub-cellular compartments  
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
		for idx in range(parent_protein_go_terms.shape[0]):
			# get the number of peptides belonging to this protein  
			num_peptides: int = peptide_count_parents.loc[peptide_count_parents.iloc[:,0]==parent_protein_go_terms.iloc[idx,0],1]
			# get the locations 
			go_terms: List[str] = parent_protein_go_terms.iloc[idx,1].split(';')
			# add the locations to the list 
			for term in go_terms:
				pep_per_term[term]+=num_peptides
		# construct a data frame from the results 
		res: pd.DataFrame = pd.DataFrame(pep_per_term).T
		# add the index as a columns 
		res['Compartment']=res.index.tolist()
		# return the results 
		return res 

	def get_mapped_proteins(self)->MappedProteins: 
		"""
		@brief: return a dictionary of all the proteins identified in the current experiment with all inferred
		peptides mapped to them. 
		"""
		results=dict()
		for prot in self._proteins:
			results[prot]=self.get_mapped_protein(prot)
		return results 
	
	def get_mono_parent_peptides(self)->Peptides:
		"""
		@brief: return a list of peptides that have only one protein
		"""
		results: List[Peptide]=[]
		for pep in self._peptides.keys():
			if self._peptides[pep].get_number_of_parents()==1:
				results.append(self._peptides[pep])
		return results
	
	def get_poly_parental_peptides(self)->Peptides:	
		"""
		@brief: return a list of peptides that have more than one parent 
		"""
		results: List[Peptide]=[]
		for pep in self._peptides.keys():
			if self._peptides[pep].get_number_of_parents()>1:
				results.append(self._peptides[pep])
		return results
	
	def get_peptide_number_parent(self, ascending: bool = False)-> pd.DataFrame:
		"""
		@brief: return a pandas dataframe with the peptide sequence in the first columns and the 
		number of parent proteins in the second column. 
		@param: ascending sort the peptide by their number of parent proteins, default is False 
		"""
		peptides=list(self._peptides.keys())
		num_parents=[]
		for second in peptides: 
			num_parents.append(self._peptides[second].get_number_of_parents())
		# combine the results into a dataframe 
		res=pd.DataFrame({'Peptides':peptides,'Number_of_parents':num_parents})
		# sort the results 
		res.sort_values(axis=0, by='Number_of_parents', ascending=ascending, inplace=True)
		# return the results 
		return res 
	
	def get_number_of_children(self,pro_id: str):
		"""
		@brief: return the number of children belonging to a parent protein
		"""
		children=0
		for pep in self._peptides.keys(): 
			if self._peptides[pep].is_child_of(pro_id):
				children+=1
		return children
	
	def get_peptides_per_protein(self, ascending: bool = False)->pd.DataFrame:
		"""
		@brief: return a pandas dataframe that contain the number of peptides belonging to each protein 
		inferred in the experiment
		"""
		proteins=list(self._proteins)
		number_peptides=[]
		for prot in proteins: 
			number_peptides.append(self.get_number_of_children(prot))
		# construct a pandas dataframe with the proteins and number of peptides 
		res=pd.DataFrame({'Proteins':proteins,'Number_of_Peptides':number_peptides})
		# sort the resulting table 
		res.sort_values(axis=0, by='Number_of_Peptides',ascending=ascending, inplace=True)
		# return the results 
		return res 

	def get_n_terminal_flanked_seqs(self, flank_length: int)->Peptides:
		"""
		@brief: return the n-terminal flanking sequences 
		@param: flank_length: the length of the flanking region upstream of the N-terminal of the peptide 
		"""
		results=[]
		for pep_id in self._peptides.keys(): 
			results.extend(self._peptides[pep_id].get_n_terminal_flank_seq(flank_length))
		return results 

	def get_c_terminal_flanked_seqs(self, flank_length: int)-> Peptides:
		"""
		@brief: return the c-terminal flanking sequences∂ 
		@param: flank_length the length of the peptide upstream of the C-terminal of the peptide 
		"""
		results=[]
		for pep_id in self._peptides.keys(): 
			results.extend(self._peptides[pep_id].get_c_terminal_flank_seq(flank_length))
		return results 

	def get_tissue_name(self)->str:
		"""
		@brief: return the tissue name 
		"""
		return self._tissue.get_name()
		
	def get_proband_name(self)->str:
		"""
		@brief: return the proband name
		"""
		return self._proband.get_name()
	
	def get_hla_class(self)->int:
		"""
		@brief: return the HLA class
		"""
		return self._hla_set.get_class()
	
	def get_hla_allele(self)-> List[str]:
		"""
		@brief: return the allele set of the experiment proband
		"""
		return self._hla_set.get_alleles()
		
	def has_hla_allele(self, individual: str)->bool:
		"""
		@brief: return whether or not the experiment contain an eluted peptides from the provided alleles 
		@param: allele: is the name of the allele as a string
		"""
		return self._hla_set.has_allele(individual)
	
	def has_gene(self, locus: str)->bool:
		"""
		@brief: return whether or not the experiment contain peptides eluted from an HLA-alleles belonging to the provided locus or not
		@param: locus: the locus of the allele to query the hla_set against
		"""
		return self._hla_set.has_gene(locus)
	
	def has_allele_group(self, gene_group:str)->bool:
		"""
		@brief: return whether or not the experiment contain peptides eluted from an HLA-alleles belonging to the provided allele group or not
		@param: gene_group: the gene group to query the hla_set against 
		"""
		return self._hla_set.has_allele_group(gene_group)
	
	def has_protein_group(self, protein_group: str)->bool:
		"""
		@brief: return whether or not the experiment contain peptides eluted from an HLA-alleles belonging to the provided protein group or not
		@param protein_group: the protein group to query the hla_set against 
		"""
		return self._hla_set.has_protein_group(protein_group)
	
	def get_peptides(self)->Peptides:
		"""
		@brief: return a set of all the peptide stored in the experimental object
		"""
		return set(self._peptides.keys())
	
	def get_proteins(self)->Proteins:
		"""
		@brief: get a set of all the proteins in the experimental object
		"""
		return set(self._proteins)
	
	def is_member(self,peptide:str)->bool:
		"""
		@brief: check if a peptide is a member of the instance peptide or not.
		@param: peptide: a string containing the peptide sequence
		"""
		if peptide in self._peptides.keys():
			return True
		return False
	
	def is_a_parent_protein(self, protein:str)->bool:
		"""
		@brief: check if a protein is a member of the instances resources 
		@param: protein: the protein name or identifer
		"""
		if protein in self._proteins:
			return True
		return False 
	
	def	__len__(self)->int: 
		"""
		@brief: a magic function for the len function, return the number of peptides in the experiment.
		"""
		return len(self._peptides)
	
	def __str__(self)->str:
		"""
		@brief the implementation of the magic method str
		"""
		return f"""an Experimental from proband: {self._proband.get_name()}, Tissue: {self._tissue.get_name()}
				   With an HLA Class: {self._hla_set.get_class()} With
				   {len(self)} peptide identification from {len(self._proteins)} Protein"""
	
	
	def __repr__(self)->str:
		"""
		@brief the implementation of the magic method repr
		"""
		return str(self)
	


				













		