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
from IPTK.Utils.Types import Sequences, MappedProtein, MappedProteins,ProteinSource
from typing import List, Dict 
# define the analysis types 
Peptides=List[Peptide]
Proteins=List[Protein]
# define the experiment class 
class Experiment: 
	"""
	@brief A representation of an immunopeptidomic experiment. 
	"""
	def __init__(self, proband:Proband, hla_set:HLASet, tissue:Tissue,  database:SeqDB,
				ident_table:pd.DataFrame, sep: str ='/t')->Experiment: 
		"""
		@brief: Construct an experiment instance. 
		@param: Proband_name:a proband instance that contain the proband, name& other  meta-data . 
		@param: Tissue: an instance of type tissue that store expression values for the corresponding tissue, @see tissue for more details
		@param: database: database to exact the sequence of the identified proteins. 
		@param: ident_table: The identification table @see IO.InFunctions for more details. 
		@param: sep: The separator to parse the provided table. 
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
		results=[]
		for pep_idx in self._peptides.keys(): 
			results.extend(self._peptides[pep_idx].get_flanked_peptide(flank_length))
	
	def get_negative_example(self, fold: int= 2)->Sequences: 
		"""
		@brief: generate negative examples, i.e., non-bounding peptides from the proteins identified in the current experiment.  
		@param: fold: the number of negative example to generate relative to the number of identified peptides. Default is 2 
		"""
		results=[]
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
		results=[]
		for pep in self._peptides.keys():
			if self._peptides[pep].get_number_of_parents()==1:
				results.append(self._peptides[pep])
		return results
	
	def get_poly_parental_peptides(self)->Peptides:	
		"""
		@brief: return a list of peptides that have more than one parent 
		"""
		results=[]
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
	


				













		