#!/usr/bin/env python 
"""
@author: Hesham ElAbd 
@brief: A representation of the eluted peptides and its identified proteins.
@version: 0.0.2
# UPDATE IN VERSION 0.0.2 relative to 0.0.1
	I- adding the function add_org_2_parent which enable the peptide to add the source organism to one of 
	the instances parents 
	II- adding the function get_parents_org to get the list of organisms from which the parent was obtained
"""
# Load the modules
from __future__ import annotations
import numpy as np 
import pandas as pd
from IPTK.Utils.Types import Range, Sequences, MappedProtein, Organisms
# define the class 
class Peptide: 
	"""
	@brief: a representation of an eluted peptide. 
	"""
	def __init__(self,pep_seq: str)->Peptide:
		"""
		@brief: the class initializer, it initialize the class with a peptide sequence.
		@param: pep_seq: the peptide sequence 
		"""
		self._peptide=pep_seq
		self._parent_proteins=dict()

	def get_length(self) -> int:
		"""
		@brief: return the length of the peptides. 
		"""
		return len(self._peptide) 

	def get_peptide_seq(self) -> str: 
		"""
		@brief: return the sequence of the peptide.
		"""
		return self._peptide

	def add_parent_protein(self, parent_protein, 
						start_index: int, end_index: int) -> None:
		"""
		@brief: add a protein instance as a parent to the current peptide. 
		@param: parent_protein: a Protein instance that act as a parent to the peptide. 
		@param: start_index: the position in the parent protein where the peptide starts
		@param: the end_index: the index of the amino acid that occurs after the last amino acid in the peptide, 
		@note: The library use Python-based indexing where its 0-indexed and ranges are treated as [start, end). 
		"""
		if parent_protein[start_index:end_index] != self._peptide: 
			raise ValueError(f"""The provided parent information can not be attributed to the instance peptide
							 As according to the provided information. Extracted peptide is: 
							 { parent_protein[start_index:end_index]}
							 While the instance peptide is: 
							 {self._peptide}
							 """)
		self._parent_proteins[parent_protein.get_id()]={
		'protein':parent_protein, 
		'start_index':start_index, 
		'end_index':end_index}

	def get_number_parent_protein(self) -> int:
		"""
		@brief: returns the current number of parent proteins the peptide has. 
		"""
		return len(self._parent_proteins) 

	def get_flanked_peptide(self,flank_len:int)-> Sequences: 
		"""
		@brief: return the peptide sequences along with the flanking region from all the parent 
		@param:  flank_len: an integer represent the flanking region of the peptide. 
		@return: Sequences a list of string containing the length of the peptide + the flanking region from 
		both the N and C terminal of the instance peptide, from all proteins. 
		"""
		# if the parent has no parents, the function returns an empty list. 
		if self.get_number_parent_protein() == 0: 
			return [''] 
		# else: 
		# iterate over all the elements to obtain the sequences 
		results=[]
		for prot_id in self._parent_proteins.keys():
			parent_protein=self._parent_proteins[prot_id]
			protein_seq=parent_protein['protein']
			# check that the flaking region + protein length is within the protein length.
			if parent_protein['start_index'] - flank_len < 0: 
				start_idx=0
			else: 
				start_idx=parent_protein['start_index'] - flank_len
			# check the C-terminal side 
			if flank_len+parent_protein['end_index'] > len(protein_seq): 
				end_idx=len(protein_seq) 
			else: 
				end_idx= flank_len+parent_protein['end_index']	
			results.append(protein_seq[start_idx:end_idx])
		return results

	def map_to_parent_protein(self) ->MappedProtein:
		"""
		@brief mapped the peptide to its parent proteins and return a list of binary
		encoded arrays represent this mapping. 
		@details Mapped the instance peptide to the parent protein and returned a 
		list of numpy arrays where each array has a size of 1 by protein length. 
		within the protein the range representing the peptide is encoded as one while
		the rest is zero.  
		""" 
		# if the peptide has no associated parent protein yet, 
		# return an empty list 
		if self.get_number_parent_protein() == 0: 
			return [np.zeros(shape=(0,0))] 
		# else we allocate a list to compute the results 
		results=[]
		# compute the results for each parent protein
		for prot in self._parent_proteins.keys(): 
			pp= self._parent_proteins[prot]
			temp_array=np.zeros(shape=(1,len(pp['protein'])))
			temp_array[0,pp['start_index']:pp['end_index']]= 1
			results.append(temp_array) 
		return results 

	def get_non_presented_peptides(self, length: int )->Sequences: 
		"""
		@brief compute non-presented peptide from all the parent protein of the current peptide instance.  
		@param: length: The length, i.e. number of amino acids, for the non-presented peptide
		"""
		results=[]
		for pro in self._parent_proteins.keys():
			results.append(self._parent_proteins[pro]['protein'].get_non_presented_peptide(self._parent_proteins[pro]['start_index'],self._parent_proteins[pro]['end_index'], length=length))
		return results
	
	def get_parent_proteins(self):
 		"""
 		@brief return the set of all the parent of the instance peptide 
 		"""
 		return set(self._parent_proteins.keys())

 	
	def is_child_of(self,pro_id:str)->bool:
		"""
		@brief check if a user provided protein-id is a parent for the instance peptide or not.
		@param: pro_id: is the protein id 
		"""
		return pro_id in self.get_parent_proteins()
	
	def get_pos_in_parent(self,pro_id)->Range:
		"""
		@brief: Return the start and end position of the instance peptide in one of its parent
		proteins.
		@param: pro_id: is the protein-id 
		"""
		if not self.is_child_of(pro_id):
			raise ValueError(f"The provided protein id: {pro_id} is not a parent of the current peptide")
		# get the start and end elements 
		start=self._parent_proteins[pro_id]['start_index']
		end=self._parent_proteins[pro_id]['end_index']
		# return the results 
		return start, end 
	
	def get_parent(self,pro_id):
		"""
		@brief: return the parent protein that has an id matching the user defined pro_id
		@param: pro_id: The protein identifer 
		"""
		if self.is_child_of(pro_id):
			return self._parent_proteins[pro_id]['protein']

	def get_number_of_parents(self)->int:
		"""
		@brief: return the number of instance parent proteins
		"""	
		return len(self._parent_proteins)
	
	def get_n_terminal_flank_seq(self, flank_len : int):
		"""
		@brief return the flanking region from the N-terminal of the peptide in all of its parent proteins
		@param: flank_len: the length of the flanking regions 
		"""
		# if the peptide has no parent -> return a list with an empty string 
		if self.get_number_of_parents()==0:
			return ['']	
		#  else, we declare a list to allocate the results to it ==>  
		results=[]
		# loop over all the parent to obtain the upstream sequences 
		for parent in self._parent_proteins.keys(): 
			parent_protein=self._parent_proteins[parent]
			start_index=parent_protein['start_index']
			if start_index-flank_len<0: 
				results.append(parent_protein['protein'][0:start_index])
			else: 
				results.append(parent_protein['protein'][start_index-flank_len:start_index])
		# return the results 
		return results 
		 
	def get_c_terminal_flank_seq(self, flank_len: int):
		"""
		@brief return the flanking region from the C-terminal of the peptide in all of its parent proteins
		@param: flank_len: the length of the flanking region 
		"""
		# if the peptide has no parent -> return a list with an empty string 
		if self.get_number_of_parents()==0:
			return ['']	
		#  else, we declare a list to allocate the results to it ==>  
		results=[]
		# loop over all the parent to obtain the upstream sequences 
		for parent in self._parent_proteins.keys(): 
			parent_protein=self._parent_proteins[parent]
			end_index=parent_protein['end_index']
			protein_length=len(parent_protein['protein'])
			if end_index+flank_len>protein_length: # ==> avoid indexing outside of the parent protein
				results.append(parent_protein['protein'][end_index:protein_length])
			else: 
				results.append(parent_protein['protein'][end_index:end_index+flank_len])
		# return the results 
		return results 

	def add_org_2_parent(self,prot_name: str,org:str)->None:
		"""
		@brief: add the source organism of one of the instance parent protein
		@param: prot_name: The name of the protein 
		@param: org: the name of the organism 
		"""
		if not self.is_child_of(prot_name):
			raise ValueError(f'The provided protein id is not a parent of the current peptide')
		self._parent_proteins[prot_name]['protein'].set_org(org)
	
	def get_parents_org(self)->Organisms:
		"""
		@brief: return a list containing the name of each parent protein source organisms
		"""
		results=[]
		for protein in self._parent_proteins.keys():
			results.append(self._parent_proteins[protein]['protein'].get_org())
		return results 
	
	def __len__(self)->int:
		"""
		@brief: return the length of the peptide 
		"""
		return len(self._peptide)
	
	def __getitem__(self,aa_index:int)->str:
		"""
		@brief: support bracket based indexing into the peptide sequence 
		@param: aa_index: the amino acid index in the peptide sequence, for example 3rd amino acid
		"""
		return self._peptide[aa_index]
	
	def __str__(self)->str:
		"""
		@brief: return the sequence of the peptide 
		"""	
		return self._peptide

	def __repr__(self)->str:
		"""
		@brief: return a string representation for the class 
		"""
		return f'A Peptide instance with {self.get_number_of_parents()} parent protein'
	
		
		
		
		