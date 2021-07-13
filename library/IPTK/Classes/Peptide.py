#!/usr/bin/env python 
"""A representation of the eluted peptides and its identified proteins.
"""
# Load the modules
from __future__ import annotations
import numpy as np 
import pandas as pd
from IPTK.Classes.Protein import Protein
from typing import List 
from IPTK.Utils.Types import Range, Sequences, MappedProtein, Organisms
# define the class 
class Peptide: 
	"""An representation of an eluted peptide. 
	"""
	def __init__(self,pep_seq: str)->Peptide:
		"""class constructor 
		:param pep_seq: the peptide sequence 
		:type pep_seq: str
		:return: a Peptide instance 
		:rtype: Peptide
		"""
		self._peptide=pep_seq
		self._parent_proteins=dict()

	def get_length(self) -> int:
		"""
		:return: the length of the peptide
		:rtype: int
		"""
		return len(self._peptide) 

	def get_peptide_seq(self) -> str: 
		"""
		:return: the sequence of the peptide.
		:rtype: str
		"""
		return self._peptide

	def add_parent_protein(self, parent_protein:Protein, 
						start_index: int, end_index: int) -> None:
		"""add a protein instance as a parent to the current peptide. 
		The library use Python-based indexing where its 0-indexed and ranges are treated as [start, end). 
		:param parent_protein: a Protein instance that act as a parent to the peptide.
		:type parent_protein: Protein 
		:param start_index: the position in the parent protein where the peptide starts
		:type start_index: int 
		:param end_index: the index of the amino acid that occurs after the last amino acid in the peptide, 
		:type start_index: int 
		"""
		self._parent_proteins[parent_protein.get_id()]={
		'protein':parent_protein, 
		'start_index':start_index, 
		'end_index':end_index}

	def get_number_parent_protein(self) -> int:
		"""
		:return: The number of parent proteins this instance has 
		:rtype: int
		"""
		return len(self._parent_proteins) 

	def get_flanked_peptide(self,flank_len:int)-> Sequences: 
		"""
		:param flank_len: The length of the flanking regions 
		:type flank_len: int
		:return: A list of string containing the length of the peptide + the flanking region from \
		both the N and C terminal of the instance peptide, from all proteins.
		:rtype: Sequences
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
		""" Mapped the instance peptide to the parent protein and returned a 
		list of Numpy arrays where each array has a size of 1 by protein length. 
		within the protein the range representing the peptide is encoded as one while
		the rest is zero. 

		:return: A list of binary encoded arrays represent this mapping. 
		:rtype: MappedProtein
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
		:param length: The length of the non-presented peptides
		:type length: int
		:return: non-presented peptides from all the parent protein of the current peptide instance.  
		:rtype: Sequences
		"""
		results=[]
		for pro in self._parent_proteins.keys():
			results.append(self._parent_proteins[pro]['protein'].get_non_presented_peptide(self._parent_proteins[pro]['start_index'],self._parent_proteins[pro]['end_index'], length=length))
		return results
	
	def get_parent_proteins(self)->List[str]:
 		return list(set(self._parent_proteins.keys())) 

 	
	def is_child_of(self,pro_id:str)->bool:
		"""
		:param pro_id: is the protein id 
		:type pro_id: str
		:return: True if a user provided protein-id is a parent for the instance peptide, False otherwise 
		:rtype: bool
		"""
		return pro_id in self.get_parent_proteins()
	
	def get_pos_in_parent(self,pro_id:str)->Range:
		"""
		:param pro_id: the id of the parent protein 
		:type pro_id: str
		:raises ValueError: If the identifier is not a parent of the instance 
		:return: the start and end position of the instance peptide in the parent with the provided identifier
		:rtype: Range
		"""
		if not self.is_child_of(pro_id):
			raise ValueError(f"The provided protein id: {pro_id} is not a parent of the current peptide")
		# get the start and end elements 
		start=self._parent_proteins[pro_id]['start_index']
		end=self._parent_proteins[pro_id]['end_index']
		# return the results 
		return start, end 
	
	def get_parent(self,pro_id:str):
		"""
		:param pro_id: The protein identifer 
		:type pro_id: str
		:return: The parent protein that has an id matching the user defined pro_id
		:rtype: Protein
		"""
		if self.is_child_of(pro_id):
			return self._parent_proteins[pro_id]['protein']

	def get_number_of_parents(self)->int:
		"""
		:return: The number of instance parent proteins
		:rtype: int
		"""
		return len(self._parent_proteins)
	
	def get_n_terminal_flank_seq(self, flank_len : int)->List[str]:
		"""
		:param flank_len: the length of the flanking regions 
		:type flank_len: int
		:return: A list of string containing the sequences located upstream of the peptide in the parent protein. 
		:rtype: List[str]
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
		 
	def get_c_terminal_flank_seq(self, flank_len: int)->List[str]:
		"""
		:param flank_len: The length of the flanking regions 
		:type flank_len: int
		:return: A list of string containing the sequences located downstream of the peptide in the parent protein. 
		:rtype: [type]
		"""
		# if the peptide has no parent -> return a list with an empty string 
		if self.get_number_of_parents()==0:
			return ['']	
		#  else, we declare a list to allocate the results to it ==>  
		results: List[str]=[]
		# loop over all the parent to obtain the upstream sequences 
		for parent in self._parent_proteins.keys(): 
			# get the parent protein 
			parent_protein=self._parent_proteins[parent]
			# get the end of the protein 
			end_index=parent_protein['end_index']
			# get the protein length  
			protein_length=len(parent_protein['protein']) 		
			if end_index+flank_len>protein_length: # ==> avoid indexing outside of the parent protein
				results.append(parent_protein['protein'][end_index:protein_length])
			else: 
				results.append(parent_protein['protein'][end_index:end_index+flank_len])
		# return the results 
		return results 

	def add_org_2_parent(self,prot_name: str,org:str)->None:
		"""adds the source organism of one of the instance parent protein

		:param prot_name: The name of the protein, i.e. the identifier of the protein 
		:type prot_name: str
		:param org: The name of the organism 
		:type org: str
		:raises ValueError: Incase the provided protein is not a parent of the provided peptide 
		"""
		if not self.is_child_of(prot_name):
			raise ValueError(f'The provided protein id is not a parent of the current peptide')
		self._parent_proteins[prot_name]['protein'].set_org(org)
	
	def get_parents_org(self)->Organisms:
		"""
		:return: A list containing the name of each parent protein source organisms
		:rtype: Organisms
		"""
		results=[]
		for protein in self._parent_proteins.keys():
			results.append(self._parent_proteins[protein]['protein'].get_org())
		return results 
	
	def has_PTM(self)->bool: 
		"""check if the current peptide has a PTM in its sequence, this is done by checking for parentheses in the sequence. 

		:return: if the peptide has a PTM, the function return true, otherwise it return false. 
		:rtype: bool
		"""
		if '(' in self._peptide:
			return True
		return False
	
	def __len__(self)->int:
		"""
		:return: The length of the peptide 
		:rtype: int
		"""
		if '(' in self._peptide:
			temp_peptide=self._peptide
			while '(' in temp_peptide or ')' in temp_peptide:
				pre_seq=temp_peptide.split('(')[0]
				post_seq=")".join(temp_peptide.split(')')[1:])
				temp_peptide=pre_seq+post_seq
			return len(temp_peptide)
		else: 
			return len(self._peptide)
	
	def __getitem__(self,aa_index:int)->str:
		"""Bracket based indexing into the peptide sequence

		:param aa_index: the amino acid index in the peptide sequence, for example 3rd amino acid
		:type aa_index: int
		:return: the amino acid correspond to the provided index  
		:rtype: str
		"""
		return self._peptide[aa_index]
	
	def __str__(self)->str:
		"""
		:return: The sequence of the peptide 
		:rtype: str
		"""
		return self._peptide

	def __repr__(self)->str:
		"""
		:return: A string representation for the class 
		:rtype: str
		"""
		return f'A peptide instance with {self.get_number_of_parents()} parent proteins.'

		
		
		
		