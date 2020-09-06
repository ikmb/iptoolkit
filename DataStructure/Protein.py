#!/usr/bin/env python 
"""
@author: Hesham ElAbd 
@brief: A representation of a protein that has been infered from an IP experiment. 
@version: 0.0.1
"""
# load the modules:
from __future__ import annotations 
import numpy as np 
import pandas as pd
from IPTK.Utils.Types import Index
# define some types 
class Protein: 
	"""
	@briefA representation of a protein that has been infered from an IP experiment. 
	"""
	def __init__(self, prot_id : str, seq: str, org: str=None) -> Protein:
		"""
		@brief: construct a portein instance. 
		@param: prot_id: the id of the protein 
		@param: seq: the sequecne of the protein 
		@param: org: the organism from which the protein originated, default is None
		""" 
		self._id=prot_id
		self._seq=seq
		self._org=org
		
	def get_id(self) -> str:
		"""
		@brief return the protein identifier.
		""" 
		return self._id 
	
	def get_seq(self) -> str:
		"""
		@brief return the protein sequence. 
		""" 
		return self._seq
	
	def get_peptides_map(self, start_idxs: Index, end_idxs:Index)->np.ndarray: 
		"""
		@brief: compute a coverage over the protein sequence
		@detail returns a numpy array with shape of 1 by the length of the protein.
		Every element in the array donates the number of times, It has been observed in the experiment.
		@param start_idxs a list of integers representing the start position
		@param end_idxs a list of integers representing the end position 
		@note start_indxs and end_idxs MUST be of equal length, where each pair of elements define
		the boundary of an identified peptide.  
		"""
		if len(start_idxs) != len(end_idxs):
			raise ValueError(f"""Annotation list MUST have an equal length where your lists have length : {len(start_idxs)}, {(end_idxs)} respectivily.""")
		prote_backbone=np.zeros(shape=(1,len(self)))
		for i,j in zip(start_idxs,end_idxs): 
			prote_backbone[0,i:j]+=1 
		return prote_backbone 
	
	def get_non_presented_peptide(self, exc_reg_s_idx: int, exc_reg_e_idx: int, length: int) -> str: 
		"""
		@brief sample a peptide from the protein sequences where the sampled peptides is not part of the 
		experimentally identified regions. 
		@param exc_reg_s_idx the start point in the reference protein sequence of the experimentally identified peptide. 
		@param exc_reg_e_idx the end point in the reference protein sequence of the experimentally identified peptide. 
		@param length the length the non-presented peptides. 
		@return a substring of the instance sequecnes. 
		"""
		if length >= len(self):
			raise ValueError(f'''The generated peptide must be shorter than the parent protein length,
		 	Provided length: {length}  While, the parent length: {len(self)}''')
		if length <=0:
			raise ValueError(f'The generated peptide length MUST be  a positive integer, your input value is: {length}')
		try_flag=True
		while try_flag: 
			anc_point=np.random.randint(low=0,high=len(self._seq)-length)
			if anc_point < exc_reg_s_idx or anc_point > exc_reg_e_idx: 
				neg_seq=self._seq[anc_point:anc_point+length] 
				try_flag=False
		return neg_seq
	
	def get_org(self)->str:
		"""
		@brief return the organims in which this instance protein belong.
		"""
		return self._org

	def set_org(self,org: str)->None:
		"""
		@brief a post-instantitation mechanims to set the organism for which the protein belong. 
		"""
		self._org=org

	def __getitem__(self,aa_index: int ) -> str:
		"""
		@brief: return the amino acid with the corresponding index as the user provided index 
		@param: aa_index: the amino acid index, for example, 3rd amino acid in the protein 5th amino acid in the protein
		""" 
		return self._seq[aa_index] 
	
	def __len__(self)->int:
		"""
		@brief: get the length of the protein instance 
		"""
		return len(self._seq)
	
	def __str__(self)->str:
		"""
		@brief: return the string of the protein length 
		"""
		return self._seq

	def __repr__(self)->str:
		"""
		@brief: a string representation of the protein  
		"""
		return f'A Protein instance with a length of {len(self)}'