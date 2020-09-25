#!/usr/bin/env python 
"""A representation of a protein that has been inferred from an IP experiment. 
"""
# load the modules:
from __future__ import annotations 
import numpy as np 
import pandas as pd
from IPTK.Utils.Types import Index
# define some types 
class Protein: 
	"""representation of a protein that has been infered from an IP experiment. 
	"""
	def __init__(self, prot_id : str, seq: str, org: str=None) -> Protein:
		"""construct a portein instance. 

		:param prot_id: the id of the protein 
		:type prot_id: str
		:param seq:  the sequence of the protein 
		:type seq: str
		:param org: the organism from which the protein originated, defaults to None
		:type org: str, optional
		:return: a Protein instance
		:rtype: Protein
		"""
		self._id=prot_id
		self._seq=seq
		self._org=org
		
	def get_id(self) -> str:
		"""
		:return: return the protein identifier.
		:rtype: str
		"""
		return self._id 
	
	def get_seq(self) -> str:
		"""
		:return: the protein sequence.
		:rtype: str
		"""
		return self._seq
	
	def get_peptides_map(self, start_idxs: Index, end_idxs:Index)->np.ndarray:
		"""compute a coverage over the protein sequence
		
		:param start_idxs: a list of integers representing the start positions
		:type start_idxs: Index 
		:param end_idxs: a list of integers representing the end positions
		:type start_idxs: Index
		:raises ValueError: if start_indxs and end_idxs MUST be of equal length are not of equal length
		:return: A numpy array with shape of 1 by the length of the protein where every element in the array donates the number of times, It has been observed in the experiment.
		:rtype: np.ndarray
		"""
		if len(start_idxs) != len(end_idxs):
			raise ValueError(f"""Annotation list MUST have an equal length where your lists have length : {len(start_idxs)}, {(end_idxs)} respectively.""")
		prote_backbone=np.zeros(shape=(1,len(self)))
		for i,j in zip(start_idxs,end_idxs): 
			prote_backbone[0,i:j]+=1 
		return prote_backbone 
	
	def get_non_presented_peptide(self, exc_reg_s_idx: int, exc_reg_e_idx: int, length: int) -> str: 
		"""sample a peptide from the protein sequences where the sampled peptides is not part of the \
		experimentally identified regions. 
		
		:param exc_reg_s_idx: the start point in the reference protein sequence of the experimentally identified peptide. 
		:type exc_reg_s_idx: int
		:param exc_reg_e_idx: the end point in the reference protein sequence of the experimentally identified peptide. 
		:type exc_reg_e_idx: int
		:param length: length the non-presented peptides. 
		:type length: int
		:raises ValueError: if the length of the peptide is bigger than the protein length 
		:raises ValueError: if the length of the peptide is smaller than or equal to zero 
		:return: a substring of the instance sequence
		:rtype: str
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
		:return: the organism in which this instance protein belong.
		:rtype: str
		"""
		return self._org

	def set_org(self,org: str)->None:
		"""a post-instantitation mechanism to set the organism for which the protein belong.
	
		:param org: the name of the organism
		:type org: str
		"""
		self._org=org

	def __getitem__(self,aa_index: int ) -> str:
		"""
		:param aa_index:  the amino acid index, for example, 3rd amino acid in the protein 5th amino acid in the protein
		:type aa_index: int
		:return: the amino acid with the corresponding index as the user provided index
		:rtype: str
		"""
		return self._seq[aa_index] 
	
	def __len__(self)->int:
		"""
		:return: the length of the instances sequence  
		:rtype: int
		"""
		return len(self._seq)
	
	def __str__(self)->str:
		"""
		:return: the instance's sequence 
		:rtype: str
		"""
		return self._seq

	def __repr__(self)->str:
		return f'A Protein instance with a length of {len(self)}'