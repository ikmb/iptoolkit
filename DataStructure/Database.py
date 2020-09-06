#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: a sequence database class handler
@date: 12.08.2020 
@version: 0.0.1
"""
# load the modules: 
from __future__ import annotations
from Bio import SeqIO 
import os
# define the class 
class SeqDB: 
	"""
	A Fasta sequence database class handler. 
	"""
	def __init__(self, path2fasta: str) -> SeqDB: 
		"""
		@brief: load a fasta file and constructs a lock up dictionary where sequence ids are 
		keys and sequences are values. 
		@param: path2fasta, the path to load the fasta database 
		"""
		if not os.path.exists(path2fasta): 
			raise ValueError(f"Your provided path: {path2fasta} does not exist!")
		# define a private variable to hold the results 
		self._seqs=dict() 
		seq_gen=SeqIO.parse(path2fasta,'fasta')
		for seq in seq_gen: 
			key=seq.id.split('|')[1]
			self._seqs[key]=str(seq.seq) 
		return

	def get_seq(self,protein_id): 
		"""
		@brief: returns the corresponding sequence if the provided protein-id is defined in the database
		otherwise it raise a keyError.
		@param: protein_id, the protein id to retrive its sequence. 
		"""
		if protein_id not in self._seqs.keys(): 
			raise KeyError(f"Your provided protein id:{protein_id} is not in the database" )
		
		return self._seqs[protein_id] 

	def __len__(self)->int:	
		"""
		@brief: a magic function to return the number of the sequence in the database 
		"""
		return len(self._seqs)

	def __getitem__(self,protein_id: str)->str:
		"""
		@brief: a magic function for accessing elements from the database 
		"""
		return self.get_seq(protein_id)
	
	def __str__(self)->str:
		"""
		@brief: return a string representation for the class 
		"""
		return f'A sequence database with {len(self)} sequence'

	def __repr__(self)->str:
		"""
		@brief: return the representation of the database instance 
		"""
		return str(self)
	
	def has_sequence(self, sequence_name: str)->bool:
		"""
		@brief: check if the provided sequence name is a part of the sequences stored in the database or not
		@param: sequence_name: the name of the sequence 
		"""
		return sequence_name in self._seqs.keys()