#!/usr/bin/env python
"""
@author: Hesham ElAbd 
@brief: A representation of the Tissue reference expression values. 
@version: 0.0.1
"""
# load the modules 
from __future__ import annotations
import pandas as pd
import numpy as np 
# define the tissue class 
class Tissue: 
	"""
	@brief: a representation of tissue reference expression value. 
	"""
	def __init__(self, name: str, expression_table: pd.DataFrame, aux_proteins: pd.DataFrame = None)->Tissue: 
		"""
		@brief: Initialize a tissue reference expression values. 
		@param: name: the name of the tissue 
		@param: expression_table: a pandas dataframe with the following three columns 
				Transcript name, protein name, expression value in TPM  
		@param: aux_proteins: A table that contain the expression table of auxillary proteins that does not belong to the tissue per sa,
		for example, pathogen-derived genes.
		"""
		if name is not None:
			self._name=name
		else:
			raise ValueError(f'The tissue name can not be None')
		
		if not isinstance(expression_table, pd.DataFrame):
			raise ValueError(f'Expression table must be a pd.DataFrame instance, however, your input is of type: {type(expression_table)}')
		
		if expression_table.shape[1] != 3:  
			raise ValueError(f"The provided expression table must have three columns, however, your table have: {expression_table.shape[1]}")
		self._exp_map=expression_table
		# set the columns name 
		self._exp_map.columns=['Transcript_name', 'Protein_name','Expression_value']
		# parse the auxillary table 
		if aux_proteins is not None:
			if aux_proteins.shape[1]!=3:
				raise ValueError(f'The provided auxillary proteins table must have three columns, however, your table have: {aux_proteins.shape[1]}')
			aux_proteins.columns=['Transcript_name', 'Protein_name','Expression_value']
			# concatenate the results 
			self._exp_map=pd.concat([self._exp_map,aux_proteins],axis=0)
		return

	def get_transcript_expression(self,transcript_id: str)-> float: 
		"""
		@brief: return the expression value of the provided transcript 
		@param: transcript_id: the transcript id to retrive its expression value from the database
		"""
		if transcript_id not in self._exp_map.iloc[:,0].tolist(): 
			raise KeyError(f"The provided transcript id is not defined in the database: {transcript_id}")
		return self._exp_map.loc[self._exp_map.iloc[:,0]==transcript_id,:].iloc[0,2]
	
	def get_protein_expression(self,protein_id: str)-> float: 
		"""
		@brief: return the expression value of the provided protein id 
		@param: protein_id: the protein id to retrive its expression value from the database
		"""
		if protein_id not in self._exp_map.iloc[:,1].tolist(): 
			raise KeyError(f"The provided protein id is not defined in the database: {protein_id}")
		return self._exp_map.loc[self._exp_map.iloc[:,1]==protein_id,:].iloc[0,2]
	
	def get_name(self)->str:
		"""
		@brief: get the name of the tissue 
		"""
		return self._name
	
	def __str__(self)->str:
		"""
		@brief: compute a string representation for the class instance 
		"""
		return f'{self._name} with {self._exp_map.shape[0]} transcript expression value'
	
	def __repr__(self)->str:
		"""
		@brief: a representation for the class 
		"""
		return str(self)

		