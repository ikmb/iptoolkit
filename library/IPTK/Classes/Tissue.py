#!/usr/bin/env python
"""A representation of the Tissue used in an IP Experiment. 
"""
# load the modules 
from __future__ import annotations
import pandas as pd
import numpy as np 
from typing import Union, List, Dict 
from IPTK.Classes.Database import CellularLocationDB, GeneExpressionDB
# define the tissue class 
class ExpressionProfile: 
	"""a representation of tissue reference expression value. 
	"""
	def __init__(self, name: str,
	 	expression_table: pd.DataFrame, 
		aux_proteins: pd.DataFrame = None)->ExpressionProfile: 
		"""Create an expression profile instance from an expression table. 

		:param name: the name of the tissue 
		:type name: str
		:param expression_table: A dataframe with the following three columns gene id, gene name, and the expression value.   
		:type expression_table: pd.DataFrame
		:param aux_proteins: A table that contain the expression table of auxillary proteins that does not belong to the tissue per sa,
		for example, pathogen-derived genes, defaults to None
		:type: pd.DataFrame, optional
		:raises ValueError: if the tissue name is None
		:raises ValueError: if the expression table is not a pandas dataframe instance 
		:raises ValueError: if the expression table does not have the correct shape 
		:raises ValueError: if the aux_proteins does not have the correct shape
		:return: an expression profile instance 
		:rtype: ExpressionProfile
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
		self._exp_map.columns=['gene', 'gene_name','exp_value']
		# parse the auxillary table 
		if aux_proteins is not None:
			if aux_proteins.shape[1]!=3:
				raise ValueError(f'The provided auxillary proteins table must have three columns, however, your table have: {aux_proteins.shape[1]}')
			aux_proteins.columns=['gene', 'gene_name','exp_value']
			# concatenate the results 
			self._exp_map=pd.concat([self._exp_map,aux_proteins],axis=0)
		return

	def get_gene_id_expression(self,gene_id: str)-> float: 
		"""
		:param gene_id: the gene id to retrive its expression value from the database
		:type gene_id: str
		:raises KeyError: if the provided id is not defined in the instance table 
		:return: the expression value of the provided gene id.
		:rtype: float
		"""
		if gene_id not in self._exp_map.iloc[:,0].tolist(): 
			raise KeyError(f"The provided gene id: {gene_id} is not defined in the database")
		return self._exp_map.loc[self._exp_map.iloc[:,0]==gene_id,:].iloc[0,2]
	
	def get_gene_name_expression(self,gene_name: str)-> float: 
		"""
		:param gene_name: the gene name to retrive its expression value from the database
		:type gene_name: str
		:raises KeyError: if the provided id is not defined in the instance table 
		:return: the expression value of the provided gene name. 
		:rtype: float
		"""
		if gene_name not in self._exp_map.iloc[:,1].tolist(): 
			raise KeyError(f"The provided gene name: {gene_name} is not defined in the database.")
		return self._exp_map.loc[self._exp_map.iloc[:,1]==gene_name,:].iloc[0,2]

	def get_name(self)->str:
		"""
		:return: the name of the tissue where the expression profile was obtained
		:rtype: str
		"""
		return self._name

	def get_table(self)->pd.DataFrame:
		"""
		:return: return a table that contain the expression of all the transcripts in the current profile \
		including core and auxiliary proteins
		:rtype: pd.DataFrame
		"""
		return self._exp_map
	
	def __len__(self)->int: 
		"""
		:return: return the number of unique transcripts in the profile 
		:rtype: int
		"""
		return len(set(self._exp_map.shape[0]))

	def __str__(self)->str:
		""" 
		:return: a string representation for the class instance  
		:rtype: str
		"""
		return f'{self._name} with an expression profile covering {self._exp_map.shape[0]} genes.'
	
	def __repr__(self)->str:
		return str(self)

class Tissue:
	def __init__(self, name: str, main_exp_value: GeneExpressionDB, 
	main_location: CellularLocationDB, aux_exp_value: GeneExpressionDB = None,
	aux_location: CellularLocationDB = None) -> Tissue:
		"""The initializer of the Tissue class 

		:param name: The name of the tissue 
		:type name: str
		:param main_exp_value: A GeneExpressionDB instace containing the gene expression accross different tissues 
		:type main_exp_value: GeneExpressionDB
		:param main_location: A CellularLocationDB instance that contain the sub cellular locations for the proteins expressed in the tissue
		:type main_location: A CellularLocationDB 
		:param aux_exp_value: A GeneExpressionDB instance that contain the expression table of auxillary proteins that does not belong to the tissue per sa,
	  	for example, pathogen-derived genes or extra-cellular matrix, defaults to None.
		:type aux_exp_value: GeneExpressionDB, optional
		:param aux_location: CellularLocationDB instance that contain the sub cellular locations for proteins that does not belong to the tissue of interest per sa
	  	for example, pathogen-derived proteins or media-added proteins, defaults to None
		:type aux_location: CellularLocationDB, optional
		:raises KeyError: if the provided tissue name is not defined in the gene expression database 
		:return: [description]
		:rtype: Tissue
		"""
		if name not in main_exp_value.get_tissues():
			raise KeyError(f'The provided tissue name: {name} is not in the provided main expression database')
		# add the expression profile 
		if aux_exp_value is not None: 
			self._exp_prof: ExpressionProfile = ExpressionProfile(name=name,
			expression_table=main_exp_value.get_expression_in_tissue(name),
			aux_proteins=aux_exp_value.get_table())
		else: 
			self._exp_prof: ExpressionProfile = ExpressionProfile(
				name=name, expression_table=main_exp_value.get_expression_in_tissue(name))	
		# add the cellular location profile:  
		self._cell_loc = main_location
		if aux_location is not None:
			self._cell_loc.add_to_database(aux_location)
	
	def get_expression_profile(self)->ExpressionProfile: 
		"""
		:return: the expresion profile of the current tissue 
		:rtype: ExpressionProfile
		"""
		return self._exp_prof
	
	def get_name(self)->str:
		"""
		:return: the name of the tissue 
		:rtype: str
		"""
		return self._exp_prof.get_name()

	def get_subCellular_locations(self) ->CellularLocationDB:
		""" 
		:return: the sub-cellular localization of all the proteins stored in current instance resources. 
		:rtype: CellularLocationDB
		"""
		return self._cell_loc
	
	def __str__(self)->str:
		"""
		:return: a string representation for the current instance
		:rtype: str
		"""
		return f'{self._exp_prof.get_name()} with an associated expression profile covering: {len(self._exp_prof)} genes and a sub-cellular location covering: {len(self._cell_loc)} genes.'
		