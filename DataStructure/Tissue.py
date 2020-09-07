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
from typing import Union, List, Dict 
from IPTK.DataStructure.Database import CellularLocationDB, GeneExpressionDB

# define the tissue class 
class ExpressionProfile: 
	"""
	@brief: a representation of tissue reference expression value. 
	"""
	def __init__(self, name: str,
	 	expression_table: pd.DataFrame, 
		aux_proteins: pd.DataFrame = None)->ExpressionProfile: 
		"""
		@brief: create an expression profile instance from an expression table. 
		@param: name: the name of the tissue 
		@param: expression_table: a pandas dataframe with the following three columns 
				gene id, gene name, and the expression value.   
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
		self._exp_map.columns=['gene', 'gene_name','exp_value']
		# parse the auxillary table 
		if aux_proteins is not None:
			if aux_proteins.shape[1]!=3:
				raise ValueError(f'The provided auxillary proteins table must have three columns, however, your table have: {aux_proteins.shape[1]}')
			aux_proteins.columns=['gene', 'gene_name','exp_value']
			# concatenate the results 
			self._exp_map=pd.concat([self._exp_map,aux_proteins],axis=0)
		return

	def get_get_id_expression(self,gene_id: str)-> float: 
		"""
		@brief: return the expression value of the provided gene id.
		@param: gene_id: the gene id to retrive its expression value from the database
		"""
		if gene_id not in self._exp_map.iloc[:,0].tolist(): 
			raise KeyError(f"The provided gene id: {gene_id} is not defined in the database")
		return self._exp_map.loc[self._exp_map.iloc[:,0]==gene_id,:].iloc[0,2]
	
	def get_gene_name_expression(self,gene_name: str)-> float: 
		"""
		@brief: return the expression value of the provided gene name. 
		@param: gene_name: the gene name to retrive its expression value from the database.
		"""
		if gene_name not in self._exp_map.iloc[:,1].tolist(): 
			raise KeyError(f"The provided gene name: {gene_name} is not defined in the database.")
		return self._exp_map.loc[self._exp_map.iloc[:,1]==gene_name,:].iloc[0,2]
	
	def get_name(self)->str:
		"""
		@brief: get the name of the tissue 
		"""
		return self._name
	
	def __str__(self)->str:
		"""
		@brief: compute a string representation for the class instance 
		"""
		return f'{self._name} with an expression profile covering {self._exp_map.shape[0]} genes.'
	
	def __repr__(self)->str:
		"""
		@brief: a string representation for the class 
		"""
		return str(self)
## define the second class, annotated tissue 
class Tissue:
	def __init__(self, name: str, main_exp_value: GeneExpressionDB, 
	main_location: CellularLocationDB, aux_exp_value: GeneExpressionDB = None,
	aux_location: CellularLocationDB = None) -> Tissue:
		"""
	  	@brief: the initializer of the AnnotatedTissue class 
	  	@param: name: the name of the tissue 
	  	@param: main_exp_value: a GeneExpressionDB instace contain the gene expression accross different tissues 
	  	@param: main_location: a CellularLocationDB instance that contain the sub cellular locations for the proteins expressed in the tissue.  
	  	@param: aux_exp_value: a GeneExpressionDB instance that contain the expression table of auxillary proteins that does not belong to the tissue per sa,
	  	for example, pathogen-derived genes or extra-cellular matrix.
	  	@param: aux_location:  CellularLocationDB instance that contain the sub cellular locations for proteins that does not belong to the tissue of interest per sa
	  	for example, pathogen-derived proteins or media-added proteins. 
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
		@brief: provide the expresion profile of the current tissue 
		"""
		return self._exp_prof

	def get_subCellular_locations(self) ->CellularLocationDB:
		"""
		@brief: provide the sub-cellular localization of the proteins stored in current instance resources. 
		"""
		return self._cell_loc
	

		