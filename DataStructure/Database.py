#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: a sequence database class handler
@version: 0.0.1
"""
# load the modules: 
from __future__ import annotations
from Bio import SeqIO 
import os
import pandas as pd 
from typing import List 
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

class CellularLocationDB:
	"""
	@brief: The class provide an API to access the human Proteome Atlas Database covering 
	the sub-cellular location of the proteins.  
	@see: https://www.proteinatlas.org/about/download
	"""
	def __init__(self, 
				path2data: str = '../resources/subcellular_location.tsv',
				sep: str = '\t')->CellularLocationDB:
		"""
		@param: path2data: the path to the subcellular locations 
		@param: sep:  the separator for the database 
		"""
		try: 
			self._table: pd.DataFrame = pd.read_csv(path2data,sep=sep)
		except Exception as exp: 
			raise IOError(f'While loading the sub cellular location dataset the following error: {exp} was encountered.')
	
	def get_genes(self)-> List[str]:
		"""
		@brief: return a list of all gene ids in the dataset 
		"""
		return self._table.iloc[:,0].tolist()
	
	def get_gene_names(self)->List[str]:
		"""
		@brief: return a list of all gene names in the dataset 
		"""
		return self._table.iloc[:,1].tolist()
	
	def get_main_location(self, gene_id: str = None, gene_name= None)->List[str]:
		"""
		@brief: return the main location(s) of the provided gene id or gene name 
		@param: gene_id: the id of the gene of interest 
		@param: gene_name: the name of the gene of interest 
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		if self._table.loc[self._table._iloc[:,0]==gene_name]['Gene name']!=gene_name:
			raise ValueError('The provided gene names: {gene_name} does not match the provided gene id: {gene_id}')
		
		# get the main location 
		try:
			if gene_id is not None:  
				main_location: str = self._table.loc[self._table.iloc[:,0]==gene_id,]['Main location'].tolist()[0]  
			else: 
				main_location: str = self._table.loc[self._table.iloc[:,1]==gene_name,]['Main location'].tolist()[0]  
		except Exception as exp: 
			raise RuntimeError(f'While getting the main cellular location the following error was encounter: {exp}')
		# split the locations 
		locations: List[str] =  main_location.split(';')
		# return the resuts 
		return locations 
	
	def get_approved_location(self, gene_id: str = None, gene_name= None) -> List[str]:
		"""
		@brief: return the location of the provided gene id or gene name 
		@param: gene_id: the id of the gene of interest 
		@param: gene_name: the name of the gene of interest 
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		if self._table.loc[self._table._iloc[:,0]==gene_name]['Gene name']!=gene_name:
			raise ValueError('The provided gene names: {gene_name} does not match the provided gene id: {gene_id}')
		
		# get the approved locations 
		try:
			if gene_id is not None:  
				approved_locations: str = self._table.loc[self._table.iloc[:,0]==gene_id,]['Approved'].tolist()[0]  
			else: 
				approved_locations: str = self._table.loc[self._table.iloc[:,1]==gene_name,]['Approved'].tolist()[0]  
		except Exception as exp: 
			raise RuntimeError(f'While getting the approved locations the following error was encounter: {exp}')
		# split the locations 
		locations: List[str] =  approved_locations.split(';')
		# return the resuts 
		return locations 
	
	def get_go_names(self,gene_id: str = None, gene_name= None ):
		"""
		@brief: return the location of the provided gene id or gene name 
		@param: gene_id: the id of the gene of interest 
		@param: gene_name: the name of the gene of interest 
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		if self._table.loc[self._table._iloc[:,0]==gene_name]['Gene name']!=gene_name:
			raise ValueError('The provided gene names: {gene_name} does not match the provided gene id: {gene_id}')
		
		# get gene ontology IDs
		try:
			if gene_id is not None:  
				go_ids : str = self._table.loc[self._table.iloc[:,0]==gene_id,]['GO id'].tolist()[0]  
			else: 
				go_ids: str = self._table.loc[self._table.iloc[:,1]==gene_name,]['GO id'].tolist()[0]  
		except Exception as exp: 
			raise RuntimeError(f'While getting the approved locations the following error was encounter: {exp}')
		# split the locations 
		locations: List[str] =  go_ids.split(';')
		# return the resuts 
		return locations

	def get_table(self) -> pd.DataFrame:
		"""
		@brief: return the instance table 
		"""
		return self._table

	def add_to_database(self,genes_to_add: CellularLocationDB )->None:
		"""
		@brief: add the the location of more proteins to the database. 
		@param: genes_to_add: a CellularLocationDB instance containing the genes that shall be added to 
		the database.   
		"""
		# check that the genes in provided database do not exist in the current instance.
		genes: List[str] = genes_to_add.get_genes()

		for gene in genes: 
			if gene in self.get_genes():
				raise ValueError(f'The provided gene: {gene} is already defined in current instance table')
		# add the genes to the database, concat the databases  
		try: 
			self._table=pd.concat([self._table, genes_to_add.get_table()],axis=1)
		except Exception as exp: 
			raise RuntimeError(f'While combing the database, the following error was encountered: {exp}')
















class GeneExpressionDB: 
	"""
	@brief: provides an API to access gene expression data stored in RNA consensus tissue gene data from Human proteome Atlas 
	@see: https://www.proteinatlas.org/about/download for more details 
	"""
	def __init__(self, path2data: str = '../resources/rna_tissue_consensus.tsv', 
					sep: str = '\t')->GeneExpressionDB: 
		"""
		@param: path2data: the path to the subcellular locations 
		@param: sep:  the separator for the input table 
		"""
		try: 
			self._table: pd.DataFrame = pd.read_csv(path2data,sep=sep)
		except Exception as exp: 
			raise IOError(f'While loading the RNA expression dataset the following error: {exp} was encountered.')

	def get_genes(self)->List[str]:
		"""
		@brief: return a list of the gene ids currently in the database 
		"""
		return list(set(self._table.iloc[:,0].tolist()))
	
	def get_gene_names(self)->List[str]:
		"""
		@brief: return a list of the gene names currently in the database 
		"""
		return list(set(self._table.iloc[:,1].tolist()))
	
	def get_tissues(self)->List[str]:
		"""
		@brief: return a list of the tissues in the current database
		"""
		return list(set(self._table.Tissue))

	def get_expression(self, gene_name: str = None, gene_id: str = None) ->pd.DataFrame: 
		"""
		@brief: return a table summarizing the expression of the provided gene name or gene id accross different tissues. 
		@param: gene_id: the id of the gene of interest 
		@param: gene_name: the name of the gene of interest 
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		if self._table.loc[self._table._iloc[:,0]==gene_name]['Gene name']!=gene_name:
			raise ValueError('The provided gene names: {gene_name} does not match the provided gene id: {gene_id}')
		# get the expression accross all tissues 			
		try: 
			if gene_id is not None:  
				return self._table.loc[self._table.iloc[:,0]==gene_id]
			else:
				return self._table.loc[self._table.iloc[:,1]==gene_name]
		except Exception as exp:
			raise RuntimeError(f'The following error has been encountered while prepearing the table: {exp}') 
			
	def get_expression_in_tissue(self, tissue_name:str)-> pd.DataFrame:
		"""
		@brief: return the expression profile of the provided gene
		@param: tissue_name: the name of the provided tissue
		"""
		if tissue_name not in self.get_tissues():
			raise KeyError(f"The provided tissue: {tissue_name} is not support, currently supported tissue are:  {','.join(self.get_tissues())} ")
		# return the table 
		try: 
			table: pd.DataFrame = self._table.loc[self._table.iloc[:,2]==tissue_name]
			res: pd.DataFrame = table[['Gene','Gene name','NX']]
			return res
		except Exception as exp:
			raise RuntimeError(f'While loading the database the following error was encountered: {exp}')
	
	def get_table(self)->pd.DataFrame:
		"""
		@brief: return a table containing the expression value of all the genes in the database
		"""
		return self._table[['Gene', 'Gene name', 'NX']]

	def __len__(self)->int: 
		"""
		@brief: return the number of tissues in the database
		"""
		return len(self.get_tissues())
	
	def __str__(self)-> str:
		"""
		@brief: a string representation of the class
		"""
		return f'an expression database covering: {len(self.get_genes())} genes accross {len(self)} tissues'
	
	def __repr__(self) ->str: 
		return str(self)


	
	