#!/usr/bin/env python 
"""This submodule defines a collection of container classes that are used through the library
"""
# load the modules: 
from __future__ import annotations
from Bio import SeqIO 
import os
import time 
import pandas as pd 
from typing import List, Dict 
from Bio import SeqIO
# define the class 
class SeqDB: 
	"""Load a FASTA file and constructs a lock up dictionary where sequence ids are  keys and sequences are values. 
	"""
	def __init__(self, path2fasta: str) -> SeqDB: 
		"""The class constructor 

		:param path2fasta: the path to load the fasta database 
		:type path2fasta: str
		:raises FileNotFoundError: In case the provided path does not exist or the file can not be found. 
		:return: a sequence database instance
		:rtype: SeqDB
		"""
		if not os.path.exists(path2fasta): 
			raise FileNotFoundError(f"Your provided path: {path2fasta} does not exist!")
		# define a private variable to hold the results 
		self._seqs=dict() 
		for seq in SeqIO.parse(path2fasta,'fasta'): 
			key=seq.id.split('|')[1]
			self._seqs[key]=str(seq.seq) 
		return

	def get_seq(self,protein_id:str)->str: 
		"""returns the corresponding sequence if the provided protein-id is defined in the database.

		:param protein_id: The protein id to retrive its sequence, CASE SENSITIVE!!. 
		:type protein_id: str
		:raises KeyError: If the provided protein does not exist in the database 
		:return: the protein sequence
		:rtype: str
		"""
		if protein_id not in self._seqs.keys(): 
			raise KeyError(f"Your provided protein id:{protein_id} is not in the database" )
		return self._seqs[protein_id] 

	def __len__(self)->int:	
		"""A magic function to return the number of the sequence in the database 

		:return: the number of the sequences in the database 
		:rtype: int
		"""
		return len(self._seqs)

	def __getitem__(self,protein_id: str)->str:
		"""a magic function for indexing through the elements of the library.   

		:param protein_id: the protein id 
		:type protein_id: str
		:return: the sequence of the provided id 
		:rtype: str
		"""
		return self.get_seq(protein_id)
	
	def __str__(self)->str:
		"""return a string representation for the class 

		:return: return a string representation for the class
		:rtype: str
		"""
		return f'A sequence database with {len(self)} sequence'

	def __repr__(self)->str:
		"""return the representation of the database instance.

		:return: a string representation of the class 
		:rtype: str
		"""
		return str(self)
	
	def has_sequence(self, sequence_id: str)->bool:
		"""check if the provided sequence id is an element of the database or not		

		:param sequence_name: The id of the sequence, CASE SENSITIVE!!. 
		:type sequence_name: str
		:return: True if the database has this id, False otherwise. 
		:rtype: bool
		"""
		return sequence_id in self._seqs.keys()

class CellularLocationDB:
	"""The class provides an API to access the cellular location information from a database \
	that follows the structure of the Human Proteome Atlas sub-cellular location database. See https://www.proteinatlas.org/about/download \
	for more details. 
	"""
	def __init__(self, 
				path2data: str = "https://www.proteinatlas.org/download/subcellular_location.tsv.zip",
				sep: str = '\t')->CellularLocationDB:
		"""	The class constructor

		:param path2data: the path to load the database or the URL to download it. defaults to https://www.proteinatlas.org/download/subcellular_location.tsv.zip
		:type path2data: str
		:param sep: The table seperator, defaults to '\t'
		:type sep: str, optional
		:raises IOError: Incase the provided table could not be loaded 
		:return: an instance of the class 
		:rtype: CellularLocationDB
		"""
		try: 
			self._table: pd.DataFrame = pd.read_csv(path2data,sep=sep)
		except Exception as exp: 
			raise IOError(f'While loading the sub cellular location dataset the following error: {exp} was encountered.')
	
	def get_genes(self)-> List[str]:
		"""return a list of all gene ids in the dataset 

		:return: all genes ids currently defined in the database 
		:rtype: List[str]
		"""
		return self._table.iloc[:,0].tolist()
	
	def get_gene_names(self)->List[str]:
		"""return a list of all gene names in the dataset 

		:return: the names of all genes in the database 
		:rtype: List[str]
		"""
		return self._table.iloc[:,1].tolist()
	
	def get_main_location(self, gene_id: str = None, corresponds= None)->List[str]:
		"""Return the main location(s) of the provided gene id or gene name. \ 
		If both gene Id and gene name are provided, gene_id has a higher precedence 

		:param gene_id: The id of the gene of interest, defaults to None
		:type gene_id: str, optional
		:param gene_name: The name of the gene of interest, defaults to None
		:type gene_name: [type], optional
		:raises ValueError: if both gene_id and gene_name are None
		:raises KeyError: if gene_id is None and gene_name is not in the database 
		:raises KeyError: if gene_name is None and gene_id is not in the database 
		:raises RuntimeError: Incase an error was encountered while retriving the element from the database
		:return: the main location where the protein that corresponds to the provided name or id is located. 
		:rtype: List[str]
		"""
		if gene_id is None and corresponds is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if corresponds is not None:
			if corresponds not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {corresponds} is not in the database')
		# get the main location 
		try:
			if gene_id is not None:  
				main_location: str = self._table.loc[self._table.iloc[:,0]==gene_id,]['Main location'] 
			else: 
				main_location: str = self._table.loc[self._table.iloc[:,1]==corresponds,]['Main location'] 
		except Exception as exp: 
			raise RuntimeError(f'While getting the main cellular location the following error was encounter: {exp}')
		# split the locations 
		if main_location.isna().any():
			locations: List[str] = ['UNK']
		else: 
			locations: List[str] = main_location.tolist()[0].split(';')
		# return the resuts 
		return locations 
	
	def get_approved_location(self, gene_id: str = None, gene_name= None) -> List[str]:
		"""return the location of the provided gene id or gene name 

		:param gene_id: the id of the gene of interest, defaults to None
		:type gene_id: str, optional
		:param gene_name: the name of gene of interest, defaults to None
		:type gene_name: [type], optional
		:raises ValueError: if both gene_id and gene_name are None
		:raises KeyError: if gene_id is None and gene_name is not in the database 
		:raises KeyError: if gene_name is None and gene_id is not in the database 
		:raises RuntimeError: Incase an error was encountered while retriving the element from the database. 
		:return: The approved location where the protein that corresponds to the provided name or id is located. 
		:rtype: List[str]
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		# get the approved locations 
		try:
			if gene_id is not None:  
				approved_locations: str = self._table.loc[self._table.iloc[:,0]==gene_id,]['Approved']  
			else: 
				approved_locations: str = self._table.loc[self._table.iloc[:,1]==gene_name,]['Approved']
		except Exception as exp: 
			raise RuntimeError(f'While getting the approved locations the following error was encounter: {exp}')
		# split the locations 
		if approved_locations.isna().any(): 
			locations: List[str] =  ['UNK']
		else: 
			locations: List[str] =  approved_locations.tolist()[0].split(';')	
		# return the resuts 
		return locations 
	
	def get_go_names(self,gene_id: str = None, gene_name= None )->List[str]:
		"""return the location of the provided gene id or gene name 

		:param gene_id: the id of the gene of interest , defaults to None
		:type gene_id: str, optional
		:param gene_name: the name of the gene of interest , defaults to None
		:type gene_name: [type], optional
		:raises ValueError: if both gene_id and gene_name are None
		:raises KeyError: if gene_id is None and gene_name is not in the database 
		:raises KeyError: if gene_name is None and gene_id is not in the database 
		:raises RuntimeError: incase an error was encountered while retriving the element from the database.
		:return: The gene ontology, GO, location where the protein that corresponds to the provided name or id is located. 
		:rtype: List[str]
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		# get gene ontology IDs
		try:
			if gene_id is not None:  
				go_ids : str = self._table.loc[self._table.iloc[:,0]==gene_id,]['GO id'] 
			else: 
				go_ids: str = self._table.loc[self._table.iloc[:,1]==gene_name,]['GO id']
		except Exception as exp: 
			raise RuntimeError(f'While getting the approved locations the following error was encounter: {exp}')
		# split the locations 
		if go_ids.isna().any(): 
			_ids: List[str] =  ['UNK']
		else: 
			_ids: List[str] =  go_ids.tolist()[0].split(';')	
		# return the resuts 
		return _ids

	def get_table(self) -> pd.DataFrame:
		"""return the instance table 

		:return: the location table of the instance. 
		:rtype: pd.DataFrame
		"""
		return self._table

	def add_to_database(self,genes_to_add: CellularLocationDB )->None:
		"""adds the the location of more proteins to the database. 

		:param genes_to_add: a CellularLocationDB instance containing the genes that shall be added to the database.   
		:type genes_to_add: CellularLocationDB
		:raises ValueError: if the genes_to_add to the database are already defined in the database
		:raises RuntimeError: incase any other error has been encountered while merging the tables.
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
	
	# define some magic functions 
	def __len__(self)->int:
		"""return the number of unique genes in the current instance
		
		:return: the number of unique genes in the current instance
		:rtype: int
		"""
		return len(set(self._table.iloc[:,0]))
		
	def __str__(self)->str: 
		""" return a string representation of the instance 

		:return: a string representation for the instance 
		:rtype: str
		"""
		return f' A sub-cellular localization database covering {len(self)} genes'

	def __repr__(self)->str: 
		return str(self)


class GeneExpressionDB: 
	"""The class provides an API to access gene expression data stored in table that follows the same structure as \ 
	the Human proteome Atlas Normalized RNA Expression see  https://www.proteinatlas.org/about/download for more details 
	"""
	def __init__(self, path2data: str = 'https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip', 
					sep: str = '\t')->GeneExpressionDB: 
		"""The class constructor 

		:param path2data: The path to the expression table, defaults to, 'https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip'. 
		:type path2data: str, optional 
		:param sep: The separator for the input table, defaults to '\t'
		:type sep: str, optional
		:raises IOError: Incase the provided table could not be load. 
		:return: an instance of the class 
		:rtype: GeneExpressionDB
		"""
		try: 
			self._table: pd.DataFrame = pd.read_csv(path2data,sep=sep)
		except Exception as exp: 
			raise IOError(f'While loading the RNA expression dataset the following error: {exp} was encountered.')

	def get_genes(self)->List[str]:
		"""returns a list of the UNIQUE gene ids currently in the database. 

		:return: The list of the UNIQUE gene ids currently in the database 
		:rtype: List[str]
		"""
		return list(set(self._table.iloc[:,0].tolist()))
	
	def get_gene_names(self)->List[str]:
		"""returns a list of the UNIQUE gene names currently in the database 

		:return: A list of the UNIQUE gene names currently in the database 
		:rtype: List[str]
		"""
		return list(set(self._table.iloc[:,1].tolist()))
	
	def get_tissues(self)->List[str]:
		"""return a list of the tissues in the current database
		
		:return: A list containing the names of the UNIQUE tissues in the database.  
		:rtype: List[str]
		"""
		return list(set(self._table.Tissue))

	def get_expression(self, gene_name: str = None, gene_id: str = None) ->pd.DataFrame: 
		"""Return a table summarizing the expression of the provided gene name or gene id accross different tissues. 

		:param gene_id: the id of the gene of interest, defaults to None
		:type gene_id: str, optional
		:param gene_name: the name of the gene of interest, defaults to None
		:type gene_name: [type], optional
		:raises ValueError: if both gene_id and gene_name are None
		:raises KeyError: if gene_id is None and gene_name is not in the database 
		:raises KeyError: if gene_name is None and gene_id is not in the database 
		:raises RuntimeError: incase an error was encountered while retriving the elements from the database
		:return: A table summarizing the expression of the provided gene accross all tissues in the database
		:rtype: pd.DataFrame
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		
		if self._table.loc[self._table.iloc[:,0]==gene_name]['Gene name']!=gene_name:
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
		"""return the expression profile of the provided tissue  

		:param tissue_name: The name of the tissue 
		:type tissue_name: str
		:raises KeyError: Incase the provided tissue is not defined in the database
		:raises RuntimeError: In case an error was encountered while generating the expression profile. 
		:return: A table summarizing the expression of all genes in the provided tissue.
		:rtype: pd.DataFrame
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
		"""return a table containing the expression value of all the genes accross all tissues in the current instance

		:return: The expression of all genes accross all tissues in the database. 
		:rtype: pd.DataFrame
		"""
		return self._table[['Gene', 'Gene name', 'NX']]
		
	def __len__(self)->int: 
		"""return the number of tissues in the database

		:return: the number of tissues in the database
		:rtype: int
		"""
		return len(self.get_tissues())
	
	def __str__(self)-> str:
		"""A string representation of the class

		:rtype: str
		"""
		return f'an expression database covering: {len(self.get_genes())} genes accross {len(self)} tissues'
	
	def __repr__(self) ->str: 
		return str(self)

class OrganismDB:
	"""Extract information about the source organsim of a collection of protein sequences\ 
	from a fasta file and provides an API to query the results. \
	The function expect the input fasta file to have headers written in the UNIPROT format. 
	"""
	def __init__(self,path2Fasta:str)->OrganismDB:
		""" The class constructor 
	
		:param path2Fasta: The path to a fasta sequence database to obtain the protein sequences
		:type path2Fasta: str
		:raises IOError: Incase the database could not be loaded 
		:return: A new OrganismDB instance 
		:rtype: OrganismDB
		"""
		try: 
			print(f"Reading the input fasta file, ..., started at: {time.ctime()}")
			seq_gene=SeqIO.parse(path2Fasta,'fasta')
		except Exception as exp: 
			raise IOError(f'While loading your input database: {path2Fasta}, the following error was encountered: {exp}')
		# define a dict that hold the map 
		self._map: Dict[str,str] = dict()
		# fill the elements in the map 
		for seq in seq_gene: 
			temp_name_org: List[str]=seq.id.split('|')
			self._map[temp_name_org[1]]= temp_name_org[2].split('_')[1]
	
	def get_unique_orgs(self)->List[str]:
		"""Get the number of unique organisms in the database
		
		:return: a list of all unique organisms in the current instance 
		:rtype: List[str]
		"""
		return list(set(self._map.values()))
	
	def get_number_protein_per_organism(self)->pd.DataFrame:
		"""Provides a table containing the number of proteins per organism. 
		
		:return: A table containing the number of proteins per organism
		:rtype: pd.DataFrame
		"""
		print(f"computing the number of proteins per organism: ..., started at {time.ctime()}")
		unique_orgs: List[str] = self.get_unique_orgs()
		# create a dict counter 
		counter: Dict[str,int]= dict()
		# Initialize the counts to 0s 
		for org in unique_orgs: 
			counter[org]=0
		# fill the dict with counts 
		for key in self._map.keys(): 
			counter[self._map[key]]+=1
		# chaning the map from int -> List[int] to fit the dataframe requirements 
		for key in counter: 
			counter[key]=[counter[key]]
		# create the dataframe 
		res: pd.DataFrame = pd.DataFrame(counter).T 
		# add the org-names as keys  
		res['Names'] = res.index.tolist()
		# return the results 
		return res 
	
	def get_org(self,prot_id:str)->str: 
		"""return the parent organism of the provided protein identifer

		:param prot_id: the id of the protein of interest 
		:type prot_id: str
		:raises KeyError: incase the provided identifier is not in the database 
		:return: the name of the parent organism, i.e. the source organism. 
		:rtype: str
		"""
		if prot_id not in self._map.keys(): 
			raise KeyError(f'The provided key is not defined in the database')
		return self._map[prot_id]
	
	def __len__(self)->int: 
		"""return the number of protein id parent organism pairs currently stored  in the database.
		
		:rtype: int
		"""
		return len(self._map)
	
	def __str__(self)->str: 
		""" return a string representation of the class 
				:rtype: str
		"""
		return str(f'An organism database with: {len(self)} entries.')

	def __repr__(self)->str: 
		return str(self)

		

	
