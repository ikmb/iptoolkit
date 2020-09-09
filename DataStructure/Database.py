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
from typing import List, Dict 
from Bio import SeqIO
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
				path2data: str ,
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
		@note: if both gene Id and gene name are provided, both gene_id has a higher precedence 
		"""
		if gene_id is None and gene_name is None:
			raise ValueError(f'gene_name and gene_id can not be None')
		
		if gene_id is not None:
			if gene_id not in self.get_genes():
				raise KeyError(f'The provided gene id: {gene_id} is not in the database')
		
		if gene_name is not None:
			if gene_name not in self.get_gene_names():
				raise KeyError(f'The provided gene name : {gene_name} is not in the database')
		# get the main location 
		try:
			if gene_id is not None:  
				main_location: str = self._table.loc[self._table.iloc[:,0]==gene_id,]['Main location'] 
			else: 
				main_location: str = self._table.loc[self._table.iloc[:,1]==gene_name,]['Main location'] 
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
		"""
		@brief: return the location of the provided gene id or gene name 
		@param: gene_id: the id of the gene of interest 
		@param: gene_name: the name of the gene of interest 
		@note: if both gene Id and gene name are provided, both gene_id has a higher precedence 
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
	
	def get_go_names(self,gene_id: str = None, gene_name= None ):
		"""
		@brief: return the location of the provided gene id or gene name 
		@param: gene_id: the id of the gene of interest 
		@param: gene_name: the name of the gene of interest 
		@note: if both gene Id and gene name are provided, both gene_id has a higher precedence 
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
	# define some magic functions 
	def __len__(self)->int:
		"""
		@brief: return the number of unique genes in the current instance
		"""
		return len(set(self._table.iloc[:,0]))
		
	def __str__(self)->str: 
		return f' A sub-cellular localization database covering {len(self)} genes'

	def __repr__(self)->str: 
		return str(self)


class GeneExpressionDB: 
	"""
	@brief: provides an API to access gene expression data stored in RNA consensus tissue gene data from Human proteome Atlas 
	@see: https://www.proteinatlas.org/about/download for more details 
	"""
	def __init__(self, path2data: str, 
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

class OrganismDB:
	"""
	@brief: extract information about the source organsim of a collection of protein sequences  
	from a fasta file and provides an API to query the results.  
	"""
	def __init__(self,path2Fasta:str)->OrganismDB:
		"""
		@param: path2fastaDB: The path to a fasta sequence database to obtain the protein sequences
		@note: The function expect the input fasta file to have header written in the FASTA format. 
		"""
		try: 
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
		"""
		@brief: get the number of unique organism in the database. 
		"""
		return list(set(self._map.values()))
	
	def get_number_protein_per_organism(self)->pd.DataFrame:
		"""
		@brief: provides a table containing the number of protein per organism. 
		"""
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
		"""
		@brief: return the parent organism of the provided proteins
		@param: prot_id: the provided protein id  
		"""
		return self._map[prot_id]
	
	def __len__(self)->int: 
		"""
		@brief: return the number of elements in the database
		"""
		return len(self._map)
	
	def __str__(self)->str: 
		"""
		@brief: return a string representation of the class 
		"""
		return str(f'An organism database with: {len(self)} entry.')

	def __repr__(self)->str: 
		return str(self)

		

	