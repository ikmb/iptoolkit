#!/usr/bin/env python
"""An interface to query and retrive data from an database designed according to the builtin database, see IPTK.Utils.CreateDB for more details 
"""
## LOAD THE MODULES 
#------------------
from __future__ import annotations
import pandas as pd 
import os
from typing import List
from datetime import datetime
import mhcnames
## DEFINE THe CLASS
#------------------
class BuildInDB:
    def __init__(self,path2file:str, sep:str='\t')->BuildInDB:
        """Loads a pandas dataframe that comply with the expected structure

        :param path2file: The path to load the CSV files. 
        :type path2file: string. 
        :param sep: The seperator for the table, defaults to tab. 
        :type sep: string. 
        :return: a buildIn database instance that can be used to handle the sequences in the database
        :rtype: BuildInDB
        """
        try:
            self._database=pd.read_csv(path2file,sep=sep)
        except Exception as exp:
            raise IOError(f'while reading the database located in: {path2file} the following error was encountered: {str(exp)}')
        try:
            self._validate_database(self._database)
        except Exception as exp:
            raise ValueError(f'The provided data is not a valid database, because: {str(exp)}')
    
    def _validate_database(self, database:pd.DataFrame)->None:
        """ validate that a provided database has the correct number, name and order of columns.

        :param database: a proposed database to validate 
        :type database: pd.DataFrame
        :raises ValueError: incase validating the database failed 
        """
        if database.shape[1]!=11:
            raise ValueError(f'Expected the provided table to have 11 eleven columns, however the provided table has : {self._database.shape[1]}')
        def_col_name=['experiment_id','peptide','protein','start_pos','end_pos','organism','allele','source','bioSample','antibody','search_engine']
        if database.columns!=def_col_name:
            raise ValueError(f"Expected the provided table to have the exact order and name of columns as follow, {', '.join(def_col_name)}, however the \
                the current table has {', '.join(database.columns)}.")
        return
    
    def __str__(self)-> str:
        return f'A build in database with {self._database.shape[0]} entries.'

    def __repr__(self) -> str: 
        return str(self)
    
    def __len__(self)->int:
        return self._database.shape[0]

    def append_to_data_base(self,content_to_add: pd.DataFrame)->None:
        """[summary]

        :param content_to_add: [description]
        :type content_to_add: pd.DataFrame
        """
        try:
            self._validate_database(content_to_add) 
        except ValueError as exp:
            raise ValueError(f'Validating your provided table because: {str(exp)}')
        self._database=pd.concat([self._database,content_to_add], axis=0)
        return 

    def write_to_disk(self, path:str, add_timestampe:bool = True)->None:
        """write the generated database to the disk, the default saperator used is tab

        :param path: the path, INCLUDING the name of the genrated database 
        :type path: str
        :param add_timestampe: wether or not to a timestamp to the end of the file_name, defaults to True
        :type add_timestampe: bool, optional
        """
        if add_timestampe:
            time_stamp=str(datetime.now()).replace(' ','_')  
            path=path.split('.')[0]
            path+=time_stamp
        try:
            self._database.to_csv(path+'.tsv',sep='\t',index=False)
        except Exception as exp:
            raise IOError(f'While writing the database the following error was encountered: {exp}')
        return

    def get_peptides_from_experiment(self,experiment_id:str)->List[str]:
        """returns a list of peptides that were identified from a specific experiment

        :param experiment_id: the experiment id 
        :type experiment_id: str
        :raises ValueError: if the experiment id is not defined in the database 
        :return: a list of peptides that were identified from a specific experiment
        :rtype: List[str]
        """
        if experiment_id not in set(list(self._database.iloc[:,0])):
            raise ValueError(f'The provided id is not in the database')
        return self._database.loc[self._database.iloc[:,0]==experiment_id].iloc[:,1].tolist() 

    def get_proteins_from_experiment(self,experiment_id:str)->List[str]:
        """returns a list of proteins that were identified from a specific experiment.

        :param experiment_id: the experiment id 
        :type experiment_id: str
        :raises ValueError: if the experiment id is not defined in the database 
        :return: a list of proteins that were identified from a specific experiment
        :rtype: List[str]
        """
        if experiment_id not in set(list(self._database.iloc[:,0])):
            raise ValueError(f'The provided id is not in the database')
        return self._database.loc[self._database.iloc[:,0]==experiment_id].iloc[:,2].tolist() 

    def get_table_from_experiment(self,experiment_id:str)->pd.DataFrame:
        """returns a subtable contain all info associated with the provided experimental-id 

        :param experiment_id: the experiment id 
        :type experiment_id: [type]
        :raises ValueError: [description]
        :return: a table contain all the experiments associated with the provided experimental-id 
        :rtype: pd.DataFrame
        """
        if experiment_id not in set(list(self._database.iloc[:,0])):
            raise ValueError(f'The provided id is not in the database')
        return self._database.loc[self._database.iloc[:,0]==experiment_id]

    def get_peptides_bound_by_allele(self,allele_name:str)->pd.DataFrame:
        """ returns a subtable contain all info associated with the provided allele_name

        :param allele_name: the allele name written in the standard notation 
        :type allele_name: str
        :raises ValueError: incase the provided allele is not in the database 
        :return: a subtable contain all the experiments associated with the provided allele_name
        :rtype: pd.DataFrame
        """
        s_allele_name=mhcnames.parse_allele_name(allele_name)
        s_allele_name=s_allele_name.gene+'*'+s_allele_name.allele_family+':'+s_allele_name.allele_code
        current_alleles=self.get_all_unique_alleles_in_db()
        if s_allele_name not in current_alleles:
            raise ValueError(f"The provided allele is not currently defined in the database, currently the following is defined: {' ,'.join(current_alleles)}")
        return  self._database.loc[self._database.allele.contains(s_allele_name),]
    
    def get_all_unique_alleles_in_db(self)->List[str]:
        """ returns a list of all unique alleles in the database
        :return: a list of all unique alleles in the database
        :rtype: List[str]
        """
        allele_pairs=list(set([alleles for alleles in self._database.allele]))
        alleles = []
        for allele_pair in allele_pairs:
            alleles.extend(allele_pair.extend(';')) 
        return list(set(alleles))
    
    def get_all_unique_proteins_in_db(self)->List[str]:
        """ return a list of all unique peptides in the database 

        :return: all unique peptides in the database 
        :rtype: List[str]
        """
        return list(set(self._database.iloc[:,1]))


    def get_all_unique_peptides_in_db(self)->List[str]:
        """ return a list of all unique proteins in the database 

        :return: all unique proteins in the database 
        :rtype: List[str]
        """
        return list(set(self._database.iloc[:,2]))
 
    def get_all_unique_organims_in_db(self)->List[str]:
        """ return a list of all unique proteins in the database 

        :return: all unique proteins in the database 
        :rtype: List[str]
        """
        return list(set(self._database.iloc[:,5]))

    def get_organism_info_from_db(self, org_name:str)->pd.DataFrame:
        """returns a sub-table contain all info associated with the provided org_name

        :param org_name: the name of the organism to query the database 
        :type org_name: string 
        :return: a sub-table contain all info associated with the provided org_name
        :rtype: pd.DataFrame
        """
        current_organims=self.get_all_unique_organims_in_db()
        if org_name not in current_organims:
            raise ValueError(f'The provided organims is not in the database: {org_name}')
        return self._database.loc[self._database.organism==org_name,]

    def get_peptide_info_from_db(self,peptide_seq)->pd.DataFrame:
        """returns a sub-table contain all info associated with the provided peptide

        :param peptide_seq: the peptide sequence to query the database with 
        :type peptide_seq: string 
        :return: a sub-table contain all info associated with the provided org_name
        :rtype: pd.DataFrame
        """ 
        return self._database.loc[self._database.peptide==peptide_seq,]

    def get_protein_info_from_db(self,protein_id)->pd.DataFrame:
        """returns a sub-table contain all info associated with the provided protein id

        :param protein_id: protein id
        :type peptide_seq: string 
        :return: a sub-table contain all info associated with the provided protein id 
        :rtype: pd.DataFrame
        """ 
        return self._database.loc[self._database.protein==protein_id,] 
        
    def get_table(self)->pd.DataFrame:
        return self._database
    
    