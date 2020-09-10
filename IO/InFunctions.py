#!/usr/bin/env Python 
"""
@author: Hesham ElAbd
@brief: parse different inputs into user a standard format used by the library. 
@version: 0.0.1
"""
# load the models 
import pandas as pd 
import numpy as np 
from typing import Dict, List
import os  
from Bio import SeqIO
from Bio.PDB import PDBList
from pyteomics.mztab import MzTab
from pyteomics import pepxml, auxiliary
from pyteomics.openms import idxml
def load_identification_table(input_path: str, sep:str) -> pd.DataFrame:
    """
    @brief: load & process an identification table 
    @param: ident_table: the path two the identification table. with the following columns: peptides which hold the identified
	peptide sequence, protein which hold the identified protein sequence, start_index, and end_index where
	the last two columns define the position of the peptide in the parent protein. 
	@param: sep: The separator to parse the provided table. 
    """
	# load the table 
    try:
	    table: pd.DataFrame = pd.read_csv(input_path, sep=sep)
    except Exception as exp:
        raise IOError(f'Loading the input table from the provide path; {input_path} gave rise to the following IOException: {exp}')
    # assert that it has the correct shape 
    if table.shape[1]!= 4:
        raise ValueError(f'The provided table does not have the expected number of columns, Expected 4 columns, however, {table.shape} columns were found')
    # set the column names
    table.columns=['peptide','protein','start_index','end_index']
    # return the table 
    return table 

def parse_mzTab_to_identification_table(path2mzTab: str, path2fastaDB: str,
    fasta_reader_param: Dict[str,str]={'filter_decoy':True, 'decoy_string':'DECOY' })->pd.DataFrame:
    """
    @brief: parse a user provided mzTab to an identification table 
    @param: path2mz_tab: the path to the input mzTab file
    @param: path2fastaDB: the path to a fasta sequence database to obtain the protein sequences 
    @param: fasta_reader_param: a dict of parameters for controlling the behavior of the fasta reader 
    """
    # load the files 
    try: 
        input_file: MzTab =MzTab(path2mzTab)
    except Exception as exp:
        raise IOError(f'While parsing your input mzTab file: {path2mzTab}, the following error was encountered: {exp}')
    # load the fasta files and extract the accession of the proteins 
    try:
        sequence_dict: Dict[str,str]= fasta2dict(path2fastaDB, **fasta_reader_param)
    except Exception as exp:
        raise IOError(f'While parsing your input fasta file: {path2fastaDB}, the following error was encountered: {exp}')
    # get the peptide tables
    peptides_table: pd.DataFrame = input_file.peptide_table
    # construct the identification table 
    peptide_seq: List[str] = peptides_table.sequence.tolist()
    protein_acc: List[str]= [acc.split('|')[1] for acc in peptides_table.accession.tolist()]
    start_index: List[int] = []
    end_index: List[int] = []
    # fill extract the start and end-index information from the library 
    for idx in range(len(protein_acc)):
        # get the protein sequence 
        try:
            prot_seq: str = sequence_dict[protein_acc[idx]]
        except KeyError as exp: 
            raise KeyError(f'Database mismatch, the current protein accession: {protein_acc[idx]} is not defined in the provided sequence database')
        # get the index of the protein sequence
        try: 
            start_index.append(prot_seq.index(peptide_seq[idx]))
        except ValueError as exp:
            raise ValueError(f'Peptide sequence: {peptide_seq[idx]} could not be extracted from protein sequence: {prot_seq}')
        # add the end index 
        end_index.append(start_index[idx]+len(peptide_seq[idx]))
    # build the data frame 
    ident_table: pd.DataFrame = pd.DataFrame({
        'peptide': peptide_seq,
        'protein': protein_acc,
        'start_index':start_index,
        'end_index':end_index
    })
    # return the results 
    return ident_table
     
def parse_xml_based_format_to_identification_table(path2XML_file: str, path2fastaDB: str,
    decoy_prefix: str ='DECOY', is_idXML: bool = False, 
    fasta_reader_param: Dict[str,str]={'filter_decoy':True, 'decoy_string':'DECOY' }, 
    remove_if_not_matched: bool = True)->pd.DataFrame:
    """
    @brief: parse either a pepXML or an idXML file to generate an identification table , 
    @param: path2XML_file: The path to the input pepXML files
    @param: path2fastaDB: The path to a fasta sequence database to obtain the protein sequences
    @param: decoy_prefix: the prefix of the decoy sequences, default is DECOY
    @param: is_idXML: Whether or not the provided file is an idXML, default is false which assume the provided file is a pepXML file 
    @param: fasta_reader_param: a dict of parameters for controlling the behavior of the fasta reader 
    @param: remove_if_not_matched: remove the peptide if it could not be matched to the parent protein, default is True
    @note: see the function fasta2dict for more information regard the fasta reader 
    """
    # load the fasta files and extract the accession of the proteins 
    try:
        sequence_dict: Dict[str,str]= fasta2dict(path2fastaDB,**fasta_reader_param)
    except Exception as exp:
        raise IOError(f'While parsing your input fasta file: {path2fastaDB}, the following error was encountered: {exp}')
    # parse that the file exists:  
    if not os.path.exists(path2XML_file):
        raise ValueError(f'The provided path: {path2XML_file} does not exist!')
    # allocate a list to hold peptide and protein list 
    peptides: List[str] = []
    protein_acc: List[str] = []
    #  parse the XML file 
    if is_idXML: 
        with idxml.IDXML(path2XML_file) as reader:
            for elem in reader:
                for hit in elem['PeptideHit']:
                    for prot in hit['protein']: 
                        if decoy_prefix not in prot['accession']:
                            peptides.append(hit['sequence'])
                            protein_acc.append(prot['accession'].split('|')[1])
    else: 
        with pepxml.read(path2XML_file) as reader: 
            for elem in reader:
                for hit in elem['search_hit']:
                    for protein in hit['proteins']: 
                        if decoy_prefix not in protein['protein']:
                            peptides.append(hit['peptide'])
                            protein_acc.append(protein['protein'].split('|')[1])
    # extract the start and end index of the peptides from the parent proteins 
    start_index: List[int] = []
    end_index: List[int] = []
    # fill extract the start and end-index information from the library 
    for idx in range(len(protein_acc)):
        # get the protein sequence 
        try:
            prot_seq: str = sequence_dict[protein_acc[idx]]
        except KeyError as exp: 
            raise KeyError(f'Database mismatch, the current protein accession: {protein_acc[idx]} is not defined in the provided sequence database')
        # get the index of the protein sequence
        try: 
            start_index.append(prot_seq.index(peptides[idx]))
        except ValueError as exp:
            if remove_if_not_matched:
                start_index.append(-1) #  add a placeholder value that will be dropped later 
            else: 
                raise ValueError(f'Peptide sequence: {peptides[idx]} could not be extracted from protein sequence: {prot_seq}')
        # add the end index 
        end_index.append(start_index[idx]+len(peptides[idx]))
    # build the data frame 
    ident_table: pd.DataFrame = pd.DataFrame({
        'peptide': peptides,
        'protein': protein_acc,
        'start_index':start_index,
        'end_index':end_index
    })
    ident_table=ident_table.loc[ident_table.iloc[:,2]!=-1,:] # filter the non-matched peptides 
    # return the results 
    return ident_table


def parse_text_table(path2file: str, 
                    path2fastaDB: str,
                    sep=',',
                    fasta_reader_param: Dict[str,str]={'filter_decoy':True, 'decoy_string':'DECOY' }, 
                    seq_column: str = 'Sequence',
                    accession_column: str = 'Protein Accessions',
                    protein_group_sep: str = ';',
                    remove_protein_version: bool = True,
                    remove_if_not_matched: bool = True
            )->pd.DataFrame:
    """
    @brief: parse a user defined table to extract columns containing the identification table 
    @param: path2file: The path to load the CSV file holding the results 
    @param: path2fastaDB: The path to a fasta sequence database to obtain the protein sequences
    @param: sep: the table separators. Default is comma.
    @param: seq_column: the name of the columns containing the peptide sequence 
    @param: index_acc_column: the name of the column containing the protein accession and index column.
    @param: accession_column: the name of the column containing the protein accession 
    @param: fasta_reader_param: a dict of parameters for controlling the behavior of the fasta reader 
    @param: protein_group_sep: the separator for the protein group, default is semi colon
    @param: remove_protein_version: a bool if true strip the version number from the protein 
    @param: remove_if_not_matched: remove the peptide if it could not be matched to the parent protein, default is True
    @note: see the function fasta2dict for more information regard the fasta reader 
    """
    # load the fasta files and extract the accession of the proteins 
    try:
        sequence_dict: Dict[str,str]= fasta2dict(path2fastaDB,**fasta_reader_param)
    except Exception as exp:
        raise IOError(f'While parsing your input fasta file: {path2fastaDB}, the following error was encountered: {exp}')
    # Load the table 
    try: 
        input_table: pd.DataFrame = pd.read_csv(path2file,sep=sep)
    except Exception as exp:
        raise IOError(f'while loading the input table: {path2file}, The following error was encountered: {exp}')
    # extract the columns of the table 
    if seq_column not in input_table.columns or accession_column not in input_table.columns :
        raise KeyError(f'The provided names for the peptides sequence: {seq_column} and/or the indexing column: {index_acc_column} and/or accession column: {accession_column} could not be found on the table')
    # allocate the lists
    peptides: List[str] = []
    protein_acc: List[str] = []
    start_index: List[int] = []
    end_index: List[int] = []
    # loop of the table to fill the lists 
    for _, row in input_table.iterrows():
        # check if the peptide is present in more than one protein 
        if protein_group_sep in row[accession_column]:
            accessions: List[str] =row[accession_column].split(protein_group_sep)
            for accession in accessions:
                # extract the macting protein sequence 
                accession=accession.strip(' ')
                try: 
                    protein_seq: str = sequence_dict[accession]
                except KeyError as exp:
                    if remove_if_not_matched:
                        accession=accession.split('-')[0]
                        try:
                            protein_seq= sequence_dict[accession]
                        except KeyError as exp: 
                            raise KeyError(f"Database mismatch, the current protein accession: {accession.strip(' ')} is not defined in the provided sequence database")
                    else: 
                        raise KeyError(f"Database mismatch, the current protein accession: {accession.strip(' ')} is not defined in the provided sequence database")
                peptides.append(row[seq_column])
                protein_acc.append(accession)
                # get the index 
                try: 
                    temp_start_idx: int =protein_seq.index(row[seq_column])
                except ValueError as exp:
                    if remove_if_not_matched: 
                        temp_start_idx=-1
                        pass
                    else: 
                        raise ValueError(f'Peptide sequence: {row[seq_column]} could not be extracted from protein sequence: {prot_seq}')
                start_index.append(temp_start_idx)
                end_index.append(temp_start_idx+len(row[seq_column]))
        else: 
            accession: int = row[accession_column].strip(' ')
            try: 
                protein_seq: str = sequence_dict[accession]
            except KeyError as exp:
                if remove_if_not_matched:
                    accession=accession.split('-')[0]
                    try:
                        protein_seq= sequence_dict[accession]
                    except KeyError as exp: 
                        raise KeyError(f"Database mismatch, the current protein accession: {accession} is not defined in the provided sequence database")
                else: 
                    raise KeyError(f"Database mismatch, the current protein accession: {accession} is not defined in the provided sequence database")
            peptides.append(row[seq_column])
            protein_acc.append(accession)
            # get the index 
            try: 
                temp_start_idx: int =protein_seq.index(row[seq_column])
            except ValueError as exp:
                if remove_if_not_matched: 
                    temp_start_idx=-1
                    pass
                else: 
                    raise ValueError(f'Peptide sequence: {row[seq_column]} could not be extracted from protein sequence: {prot_seq}')
            start_index.append(temp_start_idx)
            end_index.append(temp_start_idx+len(row[seq_column]))            
    # build the data frame 
    ident_table: pd.DataFrame = pd.DataFrame({
        'peptide': peptides,
        'protein': protein_acc,
        'start_index':start_index,
        'end_index':end_index
    })
    ident_table=ident_table.loc[ident_table.iloc[:,2]!=-1,:] # filter the non-matched peptides 
    # return the results 
    return ident_table 

def fasta2dict(path2fasta:str, filter_decoy: bool = True,
            decoy_string: str = 'DECOY')->Dict[str,str]: 
    """
    @brief: loads a fasta file and construct a dict object where ids are keys and sequences are the value
    @param: path2fasta: the path to load the fasta file 
    @param: filter_decoy: a boolean of whether or not to filter the decoy sequences from the database, default is True
    @param: decoy_string: the decoy database prefix, only valid incase filter_decoy is set to true, default is DECOY
    """
    # check the path exists 
    if not os.path.exists(path2fasta):
        raise IOError(f"The provided path to fasta: {path2fasta} does not exits!")
    # load the fasta file 
    seq_gen: SeqIO.FastaIO.FastaIterator = SeqIO.parse(path2fasta,'fasta')
    results=dict()
    for seq in seq_gen: 
        if filter_decoy:
            if decoy_string in seq.id:
               continue
            else:
                results[seq.id.split('|')[1]]=str(seq.seq)
        else:
            results[seq.id.split('|')[1]]=str(seq.seq)  
    # return the results 
    return results

def download_pdb_entry(prot_id: str) ->str: 
    """
    @brief: download the structure of a protein from protein databank form as mmCIF file 
    """
    pdb_res=PDBList()
    try: 
        return pdb_res.retrieve_pdb_file(prot_id)
    except Exception as exp: 
        raise IOError(f'While Downloading the structure of protein: {prot_id}, the following error was encountered: {exp}')

