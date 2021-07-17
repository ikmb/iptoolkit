#!/usr/bin/env Python 
"""Parse different user inputs into a standard format/tables used by the library. 
"""
# load the models 
import pandas as pd 
import numpy as np 
from typing import Dict, List
import os  
from Bio import SeqIO
from Bio.PDB import PDBList
from pyteomics.mztab import MzTab
from pyteomics import pepxml
from pyteomics.openms import idxml
from pyteomics.mzid import MzIdentML
from tqdm import tqdm 
import time 

# define the functions of the modules
def load_identification_table(input_path: str, sep:str) -> pd.DataFrame:
    """load & process an identification table 
    
    :param input_path: the path to the identification table. With the following columns: peptides which hold the identified
	peptide sequence, protein which hold the identified protein sequence, start_index, and end_index where
	the last two columns define the position of the peptide in the parent protein. 
    :type input_path: str
    :param sep: The separator to parse the provided table. 
    :type sep: str
    :rtype: pd.DataFrame
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
    fasta_reader_param: Dict[str,str]={'filter_decoy':True, 'decoy_string':'DECOY' },
    remove_if_not_matched: bool = True)->pd.DataFrame:
    """parse a user provided mzTab to an identification table 

    :param path2mzTab: the path to the input mzTab file
    :type path2mzTab: str
    :param path2fastaDB: the path to a fasta file to obtain the protein sequences 
    :type path2fastaDB: str
    :param fasta_reader_param: A dict of parameters for controlling the behavior of the fasta reader , defaults to {'filter_decoy':True, 'decoy_string':'DECOY' }
    :type fasta_reader_param: Dict[str,str], optional
    :param remove_if_not_matched: remove the peptide if it could not be matched to the parent protein, defaults to True
    :type remove_if_not_matched: bool, optional
    :raises IOError: if the mztab file could not be open and loaded or if the fasta database could not be read
    :raises KeyError: if a protein id defined in the mzTab file could not be extracted from a matched sequence database
    :raises ValueError: if the peptide can not be mapped to the identified protein 
    :return: the identification table 
    :rtype: pd.DataFrame
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
    peptides: List[str] = peptides_table.sequence.tolist()
    protein_acc: List[str]= [acc.split('|')[1] for acc in peptides_table.accession.tolist()]
    start_index: List[int] = []
    end_index: List[int] = []
    # fill extract the start and end-index information from the library 
    print(f"Parsing the provided mzTab table ... started at: {time.ctime()}")
    for idx in tqdm(range(len(protein_acc))):
        # get the protein sequence 
        try:
            prot_seq: str = sequence_dict[protein_acc[idx]]
        except KeyError as exp: 
            raise KeyError(f'Database mismatch, the current protein accession: {protein_acc[idx]} is not defined in the provided sequence database')
       # get the index of the protein sequence
        try: 
            if '(' in peptides[idx]: # that is there sequence modifications in the sequence 
                temp_peptide=peptides[idx] # that is there sequence modifications in the sequence 
                while '(' in temp_peptide or ')' in temp_peptide: 
                    pre_seq=temp_peptide.split('(')[0]
                    post_seq=")".join(temp_peptide.split(')')[1:])
                    temp_peptide=pre_seq+post_seq
                start_index.append(prot_seq.index(temp_peptide))
                peptide_len=len(temp_peptide)
            else: 
                start_index.append(prot_seq.index(peptides[idx]))
                peptide_len=len(peptides[idx])
        except ValueError as exp:
            if remove_if_not_matched:
                start_index.append(-1) #  add a placeholder value that will be dropped later 
            else: 
                raise ValueError(f'Peptide sequence: {peptides[idx]} could not be extracted from protein sequence: {prot_seq} with accession: {protein_acc[idx]}')
        # add the end index 
        end_index.append(start_index[idx]+peptide_len)
    # build the data frame 
    ident_table: pd.DataFrame = pd.DataFrame({
        'peptide': peptides,
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
    """parse either a pepXML or an idXML file to generate an identification table , 

    :param path2XML_file: The path to the input pepXML files
    :type path2XML_file: str
    :param path2fastaDB: The path to a fasta sequence database to obtain the protein sequences
    :type path2fastaDB: str
    :param decoy_prefix: the prefix of the decoy sequences, default is DECOY
    :type decoy_prefix: str, optional
    :param is_idXML: Whether or not the provided file is an idXML, default is false which assume the provided file is a pepXML file, defaults to False
    :type is_idXML: bool, optional
    :param fasta_reader_param: A dict of parameters for controlling the behavior of the fasta reader, defaults to {'filter_decoy':True, 'decoy_string':'DECOY' }
    :type fasta_reader_param: Dict[str,str], optional
    :param remove_if_not_matched: remove the peptide if it could not be matched to the parent protein, defaults to True
    :type remove_if_not_matched: bool, optional
    :raises IOError: if the fasta database could not be open 
    :raises ValueError: if the XML file can not be open 
    :raises KeyError: if a protein id defined in the mzTab file could not be extracted from a matched sequence database
    :raises ValueError: if the peptide can not be mapped to the identified protein 
    :return: the identification table 
    :rtype: pd.DataFrame
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
            print(f"Parsing the provided idXML table ..., started at: {time.ctime()}")
            for elem in tqdm(reader):
                for hit in elem['PeptideHit']:
                    for prot in hit['protein']: 
                        if decoy_prefix not in prot['accession']:
                            peptides.append(hit['sequence'])
                            if '|' in prot['accession']: 
                                protein_acc.append(prot['accession'].split('|')[1])
                            else: 
                                protein_acc.append(prot['accession'])
    else: 
        with pepxml.read(path2XML_file) as reader:
            print(f"Parsing the provided pepXML table ..., started at: {time.ctime()}") 
            for elem in tqdm(reader):
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
            if '(' in peptides[idx]: # that is there sequence modifications in the sequence 
                temp_peptide=peptides[idx] # that is there sequence modifications in the sequence 
                while '(' in temp_peptide or ')' in temp_peptide: 
                    pre_seq=temp_peptide.split('(')[0]
                    post_seq=")".join(temp_peptide.split(')')[1:])
                    temp_peptide=pre_seq+post_seq
                start_index.append(prot_seq.index(temp_peptide))
                peptide_len=len(temp_peptide)
            else: 
                start_index.append(prot_seq.index(peptides[idx]))
                peptide_len=len(peptides[idx])
        except ValueError as exp:
            if remove_if_not_matched:
                start_index.append(-1) #  add a placeholder value that will be dropped later 
            else: 
                raise ValueError(f'Peptide sequence: {peptides[idx]} could not be extracted from protein sequence: {prot_seq} with accession: {protein_acc[idx]}')
        # add the end index 
        end_index.append(start_index[idx]+peptide_len)
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
    """Parse a user defined table to extract the columns containing the identification table 
    
    :param path2file: The path to load the CSV file holding the results 
    :type path2file: str
    :param path2fastaDB: The path to a fasta sequence database to obtain the protein sequences
    :type path2fastaDB: str
    :param sep: The table separators, defaults to ','
    :type sep: str, optional
    :param fasta_reader_param: a dict of parameters for controlling the behavior of the fasta reader, defaults to {'filter_decoy':True, 'decoy_string':'DECOY' }
    :type fasta_reader_param: Dict[str,str], optional
    :param seq_column: The name of the columns containing the peptide sequence, defaults to 'Sequence'
    :type seq_column: str, optional
    :param accession_column: The name of the column containing the protein accession, defaults to 'Protein Accessions'
    :type accession_column: str, optional
    :param protein_group_sep: The separator for the protein group, defaults to ';'
    :type protein_group_sep: str, optional
    :param remove_protein_version: A bool if true strip the version number from the protein , defaults to True
    :type remove_protein_version: bool, optional
    :param remove_if_not_matched: remove the peptide if it could not be matched to the parent protein, defaults to True
    :type remove_if_not_matched: bool, optional
    :raises IOError: Incase either the sequences database or the identification table can not be open and loaded
    :raises KeyError: In case the provided column names not in the provided identification table. 
    :raises KeyError: Incase the protein sequence can not be extract from the sequence database 
    :raises ValueError: incase the peptide could not be located in the protein sequence 
    :return: an identification table
    :rtype: pd.DataFrame
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
        raise KeyError(f'The provided names for the peptides sequence: {seq_column} and/or the indexing column: {accession_column} and/or accession column: {accession_column} could not be found on the table')
    # allocate the lists
    peptides: List[str] = []
    protein_acc: List[str] = []
    start_index: List[int] = []
    end_index: List[int] = []
    # loop of the table to fill the lists 
    print(f"Parsing the provided table ... started at: {time.ctime()}")
    for _, row in tqdm(input_table.iterrows()):
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
                        raise ValueError(f'Peptide sequence: {row[seq_column]} could not be extracted from protein sequence: {protein_seq}')
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
                    raise ValueError(f'Peptide sequence: {row[seq_column]} could not be extracted from protein sequence: {protein_seq}')
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

def parse_mzIdentML_to_identification_table(file2load:MzIdentML)->pd.DataFrame:
    """Parse an input mzIdentML file and parse its content and return an identification table 

    Args:
        file2load (MzIdentML): The path to load the MzIdentML file 

    Raises:
        IOError: [description]

    Returns:
        pd.DataFrame: an identification table 
    """
    if not os.path.exists(file2load):
        raise IOError(f"The provided path to load the file: {file2load} does not exists")
    ## create a list to hold the results 
    temp_seq, accession, start_idx, end_idx=[],[],[],[]
    ## loop over all the array
    try:
        print(f"Parsing the input MzIdentML ..., started at: {time.ctime()}")
        for rec in tqdm(MzIdentML(file2load)):
            for item in rec['SpectrumIdentificationItem']:
                if item['passThreshold']==True:
                    evidences=item['PeptideEvidenceRef']
                    for evidence in evidences: 
                        if evidence['isDecoy']==False:
                          temp_seq.append(evidence['PeptideSequence'])
                          accession.append(evidence['accession'])
                          start_idx.append(evidence['start'])
                          end_idx.append(evidence['end'])        
    except Exception as exp:
        raise IOError(f"While Parsing the content of the file, the following error was encounterred: {str(exp)}")
    ## remove peptide sequence
    peptide_seq=[]
    for pep in temp_seq:
        if '(' in pep:
            loop=True
            while(loop):
                start_pos=pep.find('(')
                end_pos=pep.find(')')
                pep=pep[:start_pos]+pep[end_pos+1:]
                if '(' not in pep:
                    loop=False
        else:
            peptide_seq.append(pep)
    ## create a pandas dataframe 
    return pd.DataFrame({
        'peptide': peptide_seq,
        'protein': accession,
        'start_index':start_idx,
        'end_index':end_idx
    })

def fasta2dict(path2fasta:str, filter_decoy: bool = True,
            decoy_string: str = 'DECOY')->Dict[str,str]: 
    """loads a fasta file and construct a dict object where ids are keys and sequences are the value

    :param path2fasta: The path to load the fasta file 
    :type path2fasta: str
    :param filter_decoy: A boolean of whether or not to filter the decoy sequences from the database, defaults to True
    :type filter_decoy: bool, optional
    :param decoy_string: The decoy database prefix, only valid incase filter_decoy is set to true, defaults to 'DECOY'
    :type decoy_string: str, optional
    :raises IOError: [description]
    :return: A dict where the protein ids are the keys and the protein sequences are the values 
    :rtype: Dict[str,str]
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
                if '|' in seq.id: 
                    results[seq.id.split('|')[1]]=str(seq.seq)
                else: 
                    results[seq.id]=str(seq.seq)
        else:
            if '|' in seq.id: 
                results[seq.id.split('|')[1]]=str(seq.seq)
            else: 
                results[seq.id]=str(seq.seq)
    # return the results 
    return results

def download_pdb_entry(prot_id: str) ->str: 
    """ Download the structure of a protein from protein databank as an mmCIF file. 

    :param prot_id: the protein id 
    :type prot_id: str
    :raises IOError: incase downloading and accessing the data failed 
    :return: the path to the downloaded file 
    :rtype: str
    """
    pdb_res=PDBList()
    try: 
        return pdb_res.retrieve_pdb_file(prot_id)
    except Exception as exp: 
        raise IOError(f'While Downloading the structure of protein: {prot_id}, the following error was encountered: {exp}')

