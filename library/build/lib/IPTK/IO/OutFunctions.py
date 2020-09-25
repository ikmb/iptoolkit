#!/usr/bin/env python 
"""Write the results generated by the library into a wide variety of formats. 
"""
import numpy as np 
import h5py
from IPTK.Classes.Experiment import Experiment
from IPTK.Classes.Peptide import Peptide
from typing import List
import pandas as pd
# define some types 
Peptides= List[Peptide]
Names= List[str]
# define the function of the module 
def write_mapped_tensor_to_h5py(tensor: np.ndarray, path2write: str, dataSet_name: str ='MAPPED_TENSOR')->None:
    """ Write a mapped tensor to an hdf5 file 

    :param tensor: The provided tensor to  write it to the hdf5 file. 
    :type tensor: np.ndarray
    :param path2write: The path of the output file 
    :type path2write: str
    :param dataSet_name: The name of the dataset inside the mapped tensor, defaults to 'MAPPED_TENSOR'
    :type dataSet_name: str, optional
    :raises IOError: In case opening the file for writing failed 
    """
    try:
        file_handler=h5py.File(path2write,'w')
    except Exception as exp:
        raise IOError(f'While opening the file: {path2write} for writing, the following error was encountered: {exp}')
    # create the dataset 
    file_handler.create_dataset(dataSet_name,data=tensor)
    # flush the data 
    file_handler.close()
    return 

def write_annotated_sequences(peptides: List[str], labels: List[int], path2write: str,
                            sep: str = ',', shuffle: bool =True) ->None:
    """take a list of peptides along with it sequences and write the results to a CSV file.
    
    :param peptides: a list of peptide sequences 
    :type peptides: List[str]
    :param labels: a list of numerical labels associated with the peptides 
    :type labels: List[int]
    :param path2write: the path to write the generated file 
    :type path2write: str
    :param sep: The separator in the resulting table,  defaults to ','
    :type sep: str, optional
    :param shuffle: Whether or not to shuffle the table , defaults to True
    :type shuffle: bool, optional
    :raises ValueError: incase the length of the tables and labels is not matching
    :raises IOError: In case writing the output table failed
    """
    if len(peptides) != len(labels):
        raise ValueError(f'The provided lists do not match, peptides have: {len(peptides)} elements while labels have: {len(peptides)} elements')
    # define the dataframe results 
    table: pd.DataFrame = pd.DataFrame({
        'peptides':peptides,
        'labels': labels
    }) 
    # shuffle the database 
    if shuffle:
        table=table.sample(frac=1.0)
    # write the results 
    try: 
        table.to_csv(path2write,sep=sep, index=False, header=True)
    except Exception as exp: 
        raise IOError(f'While writing the results table, the following error was encountered: {exp}')
    

def write_named_peptides_to_fasta(names: Names, peptides: Peptides, output_file:str):
    """Takes a list of names and peptide sequences and writes them as an output file to the disk as fasta files

    :param names: A list of sequences names 
    :type names: Names
    :param peptides: A list of peptide sequences 
    :type peptides: Peptides
    :param output_file: The name of the output file to write the results to, it will OVERWRITE existing files
    :type output_file: str
    :raises ValueError: Incase the length of the tables and labels is not matching
    :raises IOError: In case writing the output file failed
    """
    if len(names) != len (peptides): 
        raise ValueError(f"The length of the names container: {len(names)} does not equal the length of the sequences container: {len(peptides)}")
    # write the files 
    try:
        with open(output_file, 'w') as writer:
            for idx in range(len(names)): 
                writer.write('>'+names[idx]+'\n')
                writer.write(peptides[idx]+'\n')
    except Exception as exp: 
        raise IOError(f'While writing the output file the following error was encountered: {exp}') 

def write_auto_named_peptide_to_fasta(peptides:Peptides, output_file:str )->None:
    """ Takes a list of peptides, generate automatic names for the peptides and write the results to the disk as FASTA files 

    :param peptides: a list of peptide sequences 
    :type peptides: Peptides
    :param output_file: the name of the output file to write the results to, it will OVERWRITE existing files
    :type output_file: str
    """
    # generate automatic names for the pepetide 
    names=['PEP_AUTO_NAME_'+str(idx) for idx in range(len(peptides))]
    # call the writer functions 
    write_named_peptides_to_fasta(names, peptides,output_file)
    
def write_pep_file(peptides:Peptides, output_file:str):
    """Takes a file and write the peptides to .pep file which is a text file that contain one peptide per line
    
    :param peptides: a list of peptide sequences 
    :type peptides: Peptides
    :param output_file: the name of the output file to write the results to, it will OVERWRITE existing files
    :type output_file: str
    :raises IOError: In case writing the output file failed
    """
    try: 
        with open(output_file,'w') as writer: 
            for pep in peptides: 
                writer.write(pep+'\n')
    except Exception as exp: 
        raise IOError(f'While writing the file, the following error was encountered: {exp}')
    