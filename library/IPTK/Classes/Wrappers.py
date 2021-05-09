"""A collection of convient classes to construct different IPTK classes 
"""
## load the modules and import other classes 
from __future__ import annotations
import os
from typing import List, Dict, Union 
from IPTK.Classes.Experiment import Experiment
from IPTK.Classes.Database import SeqDB, GeneExpressionDB, CellularLocationDB,OrganismDB
from IPTK.Classes.Tissue import Tissue
from IPTK.Classes.HLASet import HLASet
from IPTK.Classes.Proband import Proband
from IPTK.Classes.Peptide import Peptide 
from IPTK.IO.InFunctions import (parse_mzTab_to_identification_table, 
parse_xml_based_format_to_identification_table, parse_text_table) 
from tqdm import tqdm 

## Define the experiment class 
class cExperimet:
    def __init__(self,filepath:str, path2fasta:str, fileformat:str='idXML',tissue_name:str='total PMBC',
                proband_name:str='Default Proband', hla_set:List[str]=['DRB1*15:01','DRB1*15:01'])->cExperimet:
        """A Wrapper class for constracting an experimental dataset using user defined parameters\
            The class take care of initializing all classes and functions provided an easy-to-use interface\
            for working with immunopeptidomics data 
        Args:
            filepath (str): the path to load the input file, for example and idXML or an Identification table 
            path2fasta (str): the path to load Fasta database 
            fileformat (str, optional): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table (IdTable).\
                 Defaults to 'idXML'.
            tissue_name (str, optional): The name of the tissue to utilize, this is used for initializing the gene expression table\
                 Defaults to 'total PMBC'.
            proband_name (str, optional): the name of the proband from whome the data was obtained. Defaults to 'Default Proband'.
            hla_set (List[str], optional): A list of HLA alleles from whome the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].

        Returns:
            cExperimet: an IPTK.Class.Wrapper.Experiment class, an IPTK.Class.Experiment.Experiment can be extracted from the resutned instance using the get_experiment method 
        """
        ## Checking that the input is correct 
        if not os.path.exists(filepath):
            raise ValueError(f"The provided path for the identification file : {filepath} does not exist!!")
        if not os.path.exists(path2fasta):
            raise ValueError(f"The path to the proivded fasta file: {path2fasta}, does not exists!!!")
        if fileformat not in ['idXML', 'pepXML', 'IdTable','mzTab']:
            raise ValueError(f"Unknow input format, the provided format: {fileformat} is not supported, currently supported values are: {', '.join(['idXML', 'pepXML', 'IdTable','mzTab'])}")
        # define the data 
        self._proband = Proband(name=proband_name)  # the name of the proband 
        try: 
            self._hLASet = HLASet(hlas=['HLA-DRB1*15:01']) # just a place holder to represent the HLA allele, an instance of class HLASet
        except Exception as exp:
            raise RuntimeError(f"The following error was Encountered while creating an HLASet: \n{str(exp)}\n")
        try: 
            self._seqBD= SeqDB(path2fasta=path2fasta)
        except Exception as exp:
            raise IOError(f"While loading the fasta database the following error was Encountered : \n{str(exp)}\n")
        self._expresson_profile = GeneExpressionDB() # use the data on the human protein atlas @https://www.proteinatlas.org/about/download --> Normal tissue data 
        self._protein_locations = CellularLocationDB() # use the data on the human protein atlas @https://www.proteinatlas.org/about/download --> Subcellular location data
        try:
            self._tissue: Tissue(name='small intestine',main_exp_value=expresson_profile, 
                        main_location=protein_locations) # create the tissue instance
        except Exception as exp:
            raise RuntimeError(f"While creating a tissue instance, the following error was Encountered: \n{str(exp)}\n")
        try:
            if fileformat == 'idXML':
                input_table = parse_xml_based_format_to_identification_table(path2XML_file= filepath, path2fastaDB=path2fasta, is_idXML= True)
            elif fileformat == 'pepXML':
                input_table = parse_xml_based_format_to_identification_table(path2XML_file= filepath, path2fastaDB=path2fasta, is_idXML= False) 
            elif fileformat == 'mzTab':
                input_table = parse_mzTab_to_identification_table(path2mzTab= filepath, path2fastaDB=path2fasta)
            else:
                input_table = parse_text_table(filepath,path2fasta)
        except Exception as exp:
            raise ValueError(f"Loading the input table has caused to the following error: \n{str(exp)}\n")
        # constructing the experiments 
        try: 
            self._exp=Experiment(self._proband, self._hLASet, self._tissue, self._seqBD, ident_table)
        except:
            raise RuntimeError(f"Generating an Experiment instance has caused to the following error: \n{str(exp)}\n")
        self._cashed_results=dict()
        return 
    ##  
    def get_experiment(self)->Experiment:
        """return the instance Experiment

        Returns:
            Experiment: an IPTK.Class.Experiment.Experiment that can be used for analyzing the data
        """
        return self._exp 
    def cash_results(self, results_name:str, results_value):
        """add a pair of results to the current instance hashmap  

        Args:
            results_name (str): the name of the result to add 
            results_value ([type]): the result value to add 
        """
        self._cashed_results[results_name]=results_value

class PairedExperiment:
    def __init__(self,path:Union[str,os.path], anchor_name:str ,path2fasta:str,
        fileformat:str='idXML', tissue_name:str='total PMBC', proband_name:str='Default Proband',
        hla_set:List[str]=['DRB1*15:01','DRB1*15:01'])->PairedExperiment:
        """A Wrapper around an IPTK.Classes.Wrapper.Experiment instance, providing a consise, easy-to-use,  
        interface to compare different replicates. 

        Args:
            path (Union[str,os.path]): The path to the root directory where both files are defined 
            anchor_name (str): The name shared among all replicates 
            path2fasta (str): the path to load Fasta database 
            fileformat (str, optional): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table (IdTable).\
                 Defaults to 'idXML'.
            tissue_name (str, optional): The name of the tissue to utilize, this is used for initializing the gene expression table\
                 Defaults to 'total PMBC'.
            proband_name (str, optional): the name of the proband from whome the data was obtained. Defaults to 'Default Proband'.
            hla_set (List[str], optional): A list of HLA alleles from whome the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].

        Returns:
            PairedExperiment: [description]
        """
        ## Check that the provided path exist 
        if not os.path.exists(path):
            raise ValueError(f"The provided root directory: {path}, does not exist")
        ## get path to samples replicates 
        replicates=[os.path.join(path, name) for name in os.listdir(path) if anchor_name in name]
        ## check the correct number of replicates 
        if len(replicates)==0:
            raise ValueError(f"Searching using the provided path: {path} with the provided anchor name: {anchor_name} returned zero files")
        elif len(replicates)==1:
            raise ValueError(f"Searching using the provided path: {path} with the provided anchor name: {anchor_name} returned only one file: {replicates[0]}!")
        ## define a dict to hold the replicates 
        exps=dict()
        try:
            for file_name in tqdm(replicates):
                parsed_name=file_name.split('/')[-1].split('.')[0]
                experiment=cExperimet(file_name, path2fasta, fileformat, tissue_name, proband_name, hla_set)
                exps[parsed_name]=experiment
        except Exception as exp:
            raise RuntimeError(f"Creating an Experimental instance from the file: {file_name} fialed with the following error \n{str(exp)}\n")
        self._exps=ExperimentSet(exps)
        self._cashed_results=dict()
        return  
    
    def get_repro_rate(self, level='protein')->float:
        """Compute the reproducibility rate among the replicates. 
        Args:
            level (str, optional): The level of comparison, current version support, protein and peptide. Defaults to 'protein'.

        Raises:
            ValueError: If the level is not a peptide or protein

        Returns:
            float: The reproducibility rate among the replicates. 
        """
        if level not in ['peptide','protein']:
            raise ValueError(f"The provided reproducibility level: {level} is not supported, supported levels are: peptide and protein levels.")
        if level=='protein':
            return len(self.get_concise_proteins())/len(self.get_union_protein())
        elif level=='peptide':
            return len(self.get_concise_peptides())/len(self.get_union_peptides())
        
    def get_union_protein(self)->List[str]:
        """return a list of all proteins defined among the replicates, i.e. the union accross all replicates 

        Returns:
            List[str]: a list of protein ids identified among alll experiments 
        """
        return self._exps.get_unique_proteins()

    def get_concise_proteins(self)->List[str]:
        """ return the list of proteins IDss for that have been identified in all replicates 
        Returns:
            List[str]: [description]
        """
        return elf._exps.get_proteins_present_in_all()

    def get_union_peptides(self)->List[str]:
        """ return a list of all peptides detected among replicates, i.e. the union accross all replicates 
        Returns:
            List[str]: a list of all peptides identifed among all replicates
        """
        return self._exps.get_unique_peptides()
         
    def get_concise_peptides(self)->List[str]:
        """ return the list of peptide IDss for that have been identified in all replicates 
        Returns:
            List[str]: [description]
        """
        return self._exps.get_peptides_present_in_all()

    def get_experiments(self)->ExperimentSet:
        """
        Returns:
            ExperimentSet: the instance ExperimentSet object containing all the elements 
        """
        return self._exps 

    def cash_results(self, result_name:str,result_value:object):
        """append results to the instance hash map of results.
        Args:
            result_name (str): the name
            result_value (object): [description]
        """
        self._cashed_results[result_name]=result_value