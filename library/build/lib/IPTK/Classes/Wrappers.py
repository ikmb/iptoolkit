"""A collection of convient classes to construct different IPTK classes 
"""
## load the modules and import other classes 
from __future__ import annotations
import os
import re
import time
from typing import List, Dict, Tuple, Union 
from IPTK.Classes.Experiment import Experiment
from IPTK.Classes.ExperimentSet import ExperimentSet
from IPTK.Classes.Database import SeqDB, GeneExpressionDB, CellularLocationDB,OrganismDB
from IPTK.Classes.Tissue import Tissue
from IPTK.Classes.HLASet import HLASet
from IPTK.Classes.Proband import Proband
from IPTK.Classes.Peptide import Peptide 
from IPTK.IO.InFunctions import (parse_mzTab_to_identification_table, 
parse_xml_based_format_to_identification_table, parse_text_table,parse_mzIdentML_to_identification_table) 
from tqdm import tqdm 
from concurrent.futures import ProcessPoolExecutor
from concurrent import futures
import multiprocessing as mp 
## Define the experiment class 
class RExperiment:
    def __init__(self,filepath:str, path2fasta:str, fileformat:str='idXML',tissue_name:str='total PBMC',
                proband_name:str='Default Proband', hla_set:List[str]=['DRB1*15:01','DRB1*15:01'],
                parser_param:Dict[str,Union[list,set,dict,int,float]]={})->RExperiment:
        """A Wrapper class for constracting an experimental dataset using user defined parameters\
            The class take care of initializing all classes and functions provided an easy-to-use interface\
            for working with immunopeptidomics data 
        Args:
            filepath (str): the path to load the input file, for example and idXML or an Identification table 
            path2fasta (str): the path to load Fasta database 
            fileformat (str, optional): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table.\
                 Defaults to 'idXML'.
            tissue_name (str, optional): The name of the tissue to utilize, this is used for initializing the gene expression table\
                 Defaults to 'total PMBC'.
            proband_name (str, optional): the name of the proband from whome the data was obtained. Defaults to 'Default Proband'.
            hla_set (List[str], optional): A list of HLA alleles from whome the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].
            parser_param (Union[list,set,dict,int,float], optional): A list of parameters to be forwarded the file parser
        Returns:
            cExperimet: an IPTK.Class.Wrapper.Experiment class, an IPTK.Class.Experiment.Experiment can be extracted from the resutned instance using the get_experiment method 
        """
        ## Checking that the input is correct 
        if not os.path.exists(filepath):
            raise ValueError(f"The provided path for the identification file : {filepath} does not exist!!")
        if not os.path.exists(path2fasta):
            raise ValueError(f"The path to the proivded fasta file: {path2fasta}, does not exists!!!")
        if fileformat not in ['idXML', 'pepXML', 'csv','mzTab','mzid']:
            raise ValueError(f"Unknow input format, the provided format: {fileformat} is not supported, currently supported values are: {', '.join(['idXML', 'pepXML', 'csv','mzTab','mzid'])}")
        # define the data 
        self._proband = Proband(name=proband_name)  # the name of the proband 
        try: 
            self._hLASet = HLASet(hlas=hla_set) # just a place holder to represent the HLA allele, an instance of class HLASet
        except Exception as exp:
            raise RuntimeError(f"The following error was Encountered while creating an HLASet: \n{str(exp)}\n")
        try: 
            self._seqBD= SeqDB(path2fasta=path2fasta)
        except Exception as exp:
            raise IOError(f"While loading the fasta database the following error was Encountered : \n{str(exp)}\n")
        self._expresson_profile = GeneExpressionDB() # use the data on the human protein atlas @https://www.proteinatlas.org/about/download --> Normal tissue data 
        self._protein_locations = CellularLocationDB() # use the data on the human protein atlas @https://www.proteinatlas.org/about/download --> Subcellular location data
        try:
            self._tissue= Tissue(name=tissue_name,main_exp_value=self._expresson_profile, 
                        main_location=self._protein_locations) # create the tissue instance
        except Exception as exp:
            raise RuntimeError(f"While creating a tissue instance, the following error was Encountered: \n{str(exp)}\n")
        try:
            if fileformat == 'idXML':
                    input_table = parse_xml_based_format_to_identification_table(path2XML_file= filepath, path2fastaDB=path2fasta, is_idXML= True, **parser_param)
            elif fileformat == 'pepXML':
                     input_table = parse_xml_based_format_to_identification_table(path2XML_file= filepath, path2fastaDB=path2fasta, is_idXML= False, **parser_param)
            elif fileformat == 'mzTab':
                input_table = parse_mzTab_to_identification_table(path2mzTab= filepath, path2fastaDB=path2fasta, **parser_param)
            elif fileformat=='':
                input_table = parse_mzIdentML_to_identification_table(filepath)
            else:
                input_table = parse_text_table(filepath,path2fasta, **parser_param)
        except Exception as exp:
            raise ValueError(f"Loading the input table has caused to the following error: \n{str(exp)}\n")
        # constructing the experiments 
        try: 
            self._exp=Experiment(self._proband, self._hLASet, self._tissue, self._seqBD, input_table)
        except Exception as exp:
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

class ReplicatedExperiment:
    def __init__(self,path:str, anchor_name:str ,path2fasta:str,
        fileformat:str='idXML', tissue_name:str='total PBMC', proband_name:str='Default Proband',
        hla_set:List[str]=['DRB1*15:01','DRB1*15:01'],
        parser_param:Dict[str,Union[list,set,dict,int,float]]={})->ReplicatedExperiment:
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
            proband_name (str, optional): the name of the proband from whom the data was obtained. Defaults to 'Default Proband'.
            hla_set (List[str], optional): A list of HLA alleles from whom the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].
            parser_param (Union[list,set,dict,int,float], optional): A list of parameters to be forwarded the file parser.

        Returns:
            PairedExperiment: [description]
        """
        ## Check that the provided path exist 
        if not os.path.exists(path):
            raise ValueError(f"The provided root directory: {path}, does not exist")
        ## get path to samples replicates 
        replicates=[os.path.join(path, name) for name in os.listdir(path) if re.match(anchor_name, name) ]
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
                exps[parsed_name]=RExperiment(file_name, path2fasta, fileformat, tissue_name, proband_name, hla_set,parser_param).get_experiment()
        except Exception as exp:
            raise RuntimeError(f"Creating an Experimental instance from the file: {file_name} failed with the following error \n{str(exp)}\n")
        self._exps=ExperimentSet(**exps)
        self._cashed_results=dict()
        return  
    
    def get_repro_rate(self, level='protein')->float:
        """Compute the reproducibility rate among the replicates through Jaccard Index
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
            List[str]: a list of protein ids identified among all experiments 
        """
        return self._exps.get_unique_proteins()

    def get_concise_proteins(self)->List[str]:
        """ return the list of proteins IDss for that have been identified in all replicates 
        Returns:
            List[str]: [description]
        """
        return self._exps.get_proteins_present_in_all()

    def get_union_peptides(self)->List[str]:
        """ return a list of all peptides detected among replicates, i.e. the union accross all replicates 
        Returns:
            List[str]: a list of all peptides identified among all replicates
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
            result_value (object): the results to be associated with the provided name 
        """
        self._cashed_results[result_name]=result_value
    
class RExperimentSet:
    def __init__(self,path:str,path2fasta:List[str],
        fileformat:List[str]=['idXML'], tissue_name:str='total PMBC',
        proband_name:List[str]=['Default Proband'],
        hla_set:List[List[str]]=[['DRB1*15:01','DRB1*15:01']],
        num_worker:int=mp.cpu_count(),
        parser_param:Dict[str,Union[list,set,dict,int,float]]={}
        )->RExperimentSet:
        """Construction an ExperimentSet from all the input identification table

        Args:
            path (Union[str,os.path]): The path to the root directory. 
            path2fasta (str): the path to load Fasta database 
            fileformat (str, optional): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table (IdTable).\
                Defaults to 'idXML'.
            tissue_name (str, optional): The name of the tissue to utilize, this is used for initializing the gene expression table\
                Defaults to 'total PMBC'.
            proband_name (str, optional): the name of the proband from whom the data was obtained. Defaults to 'Default Proband'.
            hla_set (List[str], optional): A list of HLA alleles from whom the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].
            parser_param (Union[list,set,dict,int,float], optional): A list of parameters to be forwarded the file parser.
        Raises:
            ValueError: in case no file could not be found in the root directory or if the length of parameter mismatched

        Returns:
            RExperimentSet: An RExperimentSet containing a collection of experiments each of which correspond to a file found in the provided root directory 
        """
        ## Check that the path exists 
        if not os.path.exists(path):
            raise ValueError(f"The provided root directory: {path} does not exists!")
        # Get the filenames relative to the path base 
        filenames=[os.path.join(path,name) for name in os.listdir(path) if fileformat in name] 
        # Check that the path exists  
        if len(filenames)==0:
            raise ValueError(f"Could not find any file in the provided directory: {path} with the provided format: {fileformat}") 
        # build the collection of experiments 
        list_jobs=[]
        print(f"Creating an ExperimentSet from: {len(filenames)} files using: {num_worker}, starting at: {time.ctime()}")
        # submitting the jobs to a process pool 
        with ProcessPoolExecutor(num_worker) as exc:
            for idx in range(len(filenames)):
                parsed_name=filenames[idx].split('/')[-1].split('.')[0]
                list_jobs.append(exc.submit(build_experiments,parsed_name,filenames[idx],path2fasta,
                fileformat, tissue_name, proband_name, hla_set,parser_param
                ))
        success_process=0
        exps=dict()
        for fut in futures.as_completed(list_jobs):
            name,experiment=fut.result()
            exps[name]=experiment
            success_process+=1
            print(f"Progress: {success_process} files have been read and parsed into experiment, progress is: {(success_process/len(filenames))*100}%")
        ## create instance resources 
        self._exps=ExperimentSet(**exps)
        self._cashed_results=dict()
        return

    def get_experiments(self)->ExperimentSet:
        """ 
        Returns:
            ExperimentSet: The instance ExperimentalSet 
        """
        return self._exps

    def cash_results(self, result_name:str,result_value:object):
        """append results to the instance hash map of results.
        Args:
            result_name (str): the name
            result_value (object): the results to be associated with the provided name 
        """
        self._cashed_results[result_name]=result_value
        
class ReplicatedExperimentSet:
    def __init__(self,path:str,path2fasta:List[str],
        fileformat:List[str]=['idXML'], tissue_name:List[str]=['total PMBC'],
        proband_name:List[str]=['Default Proband'],
        hla_set:List[List[str]]=[['DRB1*15:01','DRB1*15:01']],
        num_worker:int=mp.cpu_count(),
        parser_param:Dict[str,Union[list,set,dict,int,float]]={})->ReplicatedExperimentSet:
        """
        Args:
            path (Union[str,os.path]): The path to the root directory where files are located
            path2fasta (str): the path to load Fasta database 
            fileformat (str, optional): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table (IdTable).\
                 Defaults to 'idXML'.
            tissue_name (str, optional): The name of the tissue to utilize, this is used for initializing the gene expression table\
                 Defaults to 'total PMBC'.
            proband_name (str, optional): the name of the proband from whome the data was obtained. Defaults to 'Default Proband'.
            hla_set (List[str], optional): A list of HLA alleles from whome the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].            
            parser_param (Union[list,set,dict,int,float], optional): A list of parameters to be forwarded the file parser.

        Raises:
            ValueError: in case no file could not be found in the root directory or if the length of parameter mismatched

        Returns:
            ReplicatedExperimentSet: a collection of replicated Experiments. 
        """
        ## Check that file name exists 
        if not os.path.exists(path):
            raise ValueError(f"The provided path: {path}, does not exists")
        ## list file names:
        filenames=[os.path.join(path,name) for name in os.listdir(path) if fileformat in name] 
        ## load the file names 
        if len(filenames)==0:
            raise ValueError(f"Could not find any file in the provided directory: {path} with the provided format: {fileformat}")
        ## get experimental names:
        replicates=self._group_paired_files(filenames)
        ## allocate a container to hold the results 
        list_jobs=[]
        ## Create an executioner to submit results too it 
        print(f"Creating an ExperimentSet from: {len(filenames)} files using: {num_worker}, starting at: {time.ctime()}")
        with ProcessPoolExecutor(num_worker) as exc:
            for replicate in replicates:
                anchor_name=self._get_common_prefix(replicate[0],replicate[1])
                list_jobs.append( exc.submit(
                    build_repeated_experiments,path,anchor_name,path2fasta, fileformat, tissue_name, proband_name, hla_set,parser_param
                ))
        self._results=dict()  
        success_process=0
        for fut in futures.as_completed(replicates):
            anchor_name,rep_exp=fut.result()
            self._results[anchor_name]=rep_exp
            success_process+=1
            print(f"Progress: {success_process} files have been read and parsed into experiment, progress is: {(success_process/len(replicates))*100}%")
        return 

    def _group_paired_files(self, file_names:List[str])->List[List[str]]:
        """ takes a list of input files and returns a list of lists where each child list contains replicate files

        Args:
            file_names (List[str]): a collection of file names 

        Returns:
            List[[str]]: a list of lists where each child list contains replicate files
        """
        results=[]
        for file_name_outer in file_names:
            for file_name_inter in file_names:
                if self._are_replicates(file_name_outer,file_name_inter):
                    results.append([file_name_outer,file_name_inter])
        # we need to make sure we have unique elements 
        for replicate in results:replicate.sort()
        temp_dev=set(['$'.join(replicate) for replicate in results])
        results=[replicate.split('$') for replicate in temp_dev]
        return results 

    def _are_replicates(self, filename_1:str, filename_2:str)->bool:
        """a helper function that check whether two files are replicate or not be comparing their names

        Args:
            filename_1 (str): The filename of the first experiment 
            filename_2 (str): The filename of the second experiment 

        Returns:
            bool: True if the two filenames are identical strings that only differ in one charachter and False otherwise. 
        """
        if len(filename_1)!=len(filename_2):return False 
        if filename_1==filename_2: return False 
        diff=[]
        index_mismatch=[]
        counter=0
        for a,b in zip(filename_1,filename_2):
            if a!=b:
                diff.append([a,b])
                index_mismatch.append(counter)
            if len(diff)>1:
                return False
            counter+=1
        if diff[0][0] not in '123456789' and diff[0][1] not in '123456789':
            return False
        if filename_1[index_mismatch[0]+1]!='_' or filename_2[index_mismatch[0]+1]!='_':
            return False
        return True 

    def _get_common_prefix(self,file_name1:str,file_name2:str)->str:
        """take an equal length file names and return the shared prefix between both files 

        Args:
            file_name1 (str): The name of the first file 
            file_name2 (str): The name of the second file 

        Returns:
            str: the prefix string between both files 
        """
        if len(file_name1)!=len(file_name2):
            raise ValueError(f"Provided String MUST be of equal length, found first file to have a length of {en(file_name1)} while the second has a length of {en(file_name2)}")
        idx=0
        file_name1=file_name1.split('/')[-1]
        file_name2=file_name2.split('/')[-1]
        for a, b in zip(file_name1,file_name2):
            if a!=b:
                break
            idx+=1
        return file_name1[:idx]+'(1|2|3|4|5|6|7|8|9)'+file_name1[idx+1:]

    def get_experiments(self)->Dict[str,ReplicatedExperiment]:
        """
        Returns:
            Dict[str,ReplicatedExperiment]: a dict containg a collection of ReplicatedExperiments  
        """
        return self._results
## Define helper analysis function use for multiprocessing 
#---------------------------------------------------------
def build_experiments(parsed_name:str, file_name:str, path2fasta:str, fileformat:str,
                    tissue_name:str, proband_name:str, hla_set:str,parser_param:Dict[str,Union[list,set,dict,int,float]]={})->Tuple[str,Experiment]:
    """[summary]

    Args:
        parsed_name (str): The name of the experiment in the results dictionary 
        file_name (str): The path to load the identification table. 
        path2fasta (str): The path to load the fasta identification files
        fileformat (str): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table (IdTable).\
                 Defaults to 'idXML'.
        tissue_name (str): The name of the tissue to utilize, this is used for initializing the gene expression table\
                 Defaults to 'total PMBC'.
        proband_name (str): the name of the proband from whome the data was obtained. Defaults to 'Default Proband'.
        hla_set (str): A list of HLA alleles from whome the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].
        parser_param (Union[list,set,dict,int,float], optional): A list of parameters to be forwarded the file parser.   

    Returns:
        Experiment: An experimental object constructed from the provided parameter 
    """
    return parsed_name, RExperiment(file_name, path2fasta, fileformat, tissue_name, proband_name, hla_set,parser_param).get_experiment()
###========================================================================================================================
def build_repeated_experiments(path:str, anchor_name:str, path2fasta:str, fileformat:str,
                    tissue_name:str, proband_name:str, hla_set:str,
                    parser_param:Dict[str,Union[list,set,dict,int,float]]={}
                    )->Tuple[str,Experiment]:
    """Build a new ReplicatedExperiment instance 

    Args:
        parsed_name (str): The name of the experiment in the results dictionary 
        file_name (str): The path to load the identification table. 
        path2fasta (str): The path to load the fasta identification files
        fileformat (str): type of input format, can be any of idXML, pepXML, mzTab or a CSV Table (IdTable).\
                 Defaults to 'idXML'.
        tissue_name (str): The name of the tissue to utilize, this is used for initializing the gene expression table\
                 Defaults to 'total PMBC'.
        proband_name (str): the name of the proband from whome the data was obtained. Defaults to 'Default Proband'.
        hla_set (str): A list of HLA alleles from whome the data was obtained. Defaults to ['DRB1*15:01','DRB1*15:01'].   
        parser_param (Union[list,set,dict,int,float], optional): A list of parameters to be forwarded the file parser.   

    Returns:
        ReplicatedExperimentSet: a collection of replicated Experiments. 
    """
    return anchor_name, ReplicatedExperiment(path, anchor_name, path2fasta, fileformat, tissue_name, proband_name, hla_set,parser_param)

