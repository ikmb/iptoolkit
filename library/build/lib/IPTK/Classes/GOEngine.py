"""The class provides wrapper for performing GO analysis for the identified protein hits, the engine is wrapper around the 
GOATools library.   
"""
## import the modules
from __future__ import annotations 
import os 
import time 
from typing import List, Set, Dict, Union
import pandas as pd
from IPTK.Classes.DataFeeders.gene_ncbi_9606_prot_encoding import GENEID2NT as gene2iden_human
from IPTK.Classes.DataFeeders.gene_ncbi_10090_prot_encoding import GENEID2NT as gene2iden_mouse 
from IPTK.Utils.Mapping import map_from_uniprot_to_Entrez_Gene
from IPTK.Classes.Experiment import Experiment 
import tqdm as tqdm 
from goatools.base import download_go_basic_obo,download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.go_enrichment import GOEnrichmentRecord

## define the engine class 
class GOEngine:
    def __init__(self, work_dir: str='.', clean_work_dir: bool = False,
                organism:str='human',
                study_parameters:Dict[str,Union[int,float,str,List,Dict]]={
                    'propagate_counts':False,
                    'alpha':0.05, 
                    'methods':['fdr_bh']
                })->GOEngine:
        """A GOEngine that can be used for performing analysis using GOATOOLS

        Args:
            work_dir (str, optional): The path to a temp directory were intermediate-results and raw data will be downloaded/written to. Defaults to the current working directory.
            clean_work_dir (bool, optional): Whether or not to remove data written to the work directory at class termination, default to True.
            organism (str, optional): The organism . Defaults to 'human'.
            study_parameters (Dict[str,Union[int,float,str,List,Dict]], optional): A dict of parameters to control the base function, defaults to {'propagate_counts':False,'alpha':0.05, 'methods':['fdr_bh']}
        Returns:
            GOEngine: return a GO engine that can be used for performing GO enrichment analysis GOEnrichmentStudyNS
        """
        print("Creating a GO Engine ...")
        if not os.path.exists(work_dir):
            raise ValueError(f"The provided work path: {work_dir} does not exist!!!")
        self.work_dir=work_dir
        if organism != 'human' and organism!='mouse':
            raise ValueError(f"The provided organism: {organism} is not support, current engine mainly work with human and moues only")
        print(f"\t --> Downloading data ...")
        obo_fname=download_go_basic_obo(os.path.join(work_dir,'go-basic.obo'))
        gene2go_fname=download_ncbi_associations(os.path.join(work_dir,'gene2go'))
        ## parse the GO term 
        print(f"\t --> parsing the data and intializing the base GOEA object...")
        obo_dag=GODag(obo_fname)
        if organism=='human':
            self._goea_obj=GOEnrichmentStudyNS(
                gene2iden_human.keys(), 
                Gene2GoReader(gene2go_fname,taxids=[9606]).get_ns2assc(),
                obo_dag,
                **study_parameters)
        else:
            self._goea_obj=GOEnrichmentStudyNS(
                gene2iden_human.keys(), 
                Gene2GoReader(gene2go_fname,taxids=[10090]).get_ns2assc(),
                obo_dag,
                **study_parameters)
        self._clean_work_dir=clean_work_dir
        self._gene_ids=None
        return 
    
    def load_data(self,exp:Experiment, num_proteins:int=-1)->None:
        """Load the data to the Engine, so GOEA can be conducted 

        Args:
            exp (Experiment): An Experimental object to extract uniprot ids 
            num_proteins (int, optional): The number of proteins to be included in the analysis. Defaults -1 to which mean use all proteins,\
                 otherwise it uses the number of proteins provided by the user. note that the function is sorted by number of peptides per protein,\
                      that is the first 10 protein means, getting the top 10 protein with most peptides. 
        Raises:
            ValueError: if the function called while data being already associated with the engine from a previous call
        """
        if self._gene_ids is not None:
            raise ValueError(f"There some data still in the engine, the first 10 genes are: {','.join(self._gene_ids[:10])}\
                clean your engine from previous data using the function, clean_engine and try again.")
        print(f"Getting the number of peptide per protein ..., started at: {time.ctime()}")
        num_protein_per_peptides=exp.get_peptides_per_protein()
        if num_proteins == -1: 
            list_proteins=num_protein_per_peptides.iloc[:,0].to_list()
        else:
            list_proteins=num_protein_per_peptides.iloc[:,0].to_list()[:num_proteins]
        print(f"Map uniprot to Entrez gene ids ..., starting at: {time.ctime()}")
        self._gene_ids= [int(gene_id) for gene_id in map_from_uniprot_to_Entrez_Gene(list_proteins).iloc[:,1].to_list()]
        print(f"{len(self._gene_ids)} Genes have been correctly loaded")
        return 

    def run_analysis(self, quite:bool=False, only_signifcant:bool=True,
            significance_level:float=0.05, get_list_term:bool=False)->Union[pd.DataFrame, List[GOEnrichmentRecord]]:
        if quite:
            goea_results=self._goea_obj.run_study(self._gene_ids, prt=None)
        else:
            goea_results=self._goea_obj.run_study(self._gene_ids)
        if only_signifcant:
            goea_results=[res for res in goea_results if res.p_fdr_bh<significance_level]
        if get_list_term:
            return goea_results
        else: 
            self._goea_obj.wr_tsv(os.path.join(self.work_dir,'temp_file.tsv'),goea_results)
            results_df=pd.read_csv(os.path.join(self.work_dir,'temp_file.tsv'),sep='\t')
            os.system(f"rm -f {os.path.join(self.work_dir,'temp_file.tsv')}")
        return results_df

    def clean_engine(self)->None:
        """Remove Current list of gene ids associated with the engine 
        """
        self._gene_ids=None
        return 

    def __del__(self)->None:
        """class destructor, clean work directory if  clean_work_dir is set to True 
        """
        if self.clean_work_dir: os.system(f"rm -f {self.work_dir}/*")
        return 


