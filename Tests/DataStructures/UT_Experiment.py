#!/usr/bin/env python 
"""
@brief: Uni-testing the experiment class
@author: Hesham ElAbd
"""
# load the modules 
import unittest
import pandas as pd
import os
from typing import List, Dict, Tuple, Set
import numpy as np
from IPTK.DataStructure.Experiment import Experiment
from IPTK.Utils.DevFunctions import simulate_an_experimental_ident_table_from_fasta, simulate_an_expression_table
from IPTK.DataStructure.Proband import Proband
from IPTK.DataStructure.HLASet import HLASet 
from IPTK.DataStructure.Tissue import Tissue
from IPTK.DataStructure.Database import SeqDB
from IPTK.DataStructure.Peptide import Peptide
# define some test variables 
DATA_PATH: str ='../assets/e_coli.fasta'
NUM_PEP_SMALL: int =10
NUM_PROT_SMALL: int =5
TISSUE_NAME: str = 'Test_Tissue'
PROBAND_NAME: str = 'Test_Proband'
proband: Proband =Proband(name=PROBAND_NAME)
hla_set: HLASet =HLASet(['DRB1*15:01','DRB1*13:01'])
tissue: Tissue = Tissue('Test_Tissue', simulate_an_expression_table(num_transcripts=10))
ident_table: pd.DataFrame = simulate_an_experimental_ident_table_from_fasta(DATA_PATH,NUM_PEP_SMALL,NUM_PROT_SMALL)
ident_table_mini: pd.DataFrame = simulate_an_experimental_ident_table_from_fasta(DATA_PATH,1,1)
PROTIENS=list(set(ident_table.iloc[:,1]))
database: SeqDB = SeqDB(DATA_PATH)
# define the test case 
class TestCaseExperiment(unittest.TestCase):
    def setUp(self):
        self._exp=Experiment(proband,hla_set,tissue,database,ident_table)

    def test_add_org_info(self):
        # prepear the org_info
        org_info : Dict[str,str]= dict()
        # add the org-information
        for prot in PROTIENS:
            org_info[prot]='E.Coli'
        # add the org-info 
        self._exp.add_org_info(org_info)

    #  From Here 
    def test_get_flanked_peptides(self):
        flanked_peptide: List[str] = self._exp.get_flanked_peptides(5) 
        # compute some data to assert the correctness 
    
    def test_get_negative_example(self):
        negative_peptides: List[str]= self._exp.get_negative_example(5)

    def test_get_binarized_results(self):
        mapped_protein: List[np.ndarray]=self._exp.get_binarized_results()

    def test_get_peptide(self):
        for pep in ident_table.iloc[:,0].tolist():
            self._exp.get_peptide(pep)
    
    def test_get_mapped_protein(self):
        for prot in ident_table.iloc[:,1].tolist():
            self._exp.get_mapped_protein(prot)
        
    def test_get_mapped_proteins(self): 
        mapped_array: np.ndarray = self._exp.get_mapped_proteins()

    def test_get_mono_parent_peptides(self):
        self._exp.get_mono_parent_peptides()
    
    def test_get_poly_parental_peptides(self):
        self._exp.get_poly_parental_peptides()
    # Till here, we need to check for correctness 
    
    def test_get_number_of_children(self):
        for prot in PROTIENS:
            self.assertEqual(self._exp.get_number_of_children(prot),2)

    def test_get_peptide_per_protein(self):
        test_table: pd.DataFrame =self._exp.get_peptides_per_protein()
    
    def test_get_n_terminal_flanked_seqs(self):
        seqs: List[str]=self._exp.get_n_terminal_flanked_seqs(5)

    def test_get_c_terminal_flanked_seqs(self): 
        seqs: List[str]=self._exp.get_c_terminal_flanked_seqs(5)
    
    def test_get_tissue_name(self):
        self.assertEqual(
            self._exp.get_tissue_name(),TISSUE_NAME
        )

    def test_get_proband_name(self):
        self.assertEqual(self._exp.get_proband_name(),PROBAND_NAME)

    def test_get_hla_class(self):
        self.assertEqual(self._exp.get_hla_class(),2)
    
    def test_has_hla_allele(self):
        self.assertTrue(self._exp.has_hla_allele('DRB1*15:01'))
        self.assertTrue(self._exp.has_hla_allele('DRB1*13:01'))

    def test_has_gene(self):
        self.assertTrue(self._exp.has_gene('DRB1'))
    
    def test_has_allele_group(self):
        self.assertTrue(self._exp.has_allele_group('15'))
        self.assertTrue(self._exp.has_allele_group('13'))

    def test_has_protein_group(self):
        self.assertTrue(self._exp.has_protein_group('01'))
    
    def test_get_peptides(self):
        peptides: Set[str] = self._exp.get_peptides()
        for pep in peptides:
            self.assertIn(pep, ident_table.iloc[:,0].to_list())

    def test_get_proteins(self):
        proteins: Set[str]=self._exp.get_proteins()
        for prot in proteins:
            self.assertIn(prot, PROTIENS)

    def test_is_member(self):
        peptides: Set[str] = self._exp.get_peptides()
        for pep in peptides:
            self.assertTrue(self._exp.is_member(pep))
    
    def test_is_parent(self):
        for prot in PROTIENS:
            self.assertTrue(self._exp.is_a_parent_protein(prot))

    def test_len__(self):
        self.assertEqual(len(self._exp),NUM_PEP_SMALL)
    
    def test_str__(self):
        self.assertEqual(str(self._exp),f"""an Experimental from proband: {PROBAND_NAME}, Tissue: {TISSUE_NAME}
				   With an HLA Class: {2} With
				   {NUM_PEP_SMALL} peptide identification from {NUM_PROT_SMALL} Protein""")

if __name__=='__main__':
    unittest.main()