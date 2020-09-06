#!/usr/bin/env python 
"""
@brief: unit testing the function simulate_an_experimental_ident_table_from_fasta
@author: Hesham ElAbd
"""
# load the modules 
import unittest
from IPTK.Utils.DevFunctions import simulate_an_experimental_ident_table_from_fasta
# Define Some Values 
PATH_TO_LOAD_DATA='../../assets/e_coli.fasta'
NUM_PEP=1000
NUM_PROTEIN=10
# define the test class 
class TestCaseSimulateExperiment(unittest.TestCase):
    def test_correct_input(self):
        self.assertEqual(
            simulate_an_experimental_ident_table_from_fasta(
               PATH_TO_LOAD_DATA, NUM_PEP,NUM_PROTEIN
            ).shape[0], NUM_PEP
        )
    
    def test_missing_num_peptide_and_proteins(self):
        with self.assertRaises(ValueError):
            simulate_an_experimental_ident_table_from_fasta(
                PATH_TO_LOAD_DATA, NUM_PEP-NUM_PEP,NUM_PROTEIN-NUM_PROTEIN
            )
    def test_correct_behavior(self):
        self.assertEqual(
            simulate_an_experimental_ident_table_from_fasta(
        PATH_TO_LOAD_DATA,  NUM_PEP, NUM_PROTEIN+1
        ).shape[0],
        990)



if __name__=='__main__':
    unittest.main()