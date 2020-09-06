#!/usr/bin/env python 
"""
@brief: Unit testing for the tissue class 
@author: Hesham ElAbd
"""
# import the modules 
import unittest
import pandas as pd
from IPTK.Utils.DevFunctions import simulate_an_expression_table
from IPTK.DataStructure.Tissue import Tissue
# create the test cases 
NUM_TRANS=10
NUM_AUX_TRA=5
exp_table: pd.DataFrame = simulate_an_expression_table(NUM_TRANS)
aux_table: pd.DataFrame = simulate_an_expression_table(NUM_AUX_TRA)
TISSUE_NAME: str = 'TEST_TISSUE'
class TestCasesTissue(unittest.TestCase):
    def setUp(self):
        self._tissue=Tissue(TISSUE_NAME, exp_table,aux_table)
    
    def test_get_transcript_exp(self):
        self.assertEqual(self._tissue.get_transcript_expression(exp_table.iloc[:,0].tolist()[0]), 
        exp_table.iloc[:,2].tolist()[0]
        )

    def test_get_protein_name(self):
        self.assertEqual(self._tissue.get_protein_expression(exp_table.iloc[:,1].tolist()[0]), 
        exp_table.iloc[:,2].tolist()[0]
        )

    def test_get_names(self):
        self.assertEqual(self._tissue.get_name(),TISSUE_NAME)

    def test_init__(self):
        with self.assertRaises(ValueError):
            Tissue(None, exp_table,aux_table)
            Tissue(TISSUE_NAME,None,None)
    
    def test_str__(self):
        self.assertEqual(str(self._tissue),f'{TISSUE_NAME} with {NUM_TRANS+NUM_AUX_TRA} transcript expression value')

if __name__ == '__main__':
    unittest.main()
        
