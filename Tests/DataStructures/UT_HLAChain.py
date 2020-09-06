#!/usr/bin/env python 
"""
@brief: Unit testing for HLA Class ‚
@author: Hesham ElAbd
"""
# load the modules 
import unittest
from IPTK.DataStructure.HLAChain import HLAChain
# define some test variables 
TEST_CHAIN=['A*03:01','HLA-B*04:02', 'HLA-C*02', 'HLA-DRB1*13:01']
# define the test class 
class TestHLACase(unittest.TestCase):
    def setUp(self):
        self._hla=HLAChain('DRB1*15:01')

    def test_init__(self):
        for chain in TEST_CHAIN:
            HLAChain(chain)
            
    def test_get_class(self):
        self.assertEqual(self._hla.get_class(),2)
    
    def test_get_gene(self):
        self.assertEqual(self._hla.get_gene(),'DRB1')

    def test_get_hla_group(self):
        self.assertEqual(self._hla.get_allele_group(),'15')
    
    def test_get_hla_protein(self):
        self.assertEqual(self._hla.get_protein_group(),'01')
    
    def test_str__(self):
        string_res=f"""An HLA chain of class: 2 from gene: DRB1,
                With an allele group of: 15 and a protein group of: 
                01   
                """
        self.assertEqual(str(self._hla),string_res)
    
    def test_repr__(self): 
        string_res=f"""An HLA chain of class: 2 from gene: DRB1,
                With an allele group of: 15 and a protein group of: 
                01   
                """
        self.assertEqual(repr(self._hla),string_res)
if __name__=='__main__' :
    unittest.main() 
    
