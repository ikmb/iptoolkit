#!/usr/bin/env python 
"""
@brief: The Unit testing of the HLASet class
@author: Hesham ELAbd 
"""
# load the modules 
import unittest
from IPTK.DataStructure.HLASet import HLASet
# define some Test Variable 
OK_TEST_SETS=[
    ['HLA-DRB1*15:01','DRB1*13:01'],
    ['HLA-DRB1*01:01','DQA1*05:01$DQB1*02:01'],
    ['HLA-A*01:01','HLA-B*04:01','C*03:02']
]
NOT_OK_TEST_SETS=[
    ['HLA-A*01:01','DRB1*13:01'], 
    ['HLA-DRB1*01:01','DQA1*05:01:DQB1*02:01', 'DQA1*0501:DQB1*02:01'],
]
# define the class 
class TestHLASetCase(unittest.TestCase):
    def setUp(self):
        self._hla_2=HLASet(OK_TEST_SETS[0])
        self._hla_2dq=HLASet(OK_TEST_SETS[1],gene_sep='$')
        self._hla_1=HLASet(OK_TEST_SETS[-1])
    
    def test_init__(self):
        # test initialization with a mixture of HLA-I / HLA-II  
        with self.assertRaises(ValueError):
            HLASet(NOT_OK_TEST_SETS[0])
        # test the initialization with names that can be parsed correctly 
        with self.assertRaises(RuntimeError): 
            HLASet(NOT_OK_TEST_SETS[1])
    
    def test_get_hla_count(self):
        self.assertEqual(self._hla_2.get_hla_count(),2)
        self.assertEqual(self._hla_2dq.get_hla_count(),2)
        self.assertEqual(self._hla_1.get_hla_count(),3)
    
    def test_len__(self):
        self.assertEqual(len(self._hla_2),2)
        self.assertEqual(len(self._hla_2dq),2)
        self.assertEqual(len(self._hla_1),3)
        
    def test_get_class(self):
        self.assertEqual(self._hla_2.get_class(),2)
        self.assertEqual(self._hla_2dq.get_class(),2)
        self.assertEqual(self._hla_1.get_class(),1)
    
    def test_has_alleles(self):
        self.assertTrue(self._hla_2.has_allele('DRB1*15:01'))
        self.assertTrue(self._hla_2.has_allele('HLA-DRB1*15:01'))
        self.assertTrue(self._hla_2dq.has_allele('HLA-DRB1*01:01'))
        self.assertTrue(self._hla_2dq.has_allele('HLA-DQA1*05:01$DQB1*02:01'))
        self.assertTrue(self._hla_1.has_allele('A*01:01'))
        self.assertTrue(self._hla_1.has_allele('B*04:01'))
        self.assertTrue(self._hla_1.has_allele('HLA-C*03:02'))

    def test_has_gene(self):
        self.assertTrue(self._hla_2.has_gene('DRB1'))
        self.assertTrue(self._hla_2dq.has_gene('DQA1'))
        self.assertTrue(self._hla_2dq.has_gene('DQB1'))
        self.assertTrue(self._hla_1.has_gene('B'))
        self.assertTrue(self._hla_1.has_gene('A'))
        self.assertTrue(self._hla_1.has_gene('C'))

    def test_has_allele_group(self):
        self.assertTrue(self._hla_2.has_allele_group('15'))
        self.assertTrue(self._hla_2.has_allele_group('13'))
        self.assertTrue(self._hla_2dq.has_allele_group('01'))
        self.assertTrue(self._hla_2dq.has_allele_group('05'))
        self.assertTrue(self._hla_2dq.has_allele_group('02'))
        self.assertTrue(self._hla_1.has_allele_group('01'))
        self.assertTrue(self._hla_1.has_allele_group('04'))
        self.assertTrue(self._hla_1.has_allele_group('03'))
    
    def test_has_protein_group(self):
        self.assertTrue(self._hla_2.has_protein_group('01'))
        self.assertTrue(self._hla_2.has_protein_group('01'))
        self.assertTrue(self._hla_2dq.has_protein_group('01'))
        self.assertTrue(self._hla_2dq.has_protein_group('01'))
        self.assertTrue(self._hla_2dq.has_protein_group('01'))
        self.assertTrue(self._hla_1.has_protein_group('01'))
        self.assertTrue(self._hla_1.has_protein_group('01'))
        self.assertTrue(self._hla_1.has_protein_group('02'))

    def test_str__(self):
        self.assertEqual(str(self._hla_2),'An HLASet containing 2 alleles')
        self.assertEqual(str(self._hla_2dq),'An HLASet containing 2 alleles')
        self.assertEqual(str(self._hla_1),'An HLASet containing 3 alleles')
    
    def test_repr__(self):
        self.assertEqual(repr(self._hla_2),'An HLASet containing 2 alleles')
        self.assertEqual(repr(self._hla_2dq),'An HLASet containing 2 alleles')
        self.assertEqual(repr(self._hla_1),'An HLASet containing 3 alleles')

if __name__=='__main__':
    unittest.main()