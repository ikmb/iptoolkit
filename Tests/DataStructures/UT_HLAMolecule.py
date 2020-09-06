#!/usr/bin/env python 
"""
@brief: Unit testing for HLAMolecules 
@author: Hesham ElAbd
"""
# load the modules 
import unittest
from IPTK.DataStructure.HLAMolecules import HLAMolecule
# define some test-alleles 
HLA_ALLELE_OK_OC=['HLA-DRB1*15:01','HLA-A*01:01','B*04:01']
HLA_ALLELE_OK_TC=['DQA1*0303:DQB1*0202','DQA1*0201:DQB1*0303']
HLA_ALLELE_NOT_OK=['L*15:01','HLA-A*01011','']
# define the test-case class 
class TestHLAMolecules(unittest.TestCase):
    def setUp(self):
        self._hla_molecules=HLAMolecule(chain_alpha='DQA1*02:01',chain_beta='DQB1*03:03')
    
    def test_init__(self):
        for allele in HLA_ALLELE_OK_OC:
            HLAMolecule(allele=allele)
        
        for allele in HLA_ALLELE_OK_TC:
            chains=allele.split(':')
            HLAMolecule(chain_alpha=chains[0],chain_beta=chains[1])

        with self.assertRaises(RuntimeError):
            for allele in HLA_ALLELE_NOT_OK:
                HLAMolecule(allele=allele)
        
        with self.assertRaises(ValueError):
            HLAMolecule()

        with self.assertRaises(RuntimeError):
            for allele in HLA_ALLELE_NOT_OK: 
                HLAMolecule(allele=allele)

    def test_get_names(self):
        self.assertEqual(self._hla_molecules.get_name(),'DQA1*0201:DQB1*0303')
    
    def test_get_class(self):
        self.assertEqual(self._hla_molecules.get_class(),2)

    def test_get_gene(self):
        # assert the correct number of elements is returned 
        genes=self._hla_molecules.get_gene()
        self.assertEqual(len(genes), 2)
        # assert the correctness of the content 
        self.assertIn('DQA1',genes)
        self.assertIn('DQB1',genes)
    
    def test_get_allele_group(self):
        # assert the correct number of elements is returned 
        genes=self._hla_molecules.get_allele_group()
        self.assertEqual(len(genes), 2)
        # assert the correctness of the content 
        self.assertIn('02',genes)
        self.assertIn('03',genes)
    
    def test_get_protein_group(self):
        # assert the correct number of elements is returned 
        genes=self._hla_molecules.get_protein_group()
        self.assertEqual(len(genes), 2)
        # assert the correctness of the content 
        self.assertIn('01',genes)
        self.assertIn('03',genes)

if __name__=='__main__':
    unittest.main()