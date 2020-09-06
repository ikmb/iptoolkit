#!/usr/bin/env python 
"""
@brief: A collection of test cases to test the Protein modules  
"""
# load the modules 
import unittest
from IPTK.DataStructure.Protein import Protein
import numpy as np 
# define some constants 
PROT_NAME="P11446"
PROT_SEQ="MLNTLIVGASGYAGAELVTYVNRHPHMNITALTVSAQSNDAGKLISDLHPQLKGIVDLPLQPMSDISEFSPGVDVVFLATAHEVSHDLAPQFLEAGCVVFDLSGAFRVNDATFYEKYYGFTHQYPELLEQAAYGLAEWCGNKLKEANLIAVPGCYPTAAQLALKPLIDADLLDLNQWPVINATSGVSGAGRKAAISNSFCEVSLQPYGVFTHRHQPEIATHLGADVIFTPHLGNFPRGILETITCRLKSGVTQAQVAQVLQQAYAHKPLVRLYDKGVPALKNVVGLPFCDIGFAVQGEHL"
PROT_LEN=len(PROT_SEQ)
ORG_NAME='E_COLI'
## TEST_DATA 
# _get_peptides_map data
BACKBONE=np.zeros(shape=(1,PROT_LEN))
START_IDX=[5,6,7,8]
END_IDX=[15,15,17,18]
for i,j in zip(START_IDX,END_IDX):
    BACKBONE[0,i:j]+=1
# _get_non_presented_peptide
EXCULTION_REGION_ONE=[1,25]
EXCULTION_REGION_TWO=[25,45]
EXCULTION_REGION_THREE=[270,300]
NUMBER_OF_TRIAL=100
SHORT_PEPTIDE=15
LONG_PEPTIDE=25
# define the test case 
class TestProtein(unittest.TestCase):
    def setUp(self):
        self._protein=Protein(PROT_NAME,PROT_SEQ)

    def test_get_id(self):
        self.assertEqual(self._protein.get_id(), PROT_NAME)

    def test_get_seq(self):
        self.assertEqual(self._protein.get_seq(),PROT_SEQ)
    
    def test_get_peptides_map(self):
        self.assertEqual(
            np.mean(self._protein.get_peptides_map(start_idxs=START_IDX,end_idxs=END_IDX)==BACKBONE),
            True)
        with self.assertRaises(ValueError):
            self._protein.get_peptides_map(START_IDX[1:],END_IDX)

    def test_get_non_presented_peptide(self): 
        # Test execution region one
        for _ in range(NUMBER_OF_TRIAL):
            # first with the short peptide 
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_ONE[0], EXCULTION_REGION_ONE[1], SHORT_PEPTIDE
                )
                ==PROT_SEQ[EXCULTION_REGION_ONE[0]:EXCULTION_REGION_ONE[1]]
                )
            # with longer peptides 
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_ONE[0], EXCULTION_REGION_ONE[1], LONG_PEPTIDE
                )
                ==PROT_SEQ[EXCULTION_REGION_ONE[0]:EXCULTION_REGION_ONE[1]]
                )
        # Test execution region one
        for _ in range(NUMBER_OF_TRIAL):
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_TWO[0], EXCULTION_REGION_TWO[1], SHORT_PEPTIDE
                )
                ==PROT_SEQ[EXCULTION_REGION_TWO[0]:EXCULTION_REGION_TWO[1]]
                )
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_TWO[0], EXCULTION_REGION_TWO[1], LONG_PEPTIDE
                )
                ==PROT_SEQ[EXCULTION_REGION_TWO[0]:EXCULTION_REGION_TWO[1]]
                )
        # Test execution region one
        for _ in range(NUMBER_OF_TRIAL):
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_THREE[0], EXCULTION_REGION_THREE[1], SHORT_PEPTIDE
                )
                ==PROT_SEQ[EXCULTION_REGION_THREE[0]:EXCULTION_REGION_THREE[1]]
                )
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_THREE[0], EXCULTION_REGION_THREE[1], LONG_PEPTIDE
                )
                ==PROT_SEQ[EXCULTION_REGION_THREE[0]:EXCULTION_REGION_THREE[1]]
                )
        # assert with longer peptide 
        for _ in range(NUMBER_OF_TRIAL):
            self.assertFalse(
                self._protein.get_non_presented_peptide(
                    EXCULTION_REGION_THREE[0], EXCULTION_REGION_THREE[1], 5
                )
                ==PROT_SEQ[EXCULTION_REGION_THREE[0]:EXCULTION_REGION_THREE[1]]
                )
        # assert the rise of the index error 
        with self.assertRaises(ValueError):
            self._protein.get_non_presented_peptide(EXCULTION_REGION_THREE[0],
            EXCULTION_REGION_THREE[1],PROT_LEN) 
        with self.assertRaises(ValueError):
            self._protein.get_non_presented_peptide(EXCULTION_REGION_THREE[0],
            EXCULTION_REGION_THREE[1],PROT_LEN+1)
        with self.assertRaises(ValueError):
            self._protein.get_non_presented_peptide(EXCULTION_REGION_THREE[0],
            EXCULTION_REGION_THREE[1],-1)  
            
    def test_set_org(self):
        self._protein.set_org(ORG_NAME)
        self.assertEqual(self._protein._org, ORG_NAME)
    
    def test_get_org(self): 
        self._protein.set_org(ORG_NAME)
        self.assertEqual(self._protein.get_org(),ORG_NAME)
    
    def test_getitem___(self): 
        self.assertEqual(self._protein[25:50],PROT_SEQ[25:50])
    
    def test_len__(self):
        self.assertEqual(len(self._protein),PROT_LEN)
    
    def test_str__(self): 
        self.assertEqual(str(self._protein),PROT_SEQ)

    def test_repr__(self):
        self.assertEqual(repr(self._protein),'A Protein instance with a length of 300')
    

if __name__=='__main__':
    unittest.main()