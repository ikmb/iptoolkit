#!/usr/bin/env python 
"""
@author: Hesham ElAbd
@brief: Unit-testing for the peptide class 
"""
# load the models 
import unittest
import numpy as np 
from IPTK.DataStructure.Peptide import Peptide
from IPTK.DataStructure.Protein import Protein
# define some test variables 
PEPTIDE='HMNITALTVSAQSND'
START_INDEX=25
END_INDEX=40
PAR_PROT_NAME="P11446"
PAR_PROT_SEQ="MLNTLIVGASGYAGAELVTYVNRHPHMNITALTVSAQSNDAGKLISDLHPQLKGIVDLPLQPMSDISEFSPGVDVVFLATAHEVSHDLAPQFLEAGCVVFDLSGAFRVNDATFYEKYYGFTHQYPELLEQAAYGLAEWCGNKLKEANLIAVPGCYPTAAQLALKPLIDADLLDLNQWPVINATSGVSGAGRKAAISNSFCEVSLQPYGVFTHRHQPEIATHLGADVIFTPHLGNFPRGILETITCRLKSGVTQAQVAQVLQQAYAHKPLVRLYDKGVPALKNVVGLPFCDIGFAVQGEHL"
PAR_ORG='E.coli'
NUMBER_OF_REPS=100
PEPTIDES_LENGTH=list(range(9,25,1))
# map the peptide to the protein to use the test: ==>
back_bone=np.zeros(shape=(1,len(PAR_PROT_SEQ)))
back_bone[0,START_INDEX:END_INDEX]=1
class TestCasePeptide(unittest.TestCase):
    def setUp(self):
        self._peptide=Peptide(PEPTIDE)
        self._test_parent_protein=Protein(PAR_PROT_NAME,PAR_PROT_SEQ)
    def test_get_length(self):
        self.assertEqual(self._peptide.get_length(),len(PEPTIDE))
    
    def test_len__(self):
        self.assertEqual(len(self._peptide), len(PEPTIDE))
    
    def test_get_peptide_seq(self): 
        self.assertEqual(self._peptide.get_peptide_seq(),PEPTIDE) 
    
    def test_add_parent_protein(self):
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        self.assertEqual(self._peptide.get_number_parent_protein(),1)
        ## assert the exception is thrown correctly 
        with self.assertRaises(ValueError):
            self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX-1, END_INDEX)
        # check repeated calling to the function 
        for _ in range(NUMBER_OF_REPS):
            self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        self.assertEqual(self._peptide.get_number_parent_protein(),1)
    
    def test_get_flanked_region(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # generate 3 test cases and assert they have the coorect size 
        short_region=self._peptide.get_flanked_peptide(5)
        self.assertEqual(len(short_region),1)
        long_flanked_region=self._peptide.get_flanked_peptide(15)
        self.assertEqual(len(long_flanked_region),1)
        ultra_long_flanked_region=self._peptide.get_flanked_peptide(60)
        self.assertEqual(len(ultra_long_flanked_region),1)
        # assert that the generated flanked region has the correct content 
        self.assertEqual(short_region[0],PAR_PROT_SEQ[START_INDEX-5:END_INDEX+5])
        # assert for the longer protein 
        self.assertEqual(long_flanked_region[0],PAR_PROT_SEQ[START_INDEX-15:END_INDEX+15])
        # assert for ultra-long peptide where the flanked sequence is bigger than N-terminal side of the peptide 
        self.assertEqual(ultra_long_flanked_region[0],PAR_PROT_SEQ[:END_INDEX+60]) 
    
    def test_map_to_parent_protein(self): 
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # get the parent protein and check the correct size of the output class 
        mapped_to_parent=self._peptide.map_to_parent_protein()
        self.assertEqual(len(mapped_to_parent),1)
        # assert the correctness of the output
        self.assertEqual(np.mean(mapped_to_parent[0].ravel()==back_bone.ravel()),True, )
    
    def test_get_non_present_peptides(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # assert that the returned number of peptides is correct 
        self.assertEqual(len(self._peptide.get_non_presented_peptides(15)),1)
        # assert the correctness of the returned values
        peptide_seq=self._peptide.get_peptide_seq()
        for _ in range(NUMBER_OF_REPS):
            self.assertFalse(self._peptide.get_non_presented_peptides(15)[0]==peptide_seq)
            self.assertEqual(len(self._peptide.get_non_presented_peptides(15)[0]),len(peptide_seq))
    
    def test_parent_proteins(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # check that the number of elements returned by the function 
        self.assertEqual(len(self._peptide.get_parent_proteins()),1)
        # asset the name is correct 
        self.assertEqual(list(self._peptide.get_parent_proteins())[0],PAR_PROT_NAME)
    
    def test_is_child_of(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # check that function return True
        self.assertTrue(self._peptide.is_child_of(PAR_PROT_NAME))

    def test_get_number_of_parents(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # assert that the function return the correct number of parents 
        self.assertTrue(self._peptide.get_number_of_parents()==1)
    
    def test_get_position_in_parent(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # get the start and end position 
        start, end = self._peptide.get_pos_in_parent(PAR_PROT_NAME)
        self.assertTrue(START_INDEX==start)
        self.assertTrue(END_INDEX==end)
    
    def test_get_parent(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # get the parent 
        parent_in_peptide=self._peptide.get_parent(PAR_PROT_NAME)
        self.assertEqual(len(parent_in_peptide),len(PAR_PROT_SEQ))
        self.assertEqual(str(parent_in_peptide),PAR_PROT_SEQ)
    
    def test_get_n_terminal_flank_seq(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # get the n_terminal flanking region 
        flank_region=self._peptide.get_n_terminal_flank_seq(10) # get 10 amino acids from the peptide N-terminal
        # assert that the function return a list of one elements 
        self.assertEqual(len(flank_region),1)
        # assert the correctness of the returned sequence 
        self.assertEqual(flank_region[0],PAR_PROT_SEQ[START_INDEX-10:START_INDEX])
        # assert the function behavior with out of index length
        self.assertEqual(self._peptide.get_n_terminal_flank_seq(30)[0],PAR_PROT_SEQ[0:START_INDEX])

    def test_get_c_terminal_flank_seq(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        # get the c_terminal flanking region 
        flank_region=self._peptide.get_c_terminal_flank_seq(10) # get 10 amino acids from the peptide N-terminal
        # assert that the function return a list of one elements 
        self.assertEqual(len(flank_region),1)
        # assert the correctness of the returned sequence 
        self.assertEqual(flank_region[0],PAR_PROT_SEQ[END_INDEX:END_INDEX+10])
    
    def test_add_and_get_org(self):
        # add a parent protein 
        self._peptide.add_parent_protein(self._test_parent_protein, START_INDEX, END_INDEX)
        self._peptide.add_org_2_parent(PAR_PROT_NAME,PAR_ORG)
        self.assertEqual(len(self._peptide.get_parents_org()),1)
        self.assertEqual(self._peptide.get_parents_org()[0],PAR_ORG)

    def test_getitem__(self):
        # test that the three amino acids match 
        self.assertEqual(self._peptide[6:9],PEPTIDE[6:9])

    def test_str__(self):
        self.assertEqual(str(self._peptide),PEPTIDE)
    
    def test_repr__(self):
        self.assertEqual(repr(self._peptide),'A Peptide instance with 0 parent protein')

if __name__=='__main__':
    unittest.main()