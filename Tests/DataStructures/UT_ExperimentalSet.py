#!/usr/bin/env python 
"""
@brief: Testing Unit-testing function 
@details: the unit-testing for ExperimentalSet 
"""
# import the modules 
import unittest
from typing import List, Dict
from IPTK.DataStructure.ExperimentalSet import ExperimentSet
from IPTK.DataStructure.Experiment import Experiment
from IPTK.Utils.UtilityFunction import generate_random_name
from IPTK.Utils.DevFunctions import simulate_random_experiment
import random
import numpy as np
# define some Test Variables 
NUMBER_EXPS: int =10
NAMES: List[str]=['TEST_EXP_'+str(idx) for idx in range(NUMBER_EXPS)]
ALLELEs: List[List[str]] = [
    ['DRB1*15:01','DRB1*13:01'],
    ['DRB1*15:01','DRB1*15:02'],
    ['DRB1*13:02','DRB1*08:02'],
    ['DRB1*15:01','DRB1*08:02'],
    ['DRB1*13:01','DRB1*14:02'],
    ['DRB1*11:01','DRB1*14:02'],
    ['DRB1*15:01','DRB1*14:02'],
    ['DRB1*07:01','DRB1*14:02'],
    ['DRB1*08:04','DRB1*05:02'],
    ['DRB1*11:01','DRB1*14:02']
]
TISSUES: List[str]=['Lung','Lung','Lung','gut','gut','gut','liver','liver','liver','liver']
DATA_PATH: str= '../assets/e_coli.fasta'
# define the class 
class TestCaseExperimentSet(unittest.TestCase):
    def setUp(self):
        self._exps=ExperimentSet()
        for idx in range(NUMBER_EXPS):
            temp={NAMES[idx]:simulate_random_experiment(
                ALLELEs[idx],DATA_PATH,TISSUES[idx])}
            self._exps.add_experiment(**temp)
    
    def test_len__(self):
        self.assertEqual(len(self._exps),NUMBER_EXPS)
    
    def test_get_experimental_name(self):
        names: List[str] = self._exps.get_experimental_names()
        for name in names:
            self.assertIn(name,NAMES)
    
    def test_getitem__(self):
        for name in NAMES:
            self._exps[name]

    def test_get_tissue_counts(self):
        counts: Dict[str, int]= self._exps.get_tissue_counts()
        self.assertEqual(len(counts),3)
        self.assertEqual(counts['liver'],4)
        self.assertEqual(counts['gut'],3)
        self.assertEqual(counts['Lung'],3)

    def test_get_allele_count(self):
        counts:  Dict[str, int]= self._exps.get_allele_count()
        # assert the correctness of the counts 
        self.assertEqual(counts['DRB1*15:01'],4)
        self.assertEqual(counts['DRB1*13:01'],2)
        self.assertEqual(counts['DRB1*05:02'],1)
        self.assertEqual(counts['DRB1*07:01'],1)
    
    def test_get_proband_count(self):
        counts: Dict[str, int]= self._exps.get_proband_count()
        
    def test_group_by_tissue_type(self):
        group_by_tissue: Dict[str,ExperimentSet] = self._exps.group_by_tissue()
        # define the assert in 
        self.assertIn('Lung', group_by_tissue.keys())
        self.assertIn('gut',group_by_tissue.keys())
        self.assertIn('liver',group_by_tissue.keys())
        # define the list inclusion 
        self.assertEqual(len(group_by_tissue['Lung']),3)
        self.assertEqual(len(group_by_tissue['gut']),3)
        self.assertEqual(len(group_by_tissue['liver']),4)
    
    def test_group_by_proband_name(self):
        group_by_proband: Dict[str, ExperimentSet]= self._exps.group_by_proband()
        for name in group_by_proband.keys():
            self.assertEqual(len(group_by_proband[name]),1)
    ## DEFINE SOME RUN & DEBUG TESTS
    ## Add more correctness measures later 
    def test_get_unique_peptides(self):
        peptide: List[str]= self._exps.get_unique_peptides()
        print(len(peptide)) 

    def test_get_unique_proteins(self):
        proteins: List[str]= self._exps.get_unique_proteins()
        print(len(proteins))
    
    def test_get_peptides_present_in_all(self):
        present_in_all: List[str]= self._exps.get_peptides_present_in_all()
    
    def test_get_protein_present_in_all(self):
        proteins_in_all: List[str]= self._exps.get_protein_present_in_all()

    def test_compute_peptide_overlap_matrix(self):
        overLap_matrix: np.ndarray = self._exps.compute_peptide_overlap_matrix()
    
    def test_compute_protein_overlap_matrix(self):
        overLap_matrix: np.ndarray = self._exps.compute_protein_overlap_matrix()
    
    def test_compute_peptide_representation_count(self):
        counts: Dict[str, int] = self._exps.compute_peptide_representation_count()

    def test_compute_protein_coverage_over_the_set(self):
        counts: Dict[str, int]= self._exps.compute_protein_coverage_over_the_set()


if __name__=='__main__':
    unittest.main()