#!/usr/bin/env python
"""
@brief: A script to perform unit-testing on the SeqDB classes 
@author: Hesham ElAbd
"""
# load the modules:
from IPTK.DataStructure.Database import SeqDB
import unittest
# define the testing class 
class TestDatabaseClass(unittest.TestCase):
    def setUp(self):
        self._database=SeqDB('../assets/e_coli.fasta')

    def test_len_function(self):
        self.assertEqual(len(self._database),4391)
    
    def test_str_method(self):
        self.assertEqual(str(self._database), 'A sequence database with 4391 sequence')
    
    def test_getitem_method(self):
        self.assertEqual(self._database['P11446'][:10],
        'MLNTLIVGAS'
        )
    def test_correct_exception_handle(self):
        with self.assertRaises(KeyError):
            self._database['NO_TEST_SEQUENCE']

    def test_has_sequecne(self):
        self.assertTrue(self._database.has_sequence('P02924'))
        self.assertFalse(self._database.has_sequence('NO_IN_SEQ_DB'))

# assert that the name is correct 
if __name__=='__main__':
    unittest.main()
