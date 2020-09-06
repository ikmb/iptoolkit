#!/usr/bin/env python 
"""
@brief: UT for the proband case 
@author: Hesham ElAbd
"""
# load the modules 
import unittest
from IPTK.DataStructure.Proband import Proband
# build the unit test 
class TestCaseProband(unittest.TestCase):
    def setUp(self):
        self._default_proband=Proband(name='PROBAND', status='H')

    def test_get_name(self):
        self.assertEqual(self._default_proband.get_name(), 'PROBAND')

    def test_update_info(self):
        self._default_proband.update_info(name='UPDATED_NAME')
        self.assertEqual(self._default_proband.get_name(),'UPDATED_NAME')
    
    def test_get_meta_info(self):
        test_proband=Proband(name='Test_NAME',age=27).get_meta_data()
        self.assertEqual(test_proband['name'],'Test_NAME')
        self.assertEqual(test_proband['age'],27)

if __name__=='__main__':
    unittest.main()