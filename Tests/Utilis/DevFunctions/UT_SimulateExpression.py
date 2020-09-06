#!/usr/bin/env python 
"""
@brief: unit testing Simulated gene expression 
@author: Hesham ElAbd
"""
# load the modules:
import unittest
from IPTK.Utils.DevFunctions import simulate_an_expression_table
# define the test case 
class TestCaseSimulateExpression(unittest.TestCase):
    def test_correct_output_shape(self):
        self.assertEqual(
            simulate_an_expression_table().shape[0], 100
        )
        self.assertEqual(
            simulate_an_expression_table(100).shape[1], 3
        )
    def test_assert_correct_rise(self):
        with self.assertRaises(ValueError):
            simulate_an_expression_table(0)

if __name__=='__main__':
    unittest.main()