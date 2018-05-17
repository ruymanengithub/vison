#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 11:12:02 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import unittest

from vison.point.lib import Point_CooNom
# END IMPORT



class StarTracker(object):
    
    def __init__(self):
        pass




class Test(unittest.TestCase):
    """
    Unit tests for the point-source location tools.
    """
    def setUp(self):
        
        self.startracker = StarTracker()
        

    def test_something(self):
        """
        
        :return: None
        """
        
        self.assertTrue(True)


if __name__ == '__main__':
    
    suite = unittest.TestLoader().loadTestsFromTestCase(Test)
    unittest.TextTestRunner(verbosity=3).run(suite)
    
    