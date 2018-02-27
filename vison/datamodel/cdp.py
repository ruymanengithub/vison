#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Euclid Ground Calibration Campaign


Classes to store Calibration Data Products.

Created on Tue Feb 27 10:58:42 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict

from vison import __version__
from vison.support.files import cPickleDumpDictionary, cPickleRead
# END IMPORT


def loadCDPfromPickle(pickf):
    """ """
    cdp = cPickleRead(pickf)
    return cdp

class CDP(OrderedDict):
    """ """
    
    def __init__(self,data=OrderedDict(),meta=OrderedDict(),ID=None,BLOCKID=None,CHAMBER=None):
        """ """
        
        self.ID = ID
        self.BLOCKID = BLOCKID
        self.CHAMBER = CHAMBER
        self.vison = __version__
        
        self.data = data.copy()
        self.meta = meta.copy()
        
        
    def savetopickle(self,pickf):
        """ """        
        cPickleDumpDictionary(self,pickf)
    
    
    