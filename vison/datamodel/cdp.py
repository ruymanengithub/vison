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
import  copy

from vison import __version__
from vison.support.files import cPickleDumpDictionary, cPickleRead
# END IMPORT


def loadCDPfromPickle(pickf):
    """ """
    cdp = cPickleRead(pickf)
    return cdp

class CDP(object):
    """ """
    
    #def __init__(self,data=OrderedDict(),meta=OrderedDict(),ID=None,BLOCKID=None,CHAMBER=None):
    def __init__(self,*args,**kwargs):
        """ """
        
        inputs = dict(ID=None,BLOCKID=None,CHAMBER=None)
        inputs.update(args)
        inputs.update(kwargs)
        
        self.ID = inputs['ID']
        self.BLOCKID = inputs['BLOCKID']
        self.CHAMBER = inputs['CHAMBER']
        self.vison = __version__
        
        #self.data = copy.deepcopy(data)
        #self.meta = copy.deepcopy(meta)
        
        
    def savetopickle(self,pickf):
        """ """
        cPickleDumpDictionary(self,pickf)
    
    
    