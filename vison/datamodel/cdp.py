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
import os

from vison import __version__
from vison.support.files import cPickleDumpDictionary, cPickleRead
# END IMPORT


def loadCDPfromPickle(pickf):
    """ """
    cdp = cPickleRead(pickf)
    return cdp

class CDP(object):
    """ """

    rootname = 'Unknown'
    path = ''
    header = OrderedDict()
    meta = OrderedDict()
    data = None
    
    #def __init__(self,*args,**kwargs):
    #    """ """
    def __init__(self,*args,**kwargs):
        """ """    
        #inputs = dict(ID=None,BLOCKID=None,CHAMBER=None)
        #inputs.update(args)
        #inputs.update(kwargs)
        
        #self.ID = inputs['ID']
        #self.BLOCKID = inputs['BLOCKID']
        #self.CHAMBER = inputs['CHAMBER']
        #self.vison = __version__
        
        if 'header' in args: self.header.update(['header'])
        elif 'header' in kwargs: self.header.update(kwargs['header'])
        
        if 'meta' in args: self.meta.update(['meta'])
        elif 'meta' in kwargs: self.meta.update(kwargs['meta'])
        
        
        
    def savetopickle(self,pickf=''):
        """ """
        if pickf == '':
            outf = os.path.join(self.path,'%s.pick' % self.rootname)
        else:
            outf = pickf
        #outdict = OrderedDict(Header=self.header,Data=self.data)
        cPickleDumpDictionary(self,outf)
    
    
    def savehardcopy(self,filef=''):
        """ """
        raise NotImplementedError('Subclass implements abstract method (if needed).') 