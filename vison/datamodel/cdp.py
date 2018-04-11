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
from vison.support.excel import ReportXL
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
            pickf = os.path.join(self.path,'%s.pick' % self.rootname)
        
        #outdict = OrderedDict(Header=self.header,Data=self.data)
        cPickleDumpDictionary(self,pickf)
    
    
    def savehardcopy(self,filef=''):
        """ """
        raise NotImplementedError('Subclass implements abstract method (if needed).') 
        

class Tables_CDP(CDP):
    """ """
    
    figs = OrderedDict()
    data = OrderedDict()
    
    def __init__(self,*args,**kwargs):
        """ """
        super(Tables_CDP,self).__init__(*args,**kwargs)        
        self.xlsx = None
        
    
    def ingest_inputs(self,data,meta=OrderedDict(),header=OrderedDict(),figs=OrderedDict()):
        """ """
        
        self.header.update(header)
        self.meta.update(meta)
        self.data.update(data)
        self.figs.update(figs)
        
    
    def init_workbook(self):
        """ """
        self.report = ReportXL(OrderedDict())
        self.report.wb.create_sheet('Header',0)
        self.report.wb.create_sheet('Meta',1)
        
        for i,k in enumerate(self.iterkeys()):
            self.report.wb.create_sheet(k,i)
        
    def fill_Header(self,title=''):
        """ """
        headdict = self.header.copy()
        self.report.dict_to_sheet(headdict,'Header',title=title)
        
    
    def fill_Meta(self):
        """ """
        metadict = self.meta.copy()
        self.report.dict_to_sheet(metadict,'Meta',title='MetaInfo')
    
    
    def fill_Sheet(self,sheet):
        """ """
        df = self.data[sheet]
        self.report.df_to_sheet(df,sheet,index=False,header=False)
        
    
    def fill_allSheets(self):
        """ """
        for sheet in self.data.keys():
            self.fill_Sheet(sheet)
    
    
    def savehardcopy(self,filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path,'%s.xlsx' % self.rootname)
        
        self.report.save(filef)
        
    
    