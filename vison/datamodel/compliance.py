#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Some functions to present COMPLIANCE MATRICES.

Created on Mon Apr  9 17:32:03 2018

:author: raf
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import pandas as pd
import collections
import string as st
import copy
import numpy as np

import unittest

from vison.support import vjson
from vison.support.utils import countdictlayers
# END IMPORT
        
        



def convert_compl_to_nesteditemlist(complidict):
    """ """    

    def traverse_tree(dictionary, nesteditemlist):
        Nkeys = len(dictionary.keys())
        if Nkeys > 5:
            nesteditemlist += [
                '\item %s' % dictionary.__repr__()
            ]
            return nesteditemlist

        for key, value in dictionary.items():
            if isinstance(value, (dict, OrderedDict)):
                nesteditemlist += [
                    '\item %s:' % key,
                    '\\begin{itemize}']
                nesteditemlist = traverse_tree(value, nesteditemlist)
                nesteditemlist += ['\\end{itemize}']
            else:
                nesteditemlist += [
                    '\item %s : %s' % (key, value)
                ]
        return nesteditemlist

    tex = ['\\begingroup'] +\
        ['\\scriptsize'] +\
        ['\\begin{itemize}'] +\
        traverse_tree(complidict, []) +\
        ['\\end{itemize}']
    
    return tex


def removescalars_from_dict(indict):
    """ """

    def traverse_tree(indict):
        for key, value in indict.iteritems():
            if isinstance(value, (dict, OrderedDict)):
                indict[key] = traverse_tree(value)
            elif not isinstance(value, collections.Sequence):
                indict[key] = [value]
        return indict
    
    outdict = traverse_tree(copy.deepcopy(indict))

    return outdict


def gen_compliance_tex(indict):
    """ """
    
    complidict = copy.deepcopy(indict)
    
    level = countdictlayers(indict)
    
    if level < 3:
        if level == 1: complidict = removescalars_from_dict(complidict)
        df = pd.DataFrame.from_dict(complidict)
        tex = df.to_latex(multicolumn=True, multirow=True,
                          longtable=True, index=level>1)
        tex = st.split(tex, '\n')
    elif level == 3:
        keys = []
        frames = []
        for key, d in complidict.iteritems():
            keys.append(key)
            frames.append(pd.DataFrame.from_dict(d))
        series = pd.concat(frames, keys=keys)
        
        tex = series.to_latex(multicolumn=True, multirow=True,
                          longtable=True, index=True)
        tex = st.split(tex, '\n')
    else:
        tex = convert_compl_to_nesteditemlist(complidict)
        
    return tex

class ComplianceMX(OrderedDict):
    
    def __init__(self, *args, **kwargs):
        """ """
        super(ComplianceMX,self).__init__()
        
    def check_stat(self,inparr):
        raise NotImplementedError("Subclass must implement abstract method")
    
    def get_compliance_tex(self):
        return gen_compliance_tex(dict(self))
    
    def get_compliance_txt(self):
        return vjson.dumps_to_json(self)
        #return self.__str__()
    
    def IsCompliant(self):
        
        def traverse_tree(dictionary, isOK):
            
            for key, value in dictionary.items():
                #print 'Upper: %s' % key
                if isinstance(value, (dict, OrderedDict)):
                    #print key,value
                    isOK = isOK and traverse_tree(value, isOK)
                else:
                    #print key,value
                    isOK = isOK and value[0]
            return isOK
        
        
        isOK = traverse_tree(self, True)
        
        return isOK

class ComplianceMX_CCDQ(ComplianceMX):
    
    def __init__(self,CCDs=None,Qs=None,CCDQlims=None):
        """ """
        super(ComplianceMX_CCDQ,self).__init__()
        if CCDs is None:
            self.CCDs = ['CCD1','CCD2','CCD3']
        else:
            self.CCDs = CCDs
        if Qs is None:
            self.Qs = ['E','F','G','H']
        else:
            self.Qs = Qs
            
        for CCD in self.CCDs:
            self[CCD] = OrderedDict()
            #for Q in self.Qs:
            #    self[CCD][Q] = OrderedDict()
        
        if CCDQlims is None:        
            self.CCDQlims = OrderedDict()
        else:
            self.CCDQlims = CCDQlims
        
        
    def check_stat(self,inparr):
        """ """
        for iCCD, CCDkey in enumerate(self.CCDs):
            for jQ, Q in enumerate(self.Qs):
                if isinstance(self.CCDQlims[CCDkey],dict):
                    _lims = self.CCDQlims[CCDkey][Q]
                elif isinstance(self.CCDQlims[CCDkey],list):
                    _lims = self.CCDQlims[CCDkey]
                test = (np.isnan(inparr[:, iCCD, jQ, ...]) |
                        (inparr[:, iCCD, jQ, ...] <= _lims[0]) | (inparr[:, iCCD, jQ, ...] >= _lims[1]))
                try: avvalue = float(np.nanmean(inparr[:,iCCD, jQ, ...]))
                except: avvalue = np.nan
                self[CCDkey][Q] = [not np.any(test).sum(), avvalue, _lims]


class ComplianceMX_CCD(ComplianceMX):
    
    def __init__(self,CCDs=None, CCDlims=None):
        """ """
        super(ComplianceMX_CCD,self).__init__()
        if CCDs is None:
            self.CCDs = ['CCD1','CCD2','CCD3']
        else:
            self.CCDs = CCDs
        
        if CCDlims is None:        
            self.CCDlims = OrderedDict()
        else:
            self.CCDlims = CCDlims
        
        
    def check_stat(self,inparr):
        """ """
        for iCCD, CCDkey in enumerate(self.CCDs):
            
            _lims = self.CCDQlims[CCDkey]
            test = (np.isnan(inparr[:, iCCD, ...]) |\
                    (inparr[:, iCCD, ...] <= _lims[0]) |\
                    (inparr[:, iCCD, ...] >= _lims[1]))
            
            try: avvalue = float(np.nanmean(inparr[:,iCCD, ...]))
            except: avvalue = np.nan
            
            self[CCDkey] = [not np.any(test).sum(), avvalue, _lims]


class ComplianceMX_CCDCol(ComplianceMX):
    
    def __init__(self, colnames, indexer, CCDs=None, lims=None):
        """ """
        super(ComplianceMX_CCD,self).__init__()
        
        self.colnames = colnames
        self.indexer = indexer
        
        if CCDs is None:
            self.CCDs = ['CCD1','CCD2','CCD3']
        else:
            self.CCDs = CCDs
        
        for CCD in self.CCDs:
            self[CCD] = OrderedDict()
        
        if lims is None:        
            self.lims = OrderedDict()
        else:
            self.lims = lims
        
        
    def check_stat(self,inparr):
        """ """
        for iCCD, CCDkey in enumerate(self.CCDs):
            for jcol, colname in enumerate(self.colnames):
                _lims = self.lims[CCDkey][colname]
                ixsel = np.where(self.indexer == colname)

                testv = (np.isnan(inparr[ixsel, iCCD, ...]) |
                        (inparr[ixsel, iCCD, ...] <= _lims[0]) | (inparr[ixsel, iCCD, ...] >= _lims[1]))
                
                try: avvalue = float(np.nanmean(inparr[ixsel,iCCD, ...]))
                except: avvalue = np.nan
                
                test = not (
                    np.any(testv, axis=(0, 1)).sum() | (ixsel[0].shape[0] == 0))
                
                self[CCDkey][colname] = [test, avvalue, _lims]

class ComplianceMX_CCDQCol(ComplianceMX):
    
    def __init__(self, colnames, indexer, CCDs=None, Qs=None, lims=None):
        """ """
        super(ComplianceMX_CCD,self).__init__()
        
        self.colnames = colnames
        self.indexer = indexer
        
        if CCDs is None:
            self.CCDs = ['CCD1','CCD2','CCD3']
        else:
            self.CCDs = CCDs
            
        if Qs is None:
            self.Qs = ['E','F','G','H']
        else:
            self.Qs = Qs
        
        for CCD in self.CCDs:
            self[CCD] = OrderedDict()
            for Q in self.Qs:
                self[CCD][Q] = OrderedDict()
            
        
        if lims is None:        
            self.lims = OrderedDict()
        else:
            self.lims = lims
        
        
    def check_stat(self,inparr):
        """ """
        for iCCD, CCDkey in enumerate(self.CCDs):
            for jQ, Q in enumerate(self.Qs):
                for kcol, colname in enumerate(self.colnames):
                    _lims = self.lims[CCDkey][colname]
                    ixsel = np.where(self.indexer == colname)
    
                    testv = (np.isnan(inparr[ixsel, iCCD, jQ, ...]) |\
                            (inparr[ixsel, iCCD, jQ, ...] <= _lims[0]) | \
                            (inparr[ixsel, iCCD, jQ, ...] >= _lims[1]))
                    
                    try: avvalue = float(np.nanmean(inparr[ixsel,iCCD, jQ, ...]))
                    except: avvalue = np.nan
                    
                    test = not (
                        np.any(testv, axis=(0, 1)).sum() | (ixsel[0].shape[0] == 0))
                    
                    self[CCDkey][Q][colname] = [test, avvalue, _lims]



class TestComplianceClass(unittest.TestCase):
    
    def setUp(self):
        self.comp = ComplianceMX()
    
    def test0(self):
        stop()


def testit():
    unittest.main()
    

if __name__ == '__main__':
    
    testit()
