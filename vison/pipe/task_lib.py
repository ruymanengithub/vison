#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 19:11:53 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import os
from pdb import set_trace as stop
import string as st
import numpy as np
from collections import OrderedDict
import copy

from vison.datamodel import HKtools
from vison.pipe import lib as pilib
# END IMPORT




def check_HK(self,HKKeys,reference='command',limits='P',tag='',doReport=False,doLog=True):
    """ """
    
    checkers = dict(command=HKtools.check_HK_vs_command,
                    abs=HKtools.check_HK_abs)
    
    report_HK = checkers[reference](HKKeys,self.dd,limits=limits,elvis=self.elvis) 
    
    if doReport and self.report is not None:
        nicetag = st.replace(tag,' ','\ ')
        msg_HK = ['$\\bf{HK-%s}$, Parameters not passing Tests:\n' % nicetag]
        msg_HK += [st.join(['$\\bf{\\textcolor{red}{%s}}$' % (st.replace(key,'_','\_'),) for key in report_HK if not report_HK[key]],', ')]
        if len(msg_HK[-1]) == 0: msg_HK += ['$\\bf{ALL\ PASSED}$\n']
        self.report.add_Text(st.join(msg_HK,'\n'))
        self.report.add_Text('\\')
        
    if doLog and self.log is not None:
        msg_HK = ['HK-%s, Parameters not passing Tests:\n' % tag]
        msg_HK += [st.join([key for key in report_HK if not report_HK[key]],', ')]
        if len(msg_HK[-1]) == 0: msg_HK += ['ALL PASSED']
        self.log.info(msg_HK)
    
    return report_HK


def add_checkHK_report(self,report_HK,tag):
        
    msg_HK = ['$\\bf{HK-%s}$, Listing parameters not passing Tests:\n' % tag]
    msg_HK += [st.join(['\\bf{%s}' % (st.replace(key,'_','\_'),) for key in report_HK if not report_HK[key]],', ')]
    if len(msg_HK[-1]) == 0: msg_HK += ['$\\bf{ALL\ PASSED}$\n']
    
    if self.report is not None:
        self.report.add_Text(st.join(msg_HK,'\n'))
        self.report.add_Text('\\')
    
    return msg_HK


def convert_compl_to_nesteditemlist(complidict):
    """ """
    
    def traverse_tree(dictionary,nesteditemlist):   
            Nkeys = len(dictionary.keys())
            if Nkeys >5:
                nesteditemlist +=[
                '\item %s' % dictionary.__repr__()
                        ]
                return nesteditemlist
            
            for key,value in dictionary.items():
                if isinstance(value,(dict,OrderedDict)):
                    nesteditemlist += [
                    '\item %s:' % key, 
                    '\\begin{itemize}']
                    nesteditemlist = traverse_tree(value,nesteditemlist)
                    nesteditemlist += ['\\end{itemize}']
                else:
                    nesteditemlist += [
                    '\item %s : %s' % (key,value)        
                            ]
            return nesteditemlist
    nesteditemlist = ['\\begin{itemize}']
    nesteditemlist += traverse_tree(complidict,[])
    nesteditemlist += ['\\end{itemize}']
    return nesteditemlist


def filterexposures(self,structure,explogf,datapath,OBSID_lims,colorblind=False,wavedkeys=[],
                    surrogate=''):
    """Loads a list of Exposure Logs and selects exposures from test 'test'.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    
    """

    # load exposure log(s)
    explog = pilib.loadexplogs(explogf,elvis=self.elvis,addpedigree=True,
                               datapath=datapath)
    
    if len(OBSID_lims)==0:
        OBSID_lims = [explog['ObsID'][0],explog['ObsID'][-1]]
    
    Ncols = structure['Ncols']
    
    if surrogate != '': testkey = surrogate
    else: testkey = self.inputs['test']
    #testkey = structure['col1']['test']
    
    if not colorblind:
    
        Filters = [structure['col%i' % i]['wave'] for i in range(1,Ncols+1)]
        Filter = Filters[0]
        assert np.all(np.array(Filters) == Filter)
        
        selbool = (explog['test'] == testkey) & \
            (explog['ObsID'] >= OBSID_lims[0]) & \
            (explog['ObsID'] <= OBSID_lims[1]) & \
            (explog['wave'] == Filter) # TESTS
    else:
        
        selbool = (explog['test']==testkey) & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) 
    

    explog = explog[np.where(selbool)]
    
    # Assess structure
    
    checkreport = pilib.check_test_structure(explog,structure,CCDs=[1,2,3],
                                           wavedkeys=wavedkeys)
    
    # Labeling of exposures
    
    explog = self.add_labels_to_explog(explog,structure)
    
#    explog['label'] = np.array(['NoneNoneNone']*len(explog))
#    
#    frcounter = 0
#    for ic in range(1,Ncols+1):
#        _frames = structure['col%i' % ic]['frames']
#        #print frcounter,frcounter+_frames*3
#        explog['label'][frcounter:frcounter+_frames*3] = 'col%i' % ic
#        frcounter += _frames*3
    
    return explog, checkreport

def add_labels_to_explog(self,explog,structure):
    """ """
    Ncols = structure['Ncols']
    
    explog['label'] = np.array(['NoneNoneNone']*len(explog))
    
    frcounter = 0
    for ic in range(1,Ncols+1):
        _frames = structure['col%i' % ic]['frames']
        #print frcounter,frcounter+_frames*3
        explog['label'][frcounter:frcounter+_frames*3] = 'col%i' % ic
        frcounter += _frames*3
    
    return explog
    


def addHKPlotsMatrix(self):
    """Adds to self.report a table-figure with HK [self.HKKeys] during test."""
    
    figspath = self.inputs['subpaths']['figs']
    HKKeys = copy.deepcopy(self.HKKeys)
    
    time = self.dd.mx['time'][:,0].copy()
    HKdata = OrderedDict()
    HKdata['time'] = time.copy()
    for key in HKKeys: HKdata[key] = self.dd.mx['HK_%s' % key][:].copy()
    
    HKfigs = []
    
    for key in HKKeys:
        
        HKlims = HKtools.HKlims[self.elvis]['S'][key][1:]
        
        figname = os.path.join(figspath,'HK_%s_%s.eps' % \
                               (key,self.inputs['test']))
        HKtools.doHKSinglePlot(time,HKdata[key],key,ylabel='V',
                               HKlims=HKlims,filename=figname)
        HKfigs.append(figname)
    
    self.report.add_FigsTable(HKfigs,Ncols=3,figswidth=4,
                              caption='Selected HK parameters during test.')
    
    
    
