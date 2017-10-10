# -*- coding: utf-8 -*-
"""

Support Script with Variables and Functions used for FM-Calib. analysis

Created on Wed Nov 30 11:11:27 2016


:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from pdb import set_trace as stop

import copy
from vison.datamodel.ccd import CCD
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.support import files
from vison.support import time as vistime
from glob import glob
import os
import numpy as np
# END IMPORT


elvis = '6.3.0'
Quads =['E','F','G','H']
RON = 4.5 # e-
gain = 3.5 # e-/ADU

ccdobj = CCD()

prescan = ccdobj.prescan
overscan = ccdobj.overscan
imgheight = ccdobj.NAXIS2/2
quad_width = ccdobj.NAXIS1/2
imgwidth = quad_width - prescan - overscan

sumwell = dict(fwd_bas=[9.425,4.825],fwd_e2v=[7.475,6.55],
               rwd_bas_v=[6.5,0.175],rwd_bas_s=[6.5,0.175],
               rwd_bas_vs=[6.5,0.175])

def loadexplogs(explogfs,elvis='6.3.0',addpedigree=False,datapath=None):
    """loads in memory an explog (text-file) or list of explogs (text-files)."""
    
    if isinstance(explogfs,str):
        explog = ELtools.loadExpLog(explogfs,elvis=elvis)
        
    
    elif isinstance(explogfs,list):
        expLogList = []
        for explogf in explogfs:
            expLogList.append(ELtools.loadExpLog(explogf,elvis=elvis))
        explog = ELtools.mergeExpLogs(expLogList,addpedigree)
    
    # add datapath(s)
    # The datapath becomes another column in DataDict. This helps dealing
    # with tests that run overnight and for which the input data is in several
    # date-folders.
    
    if datapath is not None:
        
        if isinstance(datapath,list):
            assert len(datapath) == len(explogfs)
            
            longestdatapathname = max([len(item) for item in datapath])
            explog['datapath'] = np.zeros(len(explog),dtype='S%i' % longestdatapathname)
            explognumber=explog['explognumber']
            for idata in range(len(datapath)):
                explog['datapath'][explognumber == idata] = datapath[idata]
        else:
            explog['datapath'] = np.array([datapath]*len(explog))
            
    
    return explog

        
def check_test_structure(explog,structure,CCDs=[1,2,3],selbool=True,wavedkeys=[]):
    """Checks whether a selected number of exposures is consistent with the 
    expected acquisition test structure.
    
    
    """
    
    from vison.datamodel.generator import _update_fromscript
    
    isconsistent = True
    
    Ncols = structure['Ncols']
    
    elogkeys = explog.keys()
    
    dummyrowdict = dict([(key,None) for key in elogkeys])
    
    failedkeys = []
    failedcols = []
    
    for iCCD in CCDs:
        
        cselbool = selbool & (explog['CCD'] == 'CCD%i' % iCCD)
    
        ix0 = 0
        for iC in range(1,Ncols+1):
            
            
            colname = 'col%i' % iC
        
            expectation = structure[colname]
            
            frames = expectation['frames']
            
            if frames >0:
                ix1 = ix0+frames
            else:
                ix1 = None

            if (ix1 is not None) and (ix1-ix0 != frames):
                failedcols.append(iC)
                continue

            
            try:
                ixsubsel = np.where(cselbool)[0][ix0:ix1]
            except IndexError:
                failedcols.append(iC)
                continue
            
            rowdict = _update_fromscript(dummyrowdict,expectation)
            
            for key in rowdict:
                val = rowdict[key]
                if val is None or key in wavedkeys: continue
                checksout = np.all(explog[key][ixsubsel] == val)                
                isconsistent &= checksout
                if ~checksout:
                    failedkeys.append(key)
                    
            ix0 += frames
    
    failedkeys = np.unique(failedkeys)
    failedcols = np.unique(failedcols)
    
    report = dict(checksout=isconsistent,failedkeys=failedkeys,failedcols=failedcols)
    
    return report
            

def DataDict_builder(explog,inputs,structure):
    """ """
    
    # Build DataDict - And Labelling
    
    DataDict = dict(meta = dict(inputs=inputs,structure=structure))
    
    for CCDindex in [1,2,3]:
        
        CCDkey = 'CCD%i' % CCDindex
        
        DataDict[CCDkey] = dict()
        
        CCDselbool = explog['CCD'] == CCDkey
        
        if len(np.where(CCDselbool)[0]) == 0:
            continue
                
        for key in explog.colnames:
            DataDict[CCDkey][key] = explog[key][CCDselbool].data.copy()
        
        
    
    return DataDict
    


def addHK(DataDict,HKKeys,elvis='5.8.X'):
    """Adds HK information to a DataDict object."""
    
    
    if len(HKKeys) == 0:
        return DataDict
    
    for ixCCD in [1,2,3]:
        CCDkey = 'CCD%i' % ixCCD
        
        if CCDkey in DataDict:
            ObsIDs = DataDict[CCDkey]['ObsID'].copy()            
            datapaths = DataDict[CCDkey]['datapath'].copy()
            
            HKlist = []
            
            for iOBS,ObsID in enumerate(ObsIDs):
                tmp = os.path.join(datapaths[iOBS],'HK_%s_*_ROE1.txt' % ObsID)
                
                HKs = glob(tmp)
                
                if len(HKs) == 1:
                    HKlist.append(HKs[0])
                elif len(HKs)>1:
                    print 'More than one HK file for ObsID %i' % ObsID
                    print HKs
                elif len(HKs) ==0:
                    print 'HK file for ObsID %i not found' % ObsID
            

            obsids, dtobjs, tdeltasec, readHKKeys, HKdata = HKtools.parseHKfiles(HKlist,elvis=elvis)
            
            for HKKey in HKKeys:
                ixkey = readHKKeys.index(HKKey)
                DataDict[CCDkey][HKKey] = HKdata[:,0,ixkey]
    
    return DataDict


def coarsefindTestinExpLog(explog,testkey,Nframes):
    """Checks whether the 'Nframes' of test with key 'testkey' are in 'explog'.
    :param explog: astropy table, exposure log object.
    :param testkey: char, 
    :param Nframes: number of frames expected from test.
    :return: bool, test was acquired or not.
    
    """
    
    indices = np.where(explog['TEST'] == testkey)
    ccds = explog['CCD'][indices]
    uccd = np.unique(ccds)
    nccds = len(uccd)
    Nobsids = len(indices[0])/nccds
    
    wasacquired = Nobsids == Nframes    
    return wasacquired


def save_progress(DataDict,reportobj,DataDictFile,reportobjFile):        
    files.cPickleDumpDictionary(DataDict,DataDictFile)
    files.cPickleDump(reportobj,reportobjFile)


def recover_progress(DataDictFile,reportobjFile):
    DataDict = files.cPickleRead(DataDictFile)
    reportobj = files.cPickleRead(reportobjFile)
    return DataDict,reportobj
