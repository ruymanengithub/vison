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
from glob import glob
import os
import numpy as np
import datetime
# END IMPORT

dtobj_default = datetime.datetime(1980,2,21,7,0,0) # early riser

Quads =['E','F','G','H']
RON = 4.5 # e-
gain = 3.5 # e-/ADU

ccdobj = CCD()

prescan = ccdobj.prescan
overscan = ccdobj.overscan
imgheight = ccdobj.NAXIS2/2
quad_width = ccdobj.NAXIS1/2
imgwidth = quad_width - prescan - overscan



def loadexplogs(explogfs,elvis='5.8.X',addpedigree=False):
    """loads in memory an explog (text-file) or list of explogs (text-files)."""
    
    if isinstance(explogfs,str):
        explog = ELtools.loadExplog(explogfs,elvis=elvis)
        return explog
    
    elif isinstance(explogfs,list):
        expLogList = []
        for explogf in explogfs:
            expLogList.append(ELtools.loadExpLog(explogf,elvis=elvis))
        explog = ELtools.mergeExpLogs(expLogList,addpedigree)
        return explog
        
def check_test_structure(explog,selbool,structure,CCDs=[1,2,3]):
    """Checks whether a selected number of exposures is consistent with the 
    expected acquisition test structure.
    
    TODO: account for 3 entries per exposure (3 CCDs). Otherwise no test
           will be compliant. [DONE?]
    
    """
    
    isconsistent = True
    
    Ncols = structure['Ncols']
    
    
    for iCCD in CCDs:
        
        cselbool = selbool & (explog['CCD'] == 'CCD%i' % iCCD)
    
        ix0 = 0
        for iC in range(1,Ncols+1):
            
            colname = 'col%i' % iC
        
            expectation = structure[colname]
            N = expectation['N']
            
            ix1 = ix0+N
            
            ixsubsel = np.where(cselbool)[0][ix0:ix1]
            
            logkeys = expectation.keys()
            logkeys.remove('N')
            
            isamatch = np.ones(N,dtype='bool')

            for key in logkeys:
                
                isamatch = isamatch & (explog[key][ixsubsel] == expectation[key])
            
            isconsistent = isconsistent and np.all(isamatch)
            
            ix0 += N
        
    
    return isconsistent        


def addHK(DataDict,HKKeys,elvis='5.8.X'):
    """Adds HK information to a DataDict object."""
    
    
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
                else:
                    print 'More than one HK file for ObsID %i' % ObsID
                    print HKs
                    stop()

            obsids, dtobjs, tdeltasec, HK_keys, HKdata = HKtools.parseHKfiles(HKlist,elvis=elvis)
            
            for HKKey in HKKeys:
                ixkey = HK_keys.index(HKKey)
                DataDict[CCDkey][HKKey] = HKdata[:,0,ixkey]
    
    return DataDict
            
def get_time_tag():
    """ """
    t = datetime.datetime.now()
    s = t.strftime('%Y%m%d_%H%M%S')
    return s
    

def get_dtobj(DT):
    """ """
    
    date = DT[0:DT.index('D')]
    y2d = int(date[4:6])
    if y2d < 20: century = 2000
    else: century = 1900
    dd,MM,yy = int(date[0:2]),int(date[2:4]),y2d+century 
    
    time = DT[DT.index('D')+1:-1]
    
    hh,mm,ss = int(time[0:2]),int(time[2:4]),int(time[4:6])

    dtobj = datetime.datetime(yy,MM,dd,hh,mm,ss)
    
    return dtobj

def save_progress(DataDict,reportobj,DataDictFile,reportobjFile):        
    files.cPickleDumpDictionary(DataDict,DataDictFile)
    files.cPickleDump(reportobj,reportobjFile)


def recover_progress(DataDictFile,reportobjFile):
    DataDict = files.cPickleRead(DataDictFile)
    reportobj = files.cPickleRead(reportobjFile)
    return DataDict,reportobj
