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

from glob import glob
import os
import numpy as np
# END IMPORT



Quads =['E','F','G','H']
RON = 4.5 # e-
gain = 3.5 # e-/ADU

ccdobj = CCD()

prescan = ccdobj.prescan
overscan = ccdobj.overscan
imgheight = ccdobj.NAXIS2/2
quad_width = ccdobj.NAXIS1/2
imgwidth = quad_width - prescan - overscan

FW = dict(Filter1='570nm',Filter2='650nm',
          Filter3='700nm',Filter4='800nm',
          Filter5='900nm',Filter6='ND')

Point_CooNom = {'CCD1':{'E':{'ALPHA':(prescan+0.2*imgwidth,imgheight*(1.+0.8)),
                              'BRAVO':(prescan+0.8*imgwidth,imgheight*(1.+0.8)),
                             'CHARLIE':(prescan+0.5*imgwidth,imgheight*(1.+0.5)),
                              'DELTA':(prescan+0.2*imgwidth,imgheight*(1.+0.2)),
                              'ECHO':(prescan+0.8*imgwidth,imgheight*(1.+0.2))},
                        'F':{'ALPHA':(quad_width+prescan+0.2*imgwidth,imgheight*(1.+0.8)),
                            'BRAVO':(quad_width+prescan+0.8*imgwidth,imgheight*(1.+0.8)),
                            'CHARLIE':(quad_width+prescan+0.5*imgwidth,imgheight*(1.+0.5)),
                            'DELTA':(quad_width+prescan+0.2*imgwidth,imgheight*(1.+0.2)),
                            'ECHO':(quad_width+prescan+0.8*imgwidth,imgheight*(1.+0.2))},
                        'H':{'ALPHA':(prescan+0.2*imgwidth,imgheight*(0.8)),
                             'BRAVO':(prescan+0.8*imgwidth,imgheight*(0.8)),
                            'CHARLIE':(prescan+0.5*imgwidth,imgheight*(0.5)),
                            'DELTA':(prescan+0.2*imgwidth,imgheight*(0.2)),
                            'ECHO':(prescan+0.8*imgwidth,imgheight*(0.2))},
                        'G':{'ALPHA':(prescan+0.2*imgwidth,imgheight*(0.8)),
                            'BRAVO':(prescan+0.8*imgwidth,imgheight*(0.8)),
                            'CHARLIE':(prescan+0.5*imgwidth,imgheight*(0.5)),
                            'DELTA':(prescan+0.2*imgwidth,imgheight*(0.2)),
                            'ECHO':(prescan+0.8*imgwidth,imgheight*(0.2))}
                            }}

for iCCD in range(2,4):
    Point_CooNom['CCD%i' % iCCD] = dict()
    for Q in Quads:
        Point_CooNom['CCD%i' % iCCD][Q] = copy.deepcopy(Point_CooNom['CCD1'][Q])

# Editions of Point_CooNom: PENDING upon characterization of OGSE




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
    
    TO-DO: account for 3 entries per exposure (3 CCDs). Otherwise no test
           will be compliant.
    
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
            
