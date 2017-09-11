#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BIAS01

Bias-structure/RON analysis script

Created on Tue Aug 29 16:53:40 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
import datetime
from copy import deepcopy


#from vison.support.report import Report


# END IMPORT

isthere = os.path.exists


HKKeys = []


BIAS01_commvalues = dict(program='CALCAMP',test='BIAS01',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,#sinvflushp=500,
  chinj=0,
  tpump=0,
  motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='BIAS')
  


def build_BIAS01_scriptdict(N,diffvalues=dict(),elvis='6.0.0'):
    """Builds BIAS01 script structure dictionary.
    
    :param N: integer, number of frames to acquire.
    :param diffvalues: dict, opt, differential values.
    :param elvis: char, ELVIS version.
    
    """
    
    BIAS01_sdict = dict(col1=dict(frames=N,exptime=0))

    Ncols = len(BIAS01_sdict.keys())    
    BIAS01_sdict['Ncols'] = Ncols
    
                
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(BIAS01_commvalues)
    
    BIAS01_sdict = sc.update_structdict(BIAS01_sdict,commvalues,diffvalues)
    
    return BIAS01_sdict



def filterexposures(structure,explogf,datapath,OBSID_lims,
                                     elvis):
    """ """
    
    DataDict = {}
    isconsistent = False

    #The filtering takes into account an expected structure for the 
    #acquisition script.

    #The datapath becomes another column in DataDict. This helps dealing
    #with tests that run overnight and for which the input data is in several
    #date-folders.

    
    
    # load exposure log(s)
    
    explog = pilib.loadexplogs(explogf,elvis=elvis,addpedigree=True)
    
    # add datapath(s)
    
    if isinstance(datapath,list):
        longestdatapathname = max([len(item) for item in datapath])
        explog['datapath'] = np.zeros(len(explog),dtype='S%i' % longestdatapathname)
        explognumber=explog['explognumber']
        for idata in range(len(datapath)):
            explog['datapath'][explognumber == idata] = datapath[idata]
    

    #rootFile_name = explog['File_name'].copy()

    
    DataDict = {}
        
    # SELECTION OF OBSIDS
    
    selbool = (['FOCUS00' in item for item in explog['TEST']]) & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) & \
        (explog['Wavelength'] == inwavelength) # TESTS

    
    # Assess structure
    
    isconsistent = pilib.check_test_structure(explog,selbool,structure)
    
    
    # Build DataDict - And Labelling
    
    DataDict = dict(meta = dict(inwavelength=inwavelength,structure=structure))
    
    for CCDindex in [1,2,3]:
        
        CCDkey = 'CCD%i' % CCDindex
        
        DataDict[CCDkey] = dict()
        
        CCDselbool = selbool & (explog['CCD'] == CCDkey)
        
        if len(np.where(CCDselbool)[0]) == 0:
            continue
        
        Nobs = len(np.where(CCDselbool)[0])
        
        for key in explog.colnames:
            DataDict[CCDkey][key] = explog[key][CCDselbool].data.copy()
        
        DataDict[CCDkey]['time'] = np.array(map(pilib.get_dtobj,DataDict[CCDkey]['DATE'])).copy()
        
        
        Mirr_pos = DataDict[CCDkey]['Mirr_pos'].copy()
        Exptime = DataDict[CCDkey]['Exptime'].copy()
        
        label = np.zeros(Nobs,dtype='40str')
        
        uMirr_pos = np.sort(np.unique(Mirr_pos))
        

        for ixMP,iMP in enumerate(uMirr_pos):
            
            ixselMP = Mirr_pos == iMP
            
            label[ixselMP & (Exptime > 0)] = 'focus_%i' % ixMP
            label[ixselMP & (Exptime ==0)] = 'BGD'
            
        
        DataDict[CCDkey]['label'] = label.copy()

        
        rootFile_name = DataDict[CCDkey]['File_name'].copy()
        
        File_name  = ['%s.fits' % item for item in rootFile_name]
        
        DataDict[CCDkey]['Files'] = np.array(File_name).copy()
        
   
    
    return DataDict, isconsistent




def prep_data(DataDict,report,inputs,log=None):
    """ """
    
def basic_analysis(DataDict,report,inputs,log=None):
    """ """
    
def meta_analysis(DataDict,report,inputs,log=None):
    """ """
    

def feeder(inputs):
    """ """
    
    subtasks = [('prep',prep_data),('basic',basic_analysis),
                ('meta',meta_analysis)]
    
    return inputs,subtasks




