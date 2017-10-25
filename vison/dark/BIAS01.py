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
import copy
from vison.datamodel import core


#from vison.support.report import Report


# END IMPORT

isthere = os.path.exists


HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']


BIAS01_commvalues = dict(program='CALCAMP',test='BIAS01',
  flushes=7,exptime=0.,shuttr=1,
  e_shuttr=0,vstart=1,vend=2086,
  siflush=0,#sinvflushp=500,
  chinj=0,
  s_tpump=0,
  v_tpump=0,
  motr_on=0,
  toi_fl=143.,toi_tp=1000.,toi_ro=1000.,toi_ch=1000.,
  wave=4,
  comments='BIAS')
  


def build_BIAS01_scriptdict(N,diffvalues=dict(),elvis='6.3.0'):
    """Builds BIAS01 script structure dictionary.
    
    :param N: integer, number of frames to acquire.
    :param diffvalues: dict, opt, differential values.
    :param elvis: char, ELVIS version.
    
    """
    
    BIAS01_sdict = dict(col1=dict(frames=N,exptime=0))

    Ncols = len(BIAS01_sdict.keys())    
    BIAS01_sdict['Ncols'] = Ncols

                
    commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(BIAS01_commvalues)
    
    BIAS01_sdict = sc.update_structdict(BIAS01_sdict,commvalues,diffvalues)

    return BIAS01_sdict



def filterexposures(structure,explogf,datapath,OBSID_lims,elvis):
    """ """
    
    wavedkeys = []


    # The filtering takes into account an expected structure for the 
    # acquisition script.

    
    # load exposure log(s)
    
    explog = pilib.loadexplogs(explogf,elvis=elvis,addpedigree=True,datapath=datapath)
    
    # SELECTION OF OBSIDS
    
    selbool = (explog['test']=='BIAS01') & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) 
    
    explog = explog[selbool]

    
    # Assess structure    
    checkreport = pilib.check_test_structure(explog,structure,CCDs=[1,2,3],
                                           wavedkeys=wavedkeys)
    
    # Labeling of exposures [optional]
    explog['label'] = np.array(['bias']*len(explog))
    
    
    return explog, checkreport


def check_data(dd,report,inputs,log=None):
    """ 
    
    BIAS01: Checks quality of ingested data.
    
    **METACODE**
    
    
    **TODO**: consider to raise an exception that
          would halt execution of task if 
          processing data could be just a waste of time.
          
    **TODO**: consider add a common binary "flags" variable as 
          input/output. It could go in DataDict, and reported in 
          log and report.
    
    ::
      
      check common HK values are within safe / nominal margins
      check voltages in HK match commanded voltages, within margins
    
      f.e.ObsID:
          f.e.CCD: 
              f.e.Q.:
                  measure offsets in pre-, img-, over-
                  measure std in pre-, img-, over-
      assess std in pre- is within allocated margins
      assess offsets in pre-, img-, over- are equal, within allocated  margins
      assess offsets are within allocated margins
    
      plot offsets vs. time
      plot std vs. time
    
      issue any warnings to log
      issue update to report
      update flags as needed
    
    """
    
    if log is not None:
        log.info('BIAS01.check_data')
    
    report.add_Section(Title='CHECK\_DATA',level=0)
    
    bypass = True # TESTS
    
    stop()
    
    # CHECK AND CROSS-CHECK HK: PENDING
    pass

    # Initialize new columns

    Xindices = copy.deepcopy(dd.indices)
    
    if 'Quad' not in Xindices.names:
        Xindices.append(core.vIndex('Quad',vals=pilib.Quads))
    
    dummyOffset = 2000.
    dummystd = 1.3
    
    newcolnames_off = ['offset_pre','offset_img','offset_ove']
    for newcolname_off in newcolnames_off:
        if newcolname_off in dd.colnames:
            dd.dropColumn(newcolname_off)
        dd.initColumn(newcolname_off,Xindices,dtype='float32',valini=dummyOffset)
    
    newcolnames_std = ['std_pre','std_img','std_ove']
    for newcolname_std in newcolnames_std:
        if newcolname_std in dd.colnames:
            dd.dropColumn(newcolname_std)
        dd.initColumn(newcolname_std,Xindices,dtype='float32',valini=dummystd)
    
    
    nObs,nCCD,nQuad = Xindices.shape
    Quads = Xindices[2].vals
    
    # Get statistics in different regions
    
    if not bypass:
    
        for iObs in range(nObs):
            
            for jCCD in range(nCCD):
                
                dpath = dd.mx['datapath'][iObs,jCCD]
                ffits = os.path.join(dpath,dd.mx['Files'][iObs,jCCD])
                
                ccdobj = ccd.CCD(ffits)
                
                for kQ in range(nQuad):
                    Quad = Quads[kQ]
                    
                    for reg in ['pre','img', 'ove']:
                        stats = ccdobj.get_stats(Quad,sector=reg,statkeys=['mean','std'],trimscan=[5,5],
                                ignore_pover=True,extension=-1)
                        dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = stats[0]
                        dd.mx['std_%s' % reg][iObs,jCCD,kQ] = stats[1]
                
    # Assess metrics are within allocated boundaries

    offset_lims = [1000.,3500.]
    
    std_lims = [0.5,2.]
    
    # absolute value of offsets
    
    compliance_offsets = dict()
    for reg in ['pre','img','ove']:
        compliance_offsets[reg] = ((dd.mx['offset_%s' % reg][:] >= offset_lims[0]) &\
                      (dd.mx['offset_%s' % reg][:] <= offset_lims[1]))
    
    # cross-check of offsets
    
    
    
    
    # absolute value of std
        
    compliance_std = dict()
    for reg in ['pre','img','ove']:
        compliance_std[reg] = ((dd.mx['std_%s' % reg][:] >= std_lims[0]) &\
                      (dd.mx['std_%s' % reg][:] <= std_lims[1]))
        
    # cross-check of stds
    
    
    
    # Do some Plots
    
    # offsets vs. time
    # std vs. time
    
    # Update Report, raise flags, fill-in log
    
    
    
    return dd, report
    

def prep_data(dd,report,inputs,log=None):
    """
    
    BIAS01: Preparation of data for further analysis.
    
    **METACODE**
    
    ::
        f.e. ObsID:
            f.e.CCD:
                f.e.Q:
                    subtract offset: save to FITS, UPDATE filename
    
    """
    
    if log is not None:
        log.info('BIAS01.prep_data')
    
    
    return dd,report
    
def basic_analysis(dd,report,inputs,log=None):
    """ 
    
    BIAS01: Basic analysis of data.
    
    **METACODE**
    
    :: 

        f. e. ObsID:
           f.e.CCD:
               f.e.Q:
                   produce a 2D poly model of bias, save coefficients
                   produce average profile along rows
                   produce average profile along cols
                   save 2D model and profiles in a pick file for each OBSID-CCD
                   measure and save RON after subtracting large scale structure
        plot RON vs. time f. each CCD and Q
        plot average profiles f. each CCD and Q (color coded by time)
 
    
    """
    
    if log is not None:
        log.info('BIAS01.basic_analysis')
   
    
    return dd,report
    
def meta_analysis(dd,report,inputs,log=None):
    """
    
    **METACODE**
    
    ::
    
        f. each CCD:
           f. e. Q:
               stack all ObsIDs to produce Master Bias
               measure average profile along rows
               measure average profile along cols
        plot average profiles of Master Bias f. each Q
        produce table with summary of results, include in report
        show Master Bias (image), include in report
        save name of MasterBias to DataDict, report
        
    """
    
    if log is not None:
        log.info('BIAS01.meta_analysis')
    
    
    return dd,report


def feeder(inputs,elvis='6.3.0'):
    """ """
    
    subtasks = [('check',check_data),('prep',prep_data),
                ('basic',basic_analysis),
                ('meta',meta_analysis)]
    
    N = inputs['N']
    if 'elvis' in inputs:
        elvis = inputs['elvis']
    if 'diffvalues' in inputs:
        diffvalues = inputs['diffvalues']
    else:
        diffvalues = {}
    
    
    scriptdict = build_BIAS01_scriptdict(N,diffvalues,elvis=elvis)
    
    inputs['structure'] = scriptdict
    inputs['subtasks'] = subtasks
    
    
    return inputs




