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
    
    wavedkeys = []


    # The filtering takes into account an expected structure for the 
    # acquisition script.

    
    # load exposure log(s)
    
    explog = pilib.loadexplogs(explogf,elvis=elvis,addpedigree=True,datapath=datapath)
    
        
    # SELECTION OF OBSIDS
    
    selbool = (explog['TEST']=='BIAS01') & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) 
    
    explog = explog[selbool]

    
    # Assess structure
    
    checkreport = pilib.check_test_structure(explog,structure,CCDs=[1,2,3],
                                           wavedkeys=wavedkeys)
    
    # Labeling of exposures [optional]
    
    explog['label'] = np.array(['bias']*len(explog))
    
    
    
    return explog, checkreport


def check_data(DataDict,report,inputs,log=None):
    """ 
    METACODE
    
    Checks quality of ingested data.
    
    
    TODO: consider to raise an exception that
          would halt execution of task if 
          processing data could be just a waste of time.
          
    TODO: consider add a common binary "flags" variable as 
          input/output. It could go in DataDict, and reported in 
          log and report.
    
    """
    
    # check common HK values are within safe / nominal margins
    # check voltages in HK match commanded voltages, within margins
    
    # f.e.ObsID, f.e.CCD, f.e.Q.:
        # measure offsets in pre-, img-, over-
        # measure std in pre-, img-, over-
    # assess std in pre- is within allocated margins
    # assess offsets in pre-, img-, over- are equal, within allocated  margins
    # assess offsets are within allocated margins
    
    # plot offsets vs. time
    # plot std vs. time
    
    # issue any warnings to log
    # issue update to report
    
    return DataDict, report
    

def prep_data(DataDict,report,inputs,log=None):
    """
    METACODE
    
    Preparation of data for further analysis.
    
    """
    
    # f.e. ObsID, f.e.CCD, f.e.Q:
        # subtract offset: save to FITS, update filename
    
    
    return DataDict,report
    
def basic_analysis(DataDict,report,inputs,log=None):
    """ 
    METACODE
    
    Basic analysis of data.
    
    """
    
    # f. e. ObsID, f.e.CCD, f.e.Q:
        # produce a 2D poly model of bias, save coefficients
        # produce average profile along rows
        # produce average profile along cols
        # save 2D model and profiles in a pick file for each OBSID-CCD
        # measure and save RON after subtracting large scale structure
    # plot RON vs. time f. each CCD and Q
    # plot average profiles f. each CCD and Q (color coded by time)
    
    
    return DataDict,report
    
def meta_analysis(DataDict,report,inputs,log=None):
    """ """
    
    # f. each CCD, f. e. Q:
        # stack all ObsIDs to produce Master Bias
        # measure average profile along rows
        # measure average profile along cols
    # plot average profiles of Master Bias f. each Q
    # produce table with summary of results, include in report
    # show Master Bias (image), include in report
    # save name of MasterBias to DataDict, report
    
    return DataDict,report


def feeder(inputs,elvis='6.1.0'):
    """ """
    
    subtasks = [('prep',prep_data),('basic',basic_analysis),
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




