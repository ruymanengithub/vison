#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: DARK01

"Dark Current" analysis script

Created on Tue Aug 29 17:21:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from vison.pipe import lib as pilib
from vison.point import  lib as polib
from vison.datamodel import scriptic as sc
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
# END IMPORT

isthere = os.path.exists




DARK01_commvalues = dict(program='CALCAMP',test='DARK01',
  flushes=7,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=1,sinvflushp=500,
  chinj=0,tpump=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='DARK')
  


def build_DARK01_scriptdict(N,exptime,diffvalues=dict(),elvis='6.0.0'):
    """Builds DARK01 script structure dictionary.
    
    :param N: integer, number of frames to acquire.
    :param exptime: integer, ms, exposure time.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    DARK01_sdict = dict(col1=dict(frames=N,exptime=exptime))
    
    Ncols = len(DARK01_sdict.keys())    
    DARK01_sdict['Ncols'] = Ncols
    
                
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(DARK01_commvalues)
    
    DARK01_sdict = sc.update_structdict(DARK01_sdict,commvalues,diffvalues)    
    
    return DARK01_sdict


def check_data(DataDict,report,inputs,log=None):
    """ 
    METACODE
    
    Checks quality of ingested data.
    

    
    """
    
    # check common HK values are within safe / nominal margins
    # check voltages in HK match commanded voltages, within margins
    
    # f.e.ObsID, f.e.CCD, f.e.Q.:
        # measure offsets/means in pre-, img-, over-
        # measure std in pre-, img-, over-
    # assess std in pre- is within allocated margins
    # assess offsets/means in pre-, img-, over- are equal, within allocated  margins
    # assess offsets/means are within allocated margins
    
    # plot offsets/means vs. time
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
        # produce mask of hot pixels
            # count hot pixels / columns
        # produce a 2D poly model of masked-image, save coefficients
        # produce average profile along rows
        # produce average profile along cols
        # save 2D model and profiles in a pick file for each OBSID-CCD
        # measure and save RON after subtracting large scale structure
    
    # plot average profiles f. each CCD and Q (color coded by time)
    
    
    return DataDict,report
    
def meta_analysis(DataDict,report,inputs,log=None):
    """ """
    
    # f. each CCD, f. e. Q:
        # stack all ObsIDs to produce Master Dark
        # produce mask of hot pixels / columns
        # count hot pixels / columns
        # measure average profile along rows
        # measure average profile along cols
        
    # plot average profiles of Master Bias f. each Q
    # show Master Dark (image), include in report
    # report stats of defects, include in report
    # save name of MasterDark to DataDict, report
    # save name of Defects in Darkness Mask to DD, report
    
    return DataDict,report

