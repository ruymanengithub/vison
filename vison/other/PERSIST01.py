#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PERSIST01

CCD Persistence test

Created on Tue Aug 29 17:39:00 2017

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

# END IMPORT


PER01_commvalues = dict(program='CALCAMP',test='PERSIST01',
  rdmode='fwd_bas',
  flushes=7,#exptime=0.,
  shuttr=1,
  siflsh=1,siflsh_p=500,
  chinj=0,
  wave=4,mirr_pos=polib.mirror_nom['F4'],
  mirr_on=1,
  comments='')
  

def build_PER01_scriptdict(exptSATUR,exptLATEN,
                diffvalues=dict(),elvis='6.0.0'):
    """ 
    Builds PERSISTENCE01 script structure dictionary.
    
    :param exptSATUR: int, saturation exposure time.
    :param exptLATEN: int, latency exposure time.
    :param diffvalues: dict, opt, differential values.
    
    
    """
    
    PER01_sdict = dict(col1=dict(frames=5,exptime=0,comments='REFER.'),
                       col2=dict(frames=1,exptime=exptSATUR,comments='EXPOSE'),
                       col3=dict(frames=3,exptime=exptLATEN,comments='LATENT'))   
    
    Ncols = len(PER01_sdict.keys())    
    PER01_sdict['Ncols'] = Ncols
    
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(PER01_commvalues)               
               
    PER01_sdict = sc.update_structdict(PER01_sdict,commvalues,diffvalues)
    
    return PER01_sdict


def check_data(DataDict,report,inputs,log=None):
    """ 

    PERSIST01: Checks quality of ingested data.
    

    **METACODE**
    
    ::

        check common HK values are within safe / nominal margins
        check voltages in HK match commanded voltages, within margins
    
        f.e.ObsID:
            f.e.CCD:
                f.e.Q.:
                    measure offsets in pre-, over-
                    measure std in pre-, over-
                    measure fluence in apertures around Point Sources
        
        assess std in pre- (~RON) is within allocated margins
        assess offsets in pre-, and over- are equal, within allocated  margins
        assess fluence is ~expected within apertures (PS) for each frame (pre-satur, satur, post-satur)
        
    
        plot point source fluence vs. OBSID, all sources
        [plot std vs. time]
    
        issue any warnings to log
        issue update to report          

    
    """

def prep_data(DataDict,report,inputs,log=None):
    """
    
    **METACODE**
    
    ::
    
        Preparation of data for further analysis:

        f.e. ObsID [images with TPing only]:
            f.e.CCD:
                f.e.Q:
                    subtract offset

    """
    
    return DataDict,report

def basic_analysis():
    """
    
    Basic analysis of data.

    **METACODE**
    
    ::

        f.e.CCD:
            f.e.Q:
                use SATURATED frame to generate pixel saturation MASK
                measure stats in pix satur MASK across OBSIDs 
                 (pre-satur, satur, post-satur)

    
    """
    
def meta_analysis():
    """
    
    Meta-analysis of data.

    **METACODE**
    
    ::
        
        
        f.e.CCD:
            f.e.Q:
                estimate delta-charge_0 and decay tau from time-series

        report:  
            persistence level (delta-charge_0) and time constant
        

    """
    
