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

