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
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,#exptime=0.,
  shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=1,sinvflushp=500,
  chinj=0,tpump=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
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
