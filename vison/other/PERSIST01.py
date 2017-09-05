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

# END IMPORT


PER01_commvalues = dict(program='CALCAMP',test='PERSIST01',
  IDL=13000,IDH=18000,IG1=5000,IG2=5000,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,#exptime=0.,
  shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=1,
  chinj=0,chinj_rows_on=20,
  chinj_rows_off=20,chinj_repeat=1,id_width=100,
  id_delay=100,tpump=0,ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  

def build_PER01_scriptdict(exptSATUR,exptLATEN,diffvalues=dict()):
    """ """
    
    PER01_sdict = dict(col1=dict(frames=5,exptime=0,comments='REFER.'),
                       col2=dict(frames=1,exptime=exptSATUR,comments='EXPOSE'),
                       col3=dict(frames=3,exptime=exptLATEN,comments='LATENT'))   
    
    Ncols = len(PER01_sdict.keys())    
    PER01_sdict['Ncols'] = Ncols
    
    PER01_sdict = sc.update_structdict(PER01_sdict,PER01_commvalues,diffvalues)
    
    return PER01_sdict
