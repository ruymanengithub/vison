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

# END IMPORT

isthere = os.path.exists




DARK01_commvalues = dict(program='CALCAMP',test='DARK01',
  IDL=13000,IDH=18000,IG1=5000,IG2=5000,
  OD=26000,RD=16000,
  iphi1='TRUE',iphi2='TRUE',iphi3='TRUE',iphi4='FALSE',
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electro_shutter='FALSE',vstart=1,vend=2066,
  sinvflush='TRUE',chinj='FALSE',chinj_rows_on=20,
  chinj_rows_off=20,chinj_repeat=1,id_width=100,
  id_delay=100,tpump='FALSE',ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,Motor='FALSE',
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 1',pos_cal_mirror=0,
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='BIAS')
  


def build_DARK01_scriptdict(N,exptime,diffvalues=dict()):
    """ """
    
    DARK01_sdict = dict(col1=dict(N=N,exptime=exptime))
    
    Ncols = len(DARK01_sdict.keys())    
    DARK01_sdict['Ncols'] = Ncols
    
    DARK01_sdict = sc.update_structdict(DARK01_sdict,DARK01_commvalues,diffvalues)    
    
    return DARK01_sdict

