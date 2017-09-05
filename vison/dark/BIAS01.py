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

# END IMPORT

isthere = os.path.exists




BIAS01_commvalues = dict(program='CALCAMP',test='BIAS01',
  IDL=13000,IDH=18000,IG1=5000,IG2=5000,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=0,chinj_rows_on=20,
  chinj_rows_off=20,chinj_repeat=1,id_width=100,
  id_delay=100,tpump=0,ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='BIAS')
  


def build_BIAS01_scriptdict(N,diffvalues=dict()):
    """Builds BIAS01 script structure dictionary.
    
    :param N: integer, number of frames to acquire.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    BIAS01_sdict = dict(col1=dict(frames=N,exptime=0))

    Ncols = len(BIAS01_sdict.keys())    
    BIAS01_sdict['Ncols'] = Ncols
    
    BIAS01_sdict = sc.update_structdict(BIAS01_sdict,BIAS01_commvalues,diffvalues)
    
    return BIAS01_sdict
