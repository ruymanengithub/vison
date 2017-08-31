#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ02

Charge injection calibration (part 2)
    Injection vs. IDL (injection threshold)

Created on Tue Aug 29 17:36:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#import datetime

# END IMPORT

isthere = os.path.exists



CHINJ02_commvalues = dict(program='CALCAMP',test='CHINJ02',
  #IDL=13000,IDH=18000,
  IG1=6000,IG2=4000,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1='TRUE',iphi2='TRUE',iphi3='TRUE',iphi4='FALSE',
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter='FALSE',vstart=1,vend=2066,
  sinvflush='FALSE',chinj='FALSE',chinj_rows_on=30,
  chinj_rows_off=100,chinj_repeat=15,id_width=57.1428571428571,
#  id_delay=100,
  tpump='FALSE',ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor='FALSE',
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  



def build_CHINJ02_scriptdict(IDLs,IDH,id_delays,diffvalues=dict()):
    """ """
    dIDL = 0.25*1.E3 # Vx1E3
    NIDL = (IDLs[1]-IDLs[0])/dIDL+1
    IDLv = np.arange(NIDL)*dIDL+IDLs[0]
    
    CHINJ02_sdict = dict()
    
    # First Injection Drain Delay
    
    colcounter = 1
    for i,IDL in enumerate(IDLv):
        colkey = 'col%i' % (i+1,)
        CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[0])
        colcounter += 1
    
    # Second Injection Drain Delay
    
    colstart  = colcounter

    for j,IDL in enumerate(IDLv):
        colkey = 'col%i' % (colstart+j,)
        CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[1])
    
    Ncols = len(CHINJ02_sdict.keys())    
    CHINJ02_sdict['Ncols'] = Ncols
    
    CHINJ02_sdict = sc.update_structdict(CHINJ02_sdict,CHINJ02_commvalues,diffvalues)
    
    return CHINJ02_sdict