#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ01

Charge injection calibration (part 1)
    Injection vs. IG1-IG2

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



CHINJ01_commvalues = dict(program='CALCAMP',test='CHINJ01',
#  IDL=13000,IDH=18000,IG1=5000,
  IDH=18000,IG2=5500,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=1,chinj_rows_on=30,
  chinj_rows_off=100,chinj_repeat=15,id_width=60,
#  id_delay=100,
  tpump=0,ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  


def build_CHINJ01_scriptdict(IDL,IDH,IG1s,id_delays,toi_chinj,diffvalues=dict()):
    """
    Builds CHINJ01 script structure dictionary.
    
    :param IDL: int, [mV], value of IDL (Inject. Drain Low).
    :param IDH: int, [mV], Injection Drain High.
    :param IG1s: list of 2 ints, [mV], [min,max] values of IG1.
    :param id_delays: list of 2 ints, [mV], injection drain delays (2).
    :param toi_chinj: int, [us], TOI-charge injection.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    assert len(IG1s) == 2
    assert len(id_delays) == 2
    
    
    dIG1 = 0.25 * 1.E3 # Vx1E3
    NIG1 = (IG1s[1]-IG1s[0])/dIG1+1
    IG1v = np.arange(NIG1)*dIG1+IG1s[0]
    
    CHINJ01_sdict = dict()
    
    # First Injection Drain Delay
    
    colcounter = 1
    for i,IG1 in enumerate(IG1v):
        colkey = 'col%i' % (i+1,)
        #print colkey
        CHINJ01_sdict[colkey] = dict(frames=1,IG1=IG1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[0],toi_chinj=toi_chinj)
        colcounter += 1
    
    # Second Injection Drain Delay
    
    colstart = colcounter

    for j,IG1 in enumerate(IG1v):
        colkey = 'col%i' % (colstart+j,)
        #print colkey
        CHINJ01_sdict[colkey] = dict(frames=1,IG1=IG1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[1],toi_chinj=toi_chinj)    

    
    Ncols = len(CHINJ01_sdict.keys())    
    CHINJ01_sdict['Ncols'] = Ncols
    
    CHINJ01_sdict = sc.update_structdict(CHINJ01_sdict,CHINJ01_commvalues,diffvalues)
    
    return CHINJ01_sdict

