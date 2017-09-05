#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (serial)

Created on Tue Aug 29 17:38:00 2017

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

# END IMPORT 

isthere = os.path.exists


TP02_commvalues = dict(program='CALCAMP',test='TP02',
  IDL=13000,IDH=18000,IG1=3500,IG2=6000,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=100,
  sinvflush=0,chinj=0,chinj_rows_on=2066,
  chinj_rows_off=0,chinj_repeat=1,id_width=60,
  id_delay=750,
  tpump=1,ser_shuffles=5000,
  ver_shuffles=0,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000,
  toi_rdout=1000,toi_chinj=250,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  


def build_TP02_scriptdict(Nshuffles_H,dwell_hv,id_delays,diffvalues=dict()):
    """MISSING: different starting points (not implemented in ELVIS yet)."""
    
    assert len(id_delays) == 2
    
    phstarts = [1,2,3]
              
    TP02_sdict = dict()
    
    
    TP02_commvalues['ser_shuffles'] = Nshuffles_H
    
    # First Injection Drain Delay
    
    TP02_sdict['col1'] = dict(frames=1,tpump=0,comments='BGD',
              id_delay=id_delays[0])
    
    colcounter = 1
    for i,dwell_h in enumerate(dwell_hv):
        
        for k,phstart in enumerate(phstarts):
            colkey = 'col%i' % colcounter
            TP02_sdict[colkey] = dict(frames=1,dwell_h=dwell_h,
                     id_delay=id_delays[0],ser_tpump_mode=phstart)
            
            colcounter += 1
    
    
    # Second Injection Drain Delay
    
    TP02_sdict['col%i' % colcounter] = dict(frames=1,tpump=0,comments='BGD',
              id_delay=id_delays[1])
    colcounter += 1 
    

    for j,dwell_h in enumerate(dwell_hv):
        
        for k,phstart in enumerate(phstarts):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP02_sdict[colkey] = dict(frames=1,dwell_h=dwell_h,
                         id_delay=id_delays[1],ser_tpump_mode=phstart)    
            
            colcounter += 1
            
    
    Ncols = len(TP02_sdict.keys())    
    TP02_sdict['Ncols'] = Ncols
    
    TP02_sdict = sc.update_structdict(TP02_sdict,TP02_commvalues,diffvalues)
    
    return TP02_sdict
