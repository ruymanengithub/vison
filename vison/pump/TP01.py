#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (vertical)

Created on Tue Aug 29 17:37:00 2017

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


TP01_commvalues = dict(program='CALCAMP',test='TP01',
  IDL=13000,IDH=18000,IG1=3500,IG2=6000,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=0,chinj_rows_on=2066,
  chinj_rows_off=0,chinj_repeat=1,id_width=57.1428571428571,
  id_delay=742.857142857144,
  tpump=1,ser_shuffles=0,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,#toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=250.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  


def build_TP01_scriptdict(Nshuffles_V,TOI_TPv,id_delays,diffvalues=dict()):
    """MISSING: different starting points (not implemented in ELVIS yet)."""
    
    assert len(id_delays) == 2
    
    phstarts = [1,2,3,4]
              
    TP01_sdict = dict()
    
    
    TP01_commvalues['ver_shuffles'] = Nshuffles_V
    
    # First Injection Drain Delay
    
    TP01_sdict['col1'] = dict(frames=1,tpump=0,comments='BGD',
              id_delay=id_delays[0])
    
    colcounter = 1
    for i,TOI_TP in enumerate(TOI_TPv):
        
        for k,phstart in enumerate(phstarts):
            colkey = 'col%i' % colcounter
            TP01_sdict[colkey] = dict(frames=1,toi_tpump=TOI_TP,
                     id_delay=id_delays[0],vert_tpump_mode=phstart)
            
            colcounter += 1
    
    # Second Injection Drain Delay
    
    
    TP01_sdict['col%i' % colcounter] = dict(frames=1,tpump=0,comments='BGD',
              id_delay=id_delays[1])
    colcounter += 1

    for j,TOI_TP in enumerate(TOI_TPv):
        
        for k,phstart in enumerate(phstarts):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP01_sdict[colkey] = dict(frames=1,toi_tpump=TOI_TP,
                         id_delay=id_delays[1],vert_tpump_mode=phstart)    
            
            colcounter += 1
            
    
    Ncols = len(TP01_sdict.keys())    
    TP01_sdict['Ncols'] = Ncols
    
    TP01_sdict = sc.update_structdict(TP01_sdict,TP01_commvalues,diffvalues)
    
    return TP01_sdict
