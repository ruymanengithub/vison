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

from copy import deepcopy
# END IMPORT 

isthere = os.path.exists

IG1comm = 3500
IG2comm = 6000

TP02_commvalues = dict(program='CALCAMP',test='TP02',
  IDL=13000,IDH=18000,
  IG1_1_T=IG1comm,IG1_2_T=IG1comm,IG1_3_T=IG1comm,
  IG1_1_B=IG1comm,IG1_2_B=IG1comm,IG1_3_B=IG1comm,  
  IG2_T=IG2comm,IG2_B=IG2comm,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=100,
  sinvflush=0,chinj=0,
  tpump=1,ser_shuffles=5000,
  ver_shuffles=0,dwell_v=0,dwell_h=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000,
  toi_rdout=1000,toi_chinj=250,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='')
  


def build_TP02_scriptdict(Nshuffles_H,dwell_hv,id_delays,
                diffvalues=dict(),elvis='6.0.0'):
    """MISSING: different starting points (not implemented in ELVIS yet)."""
    
    assert len(id_delays) == 2
    
    sermodes = ['S1&2','S3&4']
              
    TP02_sdict = dict()
    
    
    TP02_commvalues['ser_shuffles'] = Nshuffles_H
    
    # First Injection Drain Delay
    
    TP02_sdict['col1'] = dict(frames=1,tpump=0,ser_tpump=0,
              comments='BGD',id_delay=id_delays[0])
    
    colcounter = 1
    for i,dwell_h in enumerate(dwell_hv):
        
        for k,sermode in enumerate(sermodes):
            colkey = 'col%i' % colcounter
            TP02_sdict[colkey] = dict(frames=1,dwell_h=dwell_h,
                     id_delay=id_delays[0],ser_tpump_mode=sermode)
            
            colcounter += 1
    
    
    # Second Injection Drain Delay
    
    TP02_sdict['col%i' % colcounter] = dict(frames=1,tpump=0,ser_tpump=0,
              comments='BGD',id_delay=id_delays[1])
    colcounter += 1 
    

    for j,dwell_h in enumerate(dwell_hv):
        
        for k,sermode in enumerate(sermodes):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP02_sdict[colkey] = dict(frames=1,dwell_h=dwell_h,
                         id_delay=id_delays[1],ser_tpump_mode=sermode)    
            
            colcounter += 1
            
    
    Ncols = len(TP02_sdict.keys())    
    TP02_sdict['Ncols'] = Ncols
              
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(TP02_commvalues)
    
    TP02_sdict = sc.update_structdict(TP02_sdict,commvalues,diffvalues)
    
    return TP02_sdict
