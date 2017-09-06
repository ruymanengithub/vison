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

from copy import deepcopy
# END IMPORT 

isthere = os.path.exists

IG1comm = 3500
IG2comm = 6000

TP01_commvalues = dict(program='CALCAMP',test='TP01',
  IDL=13000,IDH=18000,
  IG1_1_T=IG1comm,IG1_2_T=IG1comm,IG1_3_T=IG1comm,
  IG1_1_B=IG1comm,IG1_2_B=IG1comm,IG1_3_B=IG1comm,  
  IG2_T=IG2comm,IG2_B=IG2comm,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=0,chinj_ser_wait=0,
  tpump=1,ser_shuffles=0,
  ver_shuffles=1,dwell_v=0,dwell_h=0,
  motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143,toi_tpump=1000,
  toi_rdout=1000,toi_chinj=250,
  comments='')
  


def build_TP01_scriptdict(Nshuffles_V,TOI_TPv,id_delays,
                    diffvalues=dict(),elvis='6.0.0'):
    """ """
    
    assert len(id_delays) == 2
    
    vpumpmodes = ['V1&2','V3&4','V4&1']
              
    TP01_sdict = dict()
    
    
    TP01_commvalues['ver_shuffles'] = Nshuffles_V
    
    # First Injection Drain Delay
    
    TP01_sdict['col1'] = dict(frames=1,tpump=0,comments='BGD',
              id_delay=id_delays[0])
    
    colcounter = 1
    for i,TOI_TP in enumerate(TOI_TPv):
        
        for k,vpumpmode in enumerate(vpumpmodes):
            colkey = 'col%i' % colcounter
            TP01_sdict[colkey] = dict(frames=1,toi_tpump=TOI_TP,
                     id_delay=id_delays[0],tpump_mode=vpumpmode)
            
            colcounter += 1
    
    # Second Injection Drain Delay
    
    
    TP01_sdict['col%i' % colcounter] = dict(frames=1,tpump=0,comments='BGD',
              id_delay=id_delays[1])
    colcounter += 1

    for j,TOI_TP in enumerate(TOI_TPv):
        
        for k,vpumpmode in enumerate(vpumpmode):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP01_sdict[colkey] = dict(frames=1,toi_tpump=TOI_TP,
                         id_delay=id_delays[1],tpump_mode=vpumpmode)    
            
            colcounter += 1
    
    
    Ncols = len(TP01_sdict.keys())    
    TP01_sdict['Ncols'] = Ncols

    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(TP01_commvalues)    
   
    TP01_sdict = sc.update_structdict(TP01_sdict,commvalues,diffvalues)
    
    return TP01_sdict
