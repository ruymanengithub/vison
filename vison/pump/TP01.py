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

IG1=6
IG2=5


TP01_commvalues = dict(program='CALCAMP',test='TP01',
  flushes=7,exptime=0.,shuttr=0,
  e_shuttr=0,vstart=1,vend=2066,
  siflsh=1,siflsh_p=500,
  IDL=11,IDH=18,
  IG1_1_T=IG1,IG1_2_T=IG1,IG1_3_T=IG1,
  IG1_1_B=IG1,IG1_2_B=IG1,IG1_3_B=IG1,
  IG2_T=IG2,IG2_B=IG2,
  chinj=1,chinj_on=2066,chinj_of=0,
  chin_dly=0,
  v_tpump=1,s_tpump=0,
  v_tp_cnt=1000,
  dwell_v=0,dwell_s=0,
  motr_on=0,
  comments='')
  


def build_TP01_scriptdict(Nshuffles_V,TOI_TPv,id_delays,
                    diffvalues=dict(),elvis='6.3.0'):
    """ """
    
    assert len(id_delays) == 2
    
    vpumpmodes = [123,234,341,412]
              
    TP01_sdict = dict()
    
    
    TP01_commvalues['ver_shuffles'] = Nshuffles_V
    
    # First Injection Drain Delay
    
    TP01_sdict['col1'] = dict(frames=1,v_tpump=0,comments='BGD',
              id_delay=id_delays[0])
    
    colcounter = 2
    for i,TOI_TP in enumerate(TOI_TPv):
        
        for k,vpumpmode in enumerate(vpumpmodes):
            colkey = 'col%i' % colcounter
            TP01_sdict[colkey] = dict(frames=1,toi_tp=TOI_TP,
                     id_dly=id_delays[0],v_tpmod=vpumpmode)
            
            colcounter += 1
    
    # Second Injection Drain Delay
    
    
    TP01_sdict['col%i' % colcounter] = dict(frames=1,v_tpump=0,comments='BGD',
              id_delay=id_delays[1])
    colcounter += 1

    for j,TOI_TP in enumerate(TOI_TPv):
        
        for k,vpumpmode in enumerate(vpumpmodes):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP01_sdict[colkey] = dict(frames=1,toi_tp=TOI_TP,
                         id_dly=id_delays[1],v_tpmod=vpumpmode)    
            
            colcounter += 1
    
    
    Ncols = len(TP01_sdict.keys())    
    TP01_sdict['Ncols'] = Ncols

    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(TP01_commvalues)    
   
    TP01_sdict = sc.update_structdict(TP01_sdict,commvalues,diffvalues)
    
    return TP01_sdict
