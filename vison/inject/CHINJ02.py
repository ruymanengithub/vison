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
from copy import deepcopy
# END IMPORT

isthere = os.path.exists


IG1comm = 6000
IG2comm = 4000

CHINJ02_commvalues = dict(program='CALCAMP',test='CHINJ02',
  IG1_1_T=IG1comm,IG1_2_T=IG1comm,IG1_3_T=IG1comm,
  IG1_1_B=IG1comm,IG1_2_B=IG1comm,IG1_3_B=IG1comm,
  IG2_T=IG2comm,IG2_B=IG2comm,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=1,chinj_rows_on=30,
  chinj_rows_off=100,id_width=60,
  tpump=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='')
  



def build_CHINJ02_scriptdict(IDLs,IDH,id_delays,toi_chinj,diffvalues=dict(),
                             elvis='6.0.0'):
    """ 
    Builds CHINJ02 script structure dictionary.
    
    :param IDLs: list of 2 ints, [mV], [min,max] values of IDL (Inject. Drain Low).
    :param IDH: int, [mV], Injection Drain High.
    :param id_delays: list of 2 ints, [mV], injection drain delays (2).
    :param toi_chinj: int, [us], TOI-charge injection.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    CCDs = [1,2,3]
    halves = ['T','B']
    
    assert len(IDLs) == 2
    assert len(id_delays) == 2
    
    dIDL = 0.25*1.E3 # Vx1E3
    NIDL = (IDLs[1]-IDLs[0])/dIDL+1
    IDLv = np.arange(NIDL)*dIDL+IDLs[0]
    
    CHINJ02_sdict = dict()
    
    
    # First Injection Drain Delay
    
    colcounter = 1
    for i,IDL in enumerate(IDLv):
        colkey = 'col%i' % (i+1,)
        CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[0],toi_chinj=toi_chinj)
        colcounter += 1
    
    # Second Injection Drain Delay
    
    colstart  = colcounter

    for j,IDL in enumerate(IDLv):
        colkey = 'col%i' % (colstart+j,)
        CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[1],toi_chinj=toi_chinj)
    
    Ncols = len(CHINJ02_sdict.keys())    
    CHINJ02_sdict['Ncols'] = Ncols
    
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(CHINJ02_commvalues)
                 
    CHINJ02_sdict = sc.update_structdict(CHINJ02_sdict,commvalues,diffvalues)
    
    return CHINJ02_sdict
