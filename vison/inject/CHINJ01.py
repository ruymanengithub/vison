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
from copy import deepcopy
# END IMPORT

isthere = os.path.exists



CHINJ01_commvalues = dict(program='CALCAMP',test='CHINJ01',
  IG2_T=5500,IG2_B=5500,
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='Normal',
  flushes=7,exptime=0.,
  chinj=1,chinj_on=30,
  chinj_of=100,
  id_wid=60,chin_dly=1,
#  id_delay=100,
  operator='who',
  comments='')
  


def build_CHINJ01_scriptdict(IDL,IDH,IG1s,id_delays,toi_chinj,diffvalues=dict(),
                             elvis='6.3.0'):
    """
    Builds CHINJ01 script structure dictionary.
    
    :param IDL: int, [mV], value of IDL (Inject. Drain Low).
    :param IDH: int, [mV], Injection Drain High.
    :param IG1s: list of 2 ints, [mV], [min,max] values of IG1.
    :param id_delays: list of 2 ints, [mV], injection drain delays (2).
    :param toi_chinj: int, [us], TOI-charge injection.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    CCDs = [1,2,3]
    halves = ['T','B']
    
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
        CHINJ01_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[0],toi_chinj=toi_chinj)
        
        for CCD in CCDs:
            for half in halves:
                CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD,half)] = IG1
        
        
        colcounter += 1
    
    # Second Injection Drain Delay
    
    colstart = colcounter

    for j,IG1 in enumerate(IG1v):
        colkey = 'col%i' % (colstart+j,)
        #print colkey
        CHINJ01_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[1],toi_chinj=toi_chinj)
        
        for CCD in CCDs:
            for half in halves:
                CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD,half)] = IG1

    
    Ncols = len(CHINJ01_sdict.keys())    
    CHINJ01_sdict['Ncols'] = Ncols
    
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(CHINJ01_commvalues)
    
    CHINJ01_sdict = sc.update_structdict(CHINJ01_sdict,commvalues,diffvalues)
    
    return CHINJ01_sdict

