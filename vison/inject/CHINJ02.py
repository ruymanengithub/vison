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


IG1comm = 6.
IG2comm = 4.

CHINJ02_commvalues = dict(program='CALCAMP',test='CHINJ02',
  IG1_1_T=IG1comm,IG1_2_T=IG1comm,IG1_3_T=IG1comm,
  IG1_1_B=IG1comm,IG1_2_B=IG1comm,IG1_3_B=IG1comm,
  IG2_T=IG2comm,IG2_B=IG2comm,
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,exptime=0.,
  siflsh=0,
  chinj=1,chinj_on=30,
  chinj_of=100,
  id_wid=60,
  chin_dly=1,
  comments='')
  



def build_CHINJ02_scriptdict(IDLs,IDH,id_delays,toi_chinj,diffvalues=dict(),
                             elvis='6.3.0'):
    """ 
    Builds CHINJ02 script structure dictionary.
    
    :param IDLs: list of 2 ints, [mV], [min,max] values of IDL (Inject. Drain Low).
    :param IDH: int, [mV], Injection Drain High.
    :param id_delays: list of 2 ints, [mV], injection drain delays (2).
    :param toi_chinj: int, [us], TOI-charge injection.
    :param diffvalues: dict, opt, differential values.
    
    """
        
    assert len(IDLs) == 2
    assert len(id_delays) == 2
    
    dIDL = 0.25 # V
    NIDL = (IDLs[1]-IDLs[0])/dIDL+1
    IDLv = np.arange(NIDL)*dIDL+IDLs[0]
    
    CHINJ02_sdict = dict()
    
    
    # First Injection Drain Delay
    
    colcounter = 1
    for i,IDL in enumerate(IDLv):
        colkey = 'col%i' % (i+1,)
        CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_dly=id_delays[0],toi_ch=toi_chinj)
        colcounter += 1
    
    # Second Injection Drain Delay
    
    colstart  = colcounter

    for j,IDL in enumerate(IDLv):
        colkey = 'col%i' % (colstart+j,)
        CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_dly=id_delays[1],toi_ch=toi_chinj)
    
    Ncols = len(CHINJ02_sdict.keys())    
    CHINJ02_sdict['Ncols'] = Ncols
    
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(CHINJ02_commvalues)
                 
    CHINJ02_sdict = sc.update_structdict(CHINJ02_sdict,commvalues,diffvalues)
    
    return CHINJ02_sdict
