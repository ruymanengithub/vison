#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FLAT0X

Flat-fields acquisition / analysis script

Created on Tue Aug 29 17:32:52 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
#from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
import datetime
from copy import deepcopy

# END IMPORT

isthere = os.path.exists


FLAT0X_commvalues = dict(program='CALCAMP',
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=0,
  tpump=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='')
  


def build_FLAT0X_scriptdict(exptimes,wavelength=800,testkey='FLAT0X',
                            diffvalues=dict(),elvis='6.0.0'):
    """Builds FLAT0X script structure dictionary.
    
    :param exptimes: list of ints, exposure times.
    :param wavelength: int, wavelength.
    :param testkey: char, test identifier.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    FW_ID = ogse.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    
    FLAT0X_commvalues['wavelength'] = 'Filter %i' % FW_IDX
    FLAT0X_commvalues['test'] = testkey
    
    
    FLAT0X_sdict = dict(col1=dict(frames=80,exptime=exptimes[0],comments='25pc'),
                        col2=dict(frames=60,exptime=exptimes[1],comments='50pc'),
                        col3=dict(frames=30,exptime=exptimes[2],comments='75pc'))

    Ncols = len(FLAT0X_sdict.keys())    
    FLAT0X_sdict['Ncols'] = Ncols
                
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(FLAT0X_commvalues)
    
    FLAT0X_sdict = sc.update_structdict(FLAT0X_sdict,commvalues,diffvalues)
    
    return FLAT0X_sdict
