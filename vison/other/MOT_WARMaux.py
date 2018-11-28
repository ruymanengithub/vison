#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to MOT_WARM

Created on Tue Nov 06 14:10:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import copy

from vison.flat import BF01aux
from vison.plot import figclasses
from vison.plot import trends
#from MOT_FF import extract_overscan_profiles
# END IMPORT



prof_HER_ser_dict = dict(
    figname='MOTWARM_profs_HER_ser.png',
    caption='MOT\_WARM: HER profiles, serial direction.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim=[-0.002,0.005],
              suptitle='MOT\_WARM: Serial HER.')
)



def get_basic_prof1Dver_dict(tag):
    basic_prof1Dver_dict = dict(
    figname='MOTWARM_%s_profs1D_ver_allOBSIDs.png' % tag,
    caption='MOT\_WARM-%s: Average profiles across rows.' % tag,
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              suptitle='MOT\_WARM-%s: Profiles across rows.' % tag)
    )
    return basic_prof1Dver_dict



RAMP_bits_histo_dict = dict(
    figname='',
    caption='MOT\_WARM-%s: Bits Histogram on RAMP image.',
    meta=dict()
    )

def get_MW_figs():
    MW_figs = OrderedDict()
    
    MW_figs['MOTWbasic_HER_serial'] = [figclasses.Fig_Beam2DPlot, prof_HER_ser_dict]
    
    for tag in ['RAMP','CHINJ','FLAT']:
        MW_figs['MOTWbasic_prof1D_ver_%s' % tag] = get_basic_prof1Dver_dict(tag)
    
    for jCCD in [1,2,3]:
        MW_figs['RAMP_bits_histo_CCD%i' % jCCD] = [figclasses.Fig_Husk, RAMP_bits_histo_dict]
    
    MW_figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return MW_figs
