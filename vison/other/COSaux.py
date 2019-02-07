#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to COSMETICS00

Created on Wed Dec 06 18:04:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import copy
from matplotlib import cm

from vison.datamodel import cdp
from vison.flat import BF01aux
from vison.plot import figclasses
from vison.plot import trends
#from MOT_FF import extract_overscan_profiles
# END IMPORT

def gt_meta_MASK_dict(test,masktype):
    ntest = st.replace(test,'_','\_')
    return dict(
            figname='%s_MASK_2Dimgshow_%s.png' % (test,masktype),
            caption='%s: Masks CCDs [%s]. Smoothed with gaussian kernel to highlight details.' % (ntest,masktype),
            meta=dict(doLegend=False,
              doColorbar=False,
              suptitle='%s/%s:Quadrant Images.' % (ntest,masktype),
              corekwargs=dict(cmap=cm.gray,aspect='auto',norm=None,
                                origin='lower left'))
                )



def get_COS_figs():
    """ """
    COS_figs = OrderedDict()    
    COS_figs['Masks_DARK'] = [figclasses.Fig_BeamImgShow, gt_meta_MASK_dict('COSMETICS00','DARK')]
    COS_figs['Masks_FLAT'] = [figclasses.Fig_BeamImgShow, gt_meta_MASK_dict('COSMETICS00','FLAT')]
    COS_figs['Masks_MERGE'] = [figclasses.Fig_BeamImgShow, gt_meta_MASK_dict('COSMETICS00','MERGE')]    
    COS_figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    
    return COS_figs


def get_CDP_lib():
       
    CDP_lib = OrderedDict()
    
    def_tb_cdp = cdp.Tables_CDP()
    def_tb_cdp.rootname = 'COSMETICS00_DEFECTS_TB'
    
    CDP_lib['DEF_TB'] = def_tb_cdp
    
    return CDP_lib
