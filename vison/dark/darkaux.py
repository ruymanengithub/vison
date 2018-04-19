#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary functions to BIAS01 and DARK01 test scripts.

Created on Tue Nov 14 13:25:01 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import copy
import numpy as np
from pdb import set_trace as stop
from collections import OrderedDict

from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp
# END IMPORT


def get_bgd_model(ccdobj,extension=-1):
    
    prescan = ccdobj.prescan
    
    bgd = ccdmod.CCD()
    bgd.add_extension(np.zeros_like(
        ccdobj.extensions[extension].data, dtype='float32'))

    for Quad in ccdobj.Quads:

        Qimgmodel = ccdobj.get_region2Dmodel(Quad, area='img', kind='poly2D',
                                             pdegree=5, doFilter=False, filtsize=1,
                                             vstart=0, vend=2066, canonical=True, extension=extension).imgmodel.copy()

        # pre/overscans will get nulled
        Qmodel = ccdobj.get_quad(
            Quad, canonical=True, extension=extension).copy()

        Qmodel[prescan:prescan+Qimgmodel.shape[0],
               0:Qimgmodel.shape[1]] = Qimgmodel.copy()

        bgd.set_quad(Qmodel, Quad, canonical=True, extension=-1)

    return bgd.extensions[-1].data.copy()
    

def get_DarkDefectsMask(ccdobj, threshold, extension=-1):
    """ """
    maskdata = ccdobj.extensions[extension].data
    mask = (maskdata > threshold).astype('int32')
    return mask


def get_DarkDefectsMask_CDP(self, darkccdobj, threshold, subbgd=True, bgdmodel=None, extension=-1):
    """ """
    
    if subbgd:
        if bgdmodel is None:            
            bgdmodel = get_bgd_model(darkccdobj,extension=extension)

        darkccdobj.sub_bias(bgdmodel, extension=extension)

    mask = get_DarkDefectsMask(darkccdobj, threshold, extension=extension)

    ID = self.ID
    BLOCKID = self.BLOCKID
    CHAMBER = self.CHAMBER

    data = OrderedDict(mask=copy.deepcopy(mask))
    meta = OrderedDict(threshold=threshold, subbgd=subbgd,
                       function=get_DarkDefectsMask.__name__)

    maskcdp = cdp.CDP(data=data, meta=meta, ID=ID,
                      BLOCKID=BLOCKID, CHAMBER=CHAMBER)

    return maskcdp
