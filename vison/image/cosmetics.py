#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:55:12 2018

@author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
from collections import OrderedDict

from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp
# END IMPORT


def get_bgd_model(ccdobj, extension=-1, margins=None):

    prescan = ccdobj.prescan

    bgd = ccdmod.CCD()
    bgd.add_extension(np.zeros_like(
        ccdobj.extensions[extension].data, dtype='float32'))

    for Quad in ccdobj.Quads:

        Qimgmodel = ccdobj.get_region2Dmodel(Quad, area='img', kind='poly2D',
                                             pdegree=5, doFilter=False,
                                             doBin=False, filtsize=30,
                                             vstart=0, vend=2066, canonical=True,
                                             extension=extension).imgmodel.copy()

        # pre/overscans will not be modified unles margins is not None

        Qmodel = ccdobj.get_quad(
            Quad, canonical=True, extension=extension).copy()
        if margins is not None:
            Qmodel[:] = margins

        Qmodel[prescan:prescan + Qimgmodel.shape[0],
               0:Qimgmodel.shape[1]] = Qimgmodel.copy()

        bgd.set_quad(Qmodel, Quad, canonical=True, extension=-1)

    return bgd.extensions[-1].data.copy()


def get_Thresholding_DefectsMask(maskdata, thresholds):
    """ """
    #maskdata = ccdobj.extensions[extension].data
    mask = ((maskdata < thresholds[0]) | (maskdata > thresholds[1])).astype('int32')
    return mask
