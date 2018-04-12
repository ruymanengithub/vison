#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary functions to BIAS01 and DARK01 test scripts.

Created on Tue Nov 14 13:25:01 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
from collections import OrderedDict

from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp
# END IMPORT


def get_DarkDefectsMask(ccdobj, threshold, subbgd=True, extension=-1, bgdmodel=None,
                        Full=False):
    """ """

    Quads = ccdobj.Quads
    #shape = ccdobj.shape
    #NQx = shape[0]/2
    #NQy = shape[1]/2
    prescan = ccdobj.prescan

    if subbgd:

        if bgdmodel is None:

            bgd = ccdmod.CCD()
            bgd.add_extension(np.zeros_like(
                ccdobj.extensions[extension].data, dtype='float32'))

            for Quad in Quads:

                Qimgmodel = ccdobj.get_region2Dmodel(Quad, area='img', kind='poly2D',
                                                     pdegree=5, doFilter=False, filtsize=1,
                                                     vstart=0, vend=2066, canonical=True, extension=extension).imgmodel.copy()

                # pre/overscans will get nulled
                Qmodel = ccdobj.get_quad(
                    Quad, canonical=True, extension=extension).copy()

                Qmodel[prescan:prescan+Qimgmodel.shape[0],
                       0:Qimgmodel.shape[1]] = Qimgmodel.copy()

                bgd.set_quad(Qmodel, Quad, canonical=True, extension=-1)

            bgdmodel = bgd.extensions[-1].data.copy()

        ccdobj.sub_bias(bgdmodel, extension=extension)

    maskdata = ccdobj.extensions[extension].data

    mask = (maskdata > threshold).astype('int32')

    if Full:
        return mask, bgdmodel
    else:
        return mask


def get_DarkDefectsMask_CDP(self, darkccdobj, threshold, subbgd=True, extension=-1, bgdmodel=None):
    """ """

    mask = get_DarkDefectsMask(darkccdobj, threshold, subbgd=subbgd, extension=extension, bgdmodel=bgdmodel,
                               Full=False)

    ID = self.ID
    BLOCKID = self.BLOCKID
    CHAMBER = self.CHAMBER

    data = OrderedDict(mask=mask.copy())
    meta = OrderedDict(threshold=threshold, subbgd=subbgd,
                       function=get_DarkDefectsMask.__name__)

    maskcdp = cdp.CDP(data=data, meta=meta, ID=ID,
                      BLOCKID=BLOCKID, CHAMBER=CHAMBER)

    return maskcdp
