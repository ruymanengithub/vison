#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary functions to BIAS01 and DARK01 test scripts.

Created on Tue Nov 14 13:25:01 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import copy
import numpy as np
from pdb import set_trace as stop
from collections import OrderedDict

from vison.image import cosmetics
from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp
# END IMPORT



def get_DarkDefectsMask_CDP(taskobj, darkccdobj, threshold, subbgd=True, bgdmodel=None, extension=-1):
    """ """

    if subbgd:
        if bgdmodel is None:
            bgdmodel = cosmetics.get_bgd_model(darkccdobj, extension=extension)

        darkccdobj.sub_bias(bgdmodel, extension=extension)
    
    thresholds = [-1.E-6, threshold]
    darkdata = darkccdobj.extensions[extension].data.copy()
    mask = cosmetics.get_Thresholding_DefectsMask(darkdata, thresholds)

    ID = taskobj.ID
    BLOCKID = taskobj.BLOCKID
    CHAMBER = taskobj.CHAMBER

    data = OrderedDict(mask=copy.deepcopy(mask))
    meta = OrderedDict(threshold=threshold, subbgd=subbgd,
                       function=cosmetics.get_Thresholding_DefectsMask.__name__)

    maskcdp = cdp.CDP(data=data, meta=meta, ID=ID,
                      BLOCKID=BLOCKID, CHAMBER=CHAMBER)

    return maskcdp
