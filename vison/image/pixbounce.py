#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Pixel Bounce Analysis methods.

Created on Fri Mar  9 09:50:16 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from collections import OrderedDict
# END IMPORT


def get_pixbounce_from_overscan(ccdobj, thresholds=None):
    """
    Retrieves Hard Edge Respose for all Quadrants of a CCD.
    Uses the transition from image to overscan (along rows). Averages
    across rows. Input image should have high image-area fluence but not
    saturating. Rows can be filtered by average fluence in them via "thresholds"
    keyword. Do not use on images acquired with irradiated CCDs.

    """

    if thresholds is None:
        thresholds = [0, 2.**16]

    ixjump = 10  # index of hi-to-low transition in output profile
    Npix_prof = 25  # number of columns / pixels in output profile

    #x = np.arange(Nprof) + NAXIS1_Q-20-ixjump + 1

    profiles = OrderedDict()

    profiles['Npix_prof'] = Npix_prof
    profiles['ixjump'] = ixjump

    for Q in ccdobj.Quads:

        profiles[Q] = get_pixbounce_Quad(
            ccdobj, Q, Npix_prof, ixjump, thresholds)


def get_pixbounce_Quad(ccdobj, Q, Npix_prof, ixjump, thresholds):

    NAXIS1_Q = ccdobj.NAXIS1 // 2

    if ccdobj.voverscan > 0:
        endcol = -ccdobj.voverscan
    else:
        endcol = None

    xcol = np.arange(Npix_prof) + NAXIS1_Q - ccdobj.overscan - \
        ixjump  # col. count starts from 1

    profiles = OrderedDict()

    qdata = ccdobj.get_quad(Q, canonical=True)

    strip = qdata[-ccdobj.overscan - ixjump - 1:-ccdobj.overscan -
                  ixjump - 1 + Npix_prof, :endcol].copy()
    injection = np.mean(strip[0:ixjump, :endcol], axis=0)
    bias = np.mean(strip[ixjump + 3:, :endcol], axis=0)

    ixgood = np.where((injection <= thresholds[1]) & (
        injection >= thresholds[0]))
    Nrows = len(ixgood[0])

    strip = strip[:, ixgood[0]].copy()

    nstrip = (strip - bias[ixgood]) / (injection[ixgood] - bias[ixgood])

    avBIAS = bias[ixgood].mean()

    yprofile = np.mean(nstrip, axis=1)

    profiles['profile'] = (xcol.copy(), yprofile.copy())
    profiles['Nrows'] = Nrows
    profiles['avBIAS'] = avBIAS
    profiles['FPR'] = yprofile[ixjump + 1]

    return profiles
