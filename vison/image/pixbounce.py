#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Euclid-VIS Ground Calibration Campaign

Pixel Bounce Analysis methods.

Created on Fri Mar  9 09:50:16 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from collections import OrderedDict
# END IMPORT


def get_pixbounce_from_overscan(ccdobj, thresholds=[0, 2.**16]):
    """ 
    Retrieves Hard Edge Respose for all Quadrants of a CCD. 
    Uses the transition from image to overscan (along rows). Averages
    across rows. Input image should have high image-area fluence but not 
    saturating. Rows can be filtered by average fluence in them via "thresholds"
    keyword. Do not use on images acquired with irradiated CCDs.    

    """

    ixjump = 10  # index of hi-to-low transition in output profile
    Npix_prof = 25  # number of columns / pixels in output profile

    NAXIS1_Q = ccdobj.NAXIS1/2
    #x = np.arange(Nprof) + NAXIS1_Q-20-ixjump + 1

    profiles = OrderedDict()

    profiles['Npix_prof'] = Npix_prof
    profiles['ixjump'] = ixjump

    Quads = ccdobj.Quads

    oscan = ccdobj.overscan
    voscan = ccdobj.voverscan

    if voscan > 0:
        endcol = -voscan
    else:
        endcol = None

    xcol = np.arange(Npix_prof) + NAXIS1_Q-oscan - \
        ixjump  # col. count starts from 1

    for Q in Quads:

        profiles[Q] = OrderedDict()

        qdata = ccdobj.get_quad(Q, canonical=True)

        strip = qdata[-oscan-ixjump-1:-oscan -
                      ixjump-1+Npix_prof, :endcol].copy()
        injection = np.mean(strip[0:ixjump, :endcol], axis=0)
        bias = np.mean(strip[ixjump+3:, :endcol], axis=0)

        ixgood = np.where((injection <= thresholds[1]) & (
            injection >= thresholds[0]))
        Nrows = len(ixgood[0])
        #print '%i rows averaged' % (Nrows,)

        strip = strip[:, ixgood[0]].copy()

        nstrip = (strip - bias[ixgood]) / (injection[ixgood]-bias[ixgood])

        avBIAS = bias[ixgood].mean()
        #print 'Average Bias in SubChan %s = %.1f' % (Q,avBIAS)

        yprofile = np.mean(nstrip, axis=1)

        profiles[Q]['profile'] = (xcol.copy(), yprofile.copy())
        profiles[Q]['Nrows'] = Nrows
        profiles[Q]['avBIAS'] = avBIAS
        profiles[Q]['FPR'] = yprofile[ixjump+1]

    return profiles
