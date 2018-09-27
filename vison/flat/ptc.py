#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in PTC analysis.

Created on Thu Sep 14 16:29:36 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from collections import OrderedDict

from vison.support import flags as flmod
# END IMPORT


fitPTC_flags = OrderedDict()
fitPTC_flags['EXCEPTION'] = [2**0L]
fitPTC_flags['POORFIT'] = [2**1L]
fitPTC_flags['BADERRORS'] = [2**2L]


def fitPTC(means, var):
    """Fits Photon Transfer Curve to obtain gain."""
    
    poldeg = 2
    flags = flmod.Flags(fitPTC_flags)

    order = np.argsort(means)
    means = np.array(means)[order]
    var = np.array(var)[order]

    ixmaxvar = np.argmax(var)
    maxvar = var[ixmaxvar]
    ixsel = np.where((means < means[ixmaxvar]) & (var < maxvar*0.95))
    
    try:
        
        res = np.polyfit(means[ixsel], var[ixsel], poldeg, full=False, cov=True)
        p = res[0]
        V = res[1]
    
    except:
        p = np.zeros(poldeg+1)
        p[0] = 0.01
        V = np.zeros((poldeg+1,poldeg+1),dtype='float32')
        
        flags.add('EXCEPTION')
    
    ep = np.sqrt(np.diag(V))
    
    if np.any((ep == 0.) | np.isinf(ep) | np.isnan(ep)):
        flags.add('BADERRORS')
    
    quality = flags.value
    fitresults = dict(fit=p, efit=ep, gain=1./p[1],
                      quadterm=p[0], rn=p[2], quality=quality)
    
    return fitresults


def foo_bloom(means, var):
    """DUMMY function (PLACE-HOLDER)
    (Will) Finds blooming limit (where variance drops, if it does...).

    """
    res = dict(bloom=np.nan)
    return res
