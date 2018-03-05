#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in PTC analysis.

Created on Thu Sep 14 16:29:36 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

#IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
#END IMPORT

def fitPTC(means,var): 
    """ """
    
    order = np.argsort(means)
    means = np.array(means)[order]
    var = np.array(var)[order]   
    
    ixmaxvar = np.argmax(var)
    maxvar = var[ixmaxvar]
    
    ixsel = np.where((means < means[ixmaxvar]) & (var < maxvar*0.95))
    fvar = var[ixsel]
    fmeans = means[ixsel]
    
    
    res = np.polyfit(fmeans,fvar,2,full=False,cov=True)
    
    p = res[0]
    V = res[1]
    ep = np.sqrt(np.diag(V))
    
    g = 1./p[1]
    cuadterm = p[0]
    rn = p[2]
    
    # Bad results flagging MISSING!
    quality=0
    
    fitresults = dict(fit=p,efit=ep,gain=g,cuadterm=cuadterm,rn=rn,quality=quality)
    

    return fitresults


def foo_bloom(means,var):
    """DUMMY function (PLACE-HOLDER)"""
    res = dict(bloom = np.nan)
    return res
