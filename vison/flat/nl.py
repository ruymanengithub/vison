#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in NL analysis.

Created on Mon Feb 5 15:51:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

#IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
#END IMPORT

def wrap_fitNL(raw_data,exptimes,col_labels,times=[],TrackFlux=True,subBgd=True):
    """ """
    # col1 == BGD
    # colEVEN = STAB
    # colODD = Fluences != 0
    
    col_numbers = np.array([int(item[3:]) for item in col_labels])
    
    ixboo_bgd = col_numbers == 1
    ixboo_stab = col_numbers % 2 != 0 # EVEN
    ixboo_fluences = col_labels %2 ==0 # ODD
    
    if subBgd:
        pass
    
    if TrackFlux:
        
    
    
    
    return None
    
    # COPIED BELOW HERE
    
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


