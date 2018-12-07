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

def fitPTC(means, var, debug=False):
    """Fits Photon Transfer Curve to obtain gain."""
    from sklearn.linear_model import RANSACRegressor
    from sklearn.metrics import mean_squared_error
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.pipeline import make_pipeline
    
    if debug:
        from pylab import plot,show
    
    sigmathresh=5. # sigma clipping
    poldeg = 2
    flags = flmod.Flags(fitPTC_flags)

    order = np.argsort(means)
    means = np.array(means)[order]
    var = np.array(var)[order]

    #ixmaxvar = np.argmax(var)
    #maxvar = var[ixmaxvar]
    #ixsel = np.where((means < means[ixmaxvar]) & (var < maxvar*0.95))
    
    try:
        
        pipe = make_pipeline(PolynomialFeatures(poldeg,interaction_only=True), 
               RANSACRegressor())
        pipe.fit(np.expand_dims(means,1),np.expand_dims(var,1))
        robpredict = np.squeeze(pipe.predict(np.expand_dims(means,1)))
        
        ixsel = np.where(np.abs(var-robpredict)/np.sqrt(robpredict)<sigmathresh)
        
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
    
    gain = 1./p[1]
    
    if ep[1]/gain>1.e-3:
        flags.add('POORFIT')
    
    quality = flags.value
    fitresults = dict(fit=p, efit=ep, gain=gain,
                      quadterm=p[0], rn=p[2], 
                      quality=quality)
    
    if debug:
        plot(means,var,'r.')
        plot(means[ixsel],var[ixsel],'b.')
        plot(means,np.polyval(p,means),'k-')
        show()
        stop()
    
    return fitresults

#def fitPTC_old(means, var, debug=False):
#    """Fits Photon Transfer Curve to obtain gain."""
#    if debug:
#        from pylab import plot,show
#    
#    poldeg = 2
#    flags = flmod.Flags(fitPTC_flags)
#
#    order = np.argsort(means)
#    means = np.array(means)[order]
#    var = np.array(var)[order]
#
#    ixmaxvar = np.argmax(var)
#    maxvar = var[ixmaxvar]
#    ixsel = np.where((means < means[ixmaxvar]) & (var < maxvar*0.95))
#    
#    try:
#        
#        res = np.polyfit(means[ixsel], var[ixsel], poldeg, full=False, cov=True)
#        p = res[0]
#        V = res[1]
#    
#    except:
#        p = np.zeros(poldeg+1)
#        p[0] = 0.01
#        V = np.zeros((poldeg+1,poldeg+1),dtype='float32')
#        
#        flags.add('EXCEPTION')
#    
#    ep = np.sqrt(np.diag(V))
#    
#    if np.any((ep == 0.) | np.isinf(ep) | np.isnan(ep)):
#        flags.add('BADERRORS')
#    
#    gain = 1./p[1]
#    
#    if ep[1]/gain>1.e-3:
#        flags.add('POORFIT')
#    
#    quality = flags.value
#    fitresults = dict(fit=p, efit=ep, gain=gain,
#                      quadterm=p[0], rn=p[2], 
#                      quality=quality)
#    
#    if debug:
#        plot(means,var,'r.')
#        plot(means[ixsel],var[ixsel],'b.')
#        show()
#        stop()
#    
#    return fitresults


def foo_bloom(means, var):
    """DUMMY function (PLACE-HOLDER)
    (Will) Finds blooming limit (where variance drops, if it does...).

    """
    res = dict(bloom=np.nan)
    return res
