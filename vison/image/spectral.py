#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Spectral Analysis of VIS CCD Images.

Created on Tue Mar 19 10:53:30 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np

from scipy.fftpack import fft, fftfreq, fftshift
# END IMPORT


def serialize_ccdquadrant(ccdobj,Q,extension=-1,blankfiller='median'):
    """ """
    
    hdr = ccdobj.extensions[extension].header
    vstart = int(hdr['VSTART'])
    vend = int(hdr['VEND'])
    
    toi_p = 1000.E-6 #
    tpix = 4.75E-6 # s
    
    Nbuffer = int(toi_p * 4 / tpix)
    
    qdata = ccdobj.get_quad(Q,canonical=True,extension=extension)
    
    if blankfiller == 'median':
        filler = np.nanmedian(qdata[:,vstart:vend])
    elif blankfiller == 'mean':
        filler = np.nanmean(qdata[:,vstart:vend])
    elif blankfiller == 'zero':
        filler = 0.
    
    timeseries = np.zeros(shape=(1,),dtype='float32')
    
    Nrline = qdata.shape[0]
    
    for i in range(vstart,vend):
        line = np.zeros(Nrline+Nbuffer,dtype='float32')+filler
        line[0:Nrline] = qdata[:,i].copy()
        
        timeseries = np.concatenate((timeseries,line),axis=0)
    
    timeseries = timeseries[1:].copy()
    
    return timeseries

Tdef = 4.75E-6*3

def get_spectrum(timeseries,T=Tdef,Full=True):
    """ """
    N = len(timeseries)
    
    yspec = fft(timeseries)
    xspec = fftfreq(N,T)
    
    yspec = np.abs(yspec[0:N//2])
    xspec = xspec[0:N//2]
    
    #xspec = fftshift(xspec)
    #yspec = fftshift(yspec)    
        
    if Full:
        return xspec, yspec
    else:
        return yspec
    