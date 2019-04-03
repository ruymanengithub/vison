#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module to provide common tools for analysis of Charge Injection acquisitions.

Created on Thu Sep 14 15:32:10 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import copy

from astropy.io import fits as fts
#from pylab import plot,imshow,show
from scipy.optimize import curve_fit
from scipy import special
from matplotlib import pyplot as plt
# END IMPORT

def_bgd_drop = [0., 15.]
lineoffsets = dict(E=0, F=0, G=0, H=0)

def extract_injection_lines(ccdobj, Q, pattern, VSTART=0,
                            VEND=2066, suboffmean=False,debug=False):
    """     
    ccdobj: ccd.CCD object
    pattern: non,noff,nrep (lines on, off, repeatitions)
    VSTART: VSTART
    VEND: VEND
    suboffmean: bool, subtract median of non-injected lines
    lineoffset: integer, to account for a shift between readout lines and charge
                   injection pattern.
    """
    
    lineoffset = lineoffsets[Q]

    npre = ccdobj.prescan
    npost = ccdobj.overscan
    
    quaddata = ccdobj.get_quad(Q, canonical=True, extension=-1)
    
    masked=False
    if isinstance(quaddata,np.ma.masked_array):
        masked=True

    non, noff, nrep = pattern

    npercycle = non+noff

    NX = quaddata.shape[0]

    nlines = VEND-VSTART

    rowmedians = np.zeros(nlines, dtype='float32')
    rowstds = np.zeros(nlines, dtype='float32')
    rowix = np.arange(VSTART, VEND)

    stack_2d = np.zeros((nrep, NX-(npre+npost), npercycle),
                        dtype='float32') + np.nan
    stack_2d = np.ma.masked_array(data=stack_2d,
                                  mask=np.zeros_like(stack_2d,dtype='bool'))
               

    for ii in (rowix):

        row = quaddata[npre:-1*npost, ii].copy()
        
        rowmedians[ii-VSTART] = np.nanmedian(row)
        rowstds[ii-VSTART] = np.nanstd(row)

        icycle = max((ii+lineoffset) / npercycle, 0)
        ix_in_cycle = (ii+lineoffset) - npercycle * icycle
                      
        if icycle >= nrep:
            break

        stack_2d[icycle, :, ix_in_cycle] = row.copy()
    
    
    stack_2d.mask[np.where(np.isnan(stack_2d.data))]=True

    if suboffmean:
        for icycle in range(nrep):
            stack_2d[icycle, :, :] -= np.nanmean(stack_2d[icycle, :, non:])
    
    
    stacked_2d = np.ma.mean(stack_2d, axis=0)
    
    #mask2d = np.sum(stack_2d.mask,axis=0) # TESTS    
    #fts.writeto('stacked_2d.fits',stacked_2d.data.transpose(),overwrite=True) # TESTS
    #fts.writeto('mask2d.fits',mask2d.transpose().astype('float32'),overwrite=True) # TESTS

    avprof_alcol = np.ma.mean(stacked_2d, axis=0).data.copy()
    avprof_alrow = np.ma.mean(stacked_2d[:, 0:non], axis=1).data.copy()
    
    #ixnan = np.where(np.isnan(stacked_2d))
    #if len(ixnan[0])>0:
    #    stacked_2d.mask[ixnan] = True
    
    values = stacked_2d[:,0:non]
    gvalues = values.data[np.where(~values.mask)]
    
    stats_injection = dict(mean=np.mean(gvalues),
                       std=np.std(gvalues),
                       min=np.min(gvalues),
                       max=np.max(gvalues),
                       p5=np.percentile(gvalues, 5),
                       p25=np.percentile(gvalues, 25),
                       p50=np.percentile(gvalues, 50),
                       p75=np.percentile(gvalues, 75),
               p95=np.percentile(gvalues, 95))

    results = dict(avprof_alrow=avprof_alrow,
                   avprof_alcol=avprof_alcol, 
                   stats_injection=stats_injection)
    
    if debug:
        stop()
    
    return results


def get_spill(avprof_alcol,pattern):
    """ """
    non, noff, _ = pattern
    
    bgd = np.nanmean(avprof_alcol[non:])    
    y = avprof_alcol - bgd
    
    injection = np.nanmedian(y[0:non])
    maxspill = y[non]
    spill = maxspill / injection
    
    return spill

def msoftplus(IG1, a, xt):
    """ """
    return np.log10(1.+np.exp(-a*(IG1-xt)))

def invert_msoftplus(f,xt,a):
    return xt - np.log(10.**f-1.)/a

def der_msoftplus(IG1,xt,a):
    """ """
    der = -a/(np.log(10.))*(1./(1.+np.exp(a*(IG1-xt))))
    return der


def relu(IG1, a, xt):
    relu = np.zeros_like(IG1)
    relu[IG1<xt] = -a * (IG1[IG1<xt]-xt)
    return relu

def f_Inj_vs_IG1(IG1,b,k,xt,xN,a,N):
    """ """
    M = b + special.expit(k*(IG1-xt)) * (msoftplus(IG1,a,xN) + N)
    return M

def f_Inj_vs_IG1_ReLU(IG1, b, k, xt, xN, a, N):
    M = b + special.expit(k*(IG1-xt)) * (relu(IG1, a, xN) + N)
    return M


def redf_Inj_vs_IG1(IG1, xt, xN, a, N):
    bgd, drop = def_bgd_drop
    return f_Inj_vs_IG1(IG1, bgd, drop, xt, xN, a, N)

def redf_Inj_vs_IG1_ReLU(IG1, xt, xN, a, N):
    bgd, drop = def_bgd_drop
    return f_Inj_vs_IG1_ReLU(IG1, bgd, drop, xt, xN, a, N)


def f_Inj_vs_IDL(IDL,b,k,a,xt):
    """ """
    M = b + a * special.expit(-k*(IDL-xt))
    return M

def redf_Inj_vs_IDL(IDL,a,xt):
    """ """
    bgd, drop = def_bgd_drop
    return f_Inj_vs_IDL(IDL,bgd,drop,a,xt)


def fit_Inj_vs_IG1(IG1,med_inj,submodel='ReLU',doPlot=False,debug=False):
    """ """
    
    Npoints = len(IG1)
    IG1half = np.median(IG1)
    
    models_dict = dict(ReLU=f_Inj_vs_IG1_ReLU,
                       Softmax=f_Inj_vs_IG1,
                       reduced=dict(
                               ReLU=redf_Inj_vs_IG1_ReLU,
                               Softmax=redf_Inj_vs_IG1))
    
    if 4 <= Npoints <= 6:
        reduced = True
        fmodel = models_dict['reduced'][submodel]
        p0 = [IG1half, IG1half+3., 1., 0.01]
        
        bounds = ([2.,2.,0.01,0.],
                  [8.,10.,1.,0.2])
        
    elif Npoints > 6:
        reduced = False
        fmodel = models_dict[submodel]
        p0 = def_bgd_drop+[IG1half, IG1half +3., 1., 0.01]
        
        
        maxnmbgd = 200./2.**16
        bounds = ([-maxnmbgd,2.,2.,2.,0.01,0.],
                  [maxnmbgd,50.,8.,10.,1.,0.2])
        
    elif Npoints < 4:
        return dict(didfit=False)
    
    nmed_inj = med_inj / 2.**16  # Handy scaling
    
    xIG1 = np.linspace(IG1.min(),IG1.max(),1000)
    
        
    try: 
        popt, pcov = curve_fit(fmodel,IG1,nmed_inj,p0=p0,
                               method='trf',
                               bounds=bounds,
                               absolute_sigma=False)
        
        
        Inj_bf = fmodel(xIG1, *popt)
        
        didfit = True
        
    except RuntimeError:
        
        didfit = False
        
        popt = np.zeros(len(p0)) + np.nan
        Inj_bf = np.zeros_like(xIG1)
    
    
    if doPlot:
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(IG1,nmed_inj,'bo')
        ax.plot(xIG1,Inj_bf,'r--')
        plt.show()
    
    if reduced:
        arrsolution = def_bgd_drop+popt.tolist()
    else:
        arrsolution = copy.deepcopy(popt.tolist())
    
    solution = dict(zip(['BGD','K','XT', 'XN', 'A', 'N'],
                    arrsolution)
                    )
    solution['didfit'] = didfit
            
    solution['IG1_BF'] = xIG1.copy()
    solution['NORMINJ_BF'] = Inj_bf.copy()
    
    if debug:
        stop()
    
    return solution



def fit_Inj_vs_IDL(IDL,med_inj,doPlot=False,debug=False):
    """ """
    
    Npoints = len(IDL)
    IDLhalf = np.median(IDL)
    
    
    if 4 <= Npoints <= 6:
        reduced = True
        fmodel = redf_Inj_vs_IDL
        p0 = [0.1, IDLhalf]
        
        bounds = ([0.005,np.min(IDL)+0.25],
                  [1.,np.max(IDL)+0.25])
        
    elif Npoints > 6:
        reduced = False
        fmodel = f_Inj_vs_IDL
        p0 = def_bgd_drop+[0.1,IDLhalf]
        
        
        maxnmbgd = 200./2.**16
        bounds = ([-maxnmbgd,2.,0.005,np.min(IDL)+0.25],
                  [maxnmbgd,50.,10.,np.max(IDL)-0.25])
        
    elif Npoints < 4:
        return dict(didfit=False)
    
    nmed_inj = med_inj / 2.**16  # Handy scaling
    
    xIDL = np.linspace(IDL.min(),IDL.max(),1000)
    
        
    try:
        popt, pcov = curve_fit(fmodel,IDL,nmed_inj,p0=p0,
                               method='trf',
                               bounds=bounds,
                               absolute_sigma=False)
        
        Inj_bf = fmodel(xIDL, *popt)
        
        didfit = True
        
    except RuntimeError:
        
        didfit = False
        
        popt = np.zeros(len(p0)) + np.nan
        Inj_bf = np.zeros_like(xIDL)
    
    
    if doPlot:
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(IDL,nmed_inj,'bo')
        ax.plot(xIDL,Inj_bf,'r--')
        plt.show()
    
    if reduced:
        arrsolution = def_bgd_drop+popt.tolist()
    else:
        arrsolution = copy.deepcopy(popt.tolist())
    
    solution = dict(zip(['BGD','K','A','XT'],
                    arrsolution)
                    )
    solution['didfit'] = didfit
            
    solution['IDL_BF'] = xIDL.copy()
    solution['NORMINJ_BF'] = Inj_bf.copy()
    
    if debug:
        stop()
    
    return solution


def predict_inj_level(ID, IGs, id_timing, toi_ch, sectag):
    """ """

    IDL, IDH = ID
    IG1, IG2 = IGs
    id_wid, id_dly = id_timing

    inj_threshold = 7.5

    if id_wid < 10.:
        return 0.

    discrim = id_dly / float(toi_ch)
    if 2. < discrim < 3.:
        injsect = 'B'
    elif 1. < discrim < 2.:
        injsect = 'T'
    else:
        injsect = 'None'

    if injsect != sectag:
        return np.nan

    if (IDL < IG1+inj_threshold):
        injlevel = 600. + max(0, IG2-IG1)*6000.
    else:
        injlevel = np.nan

    return injlevel
