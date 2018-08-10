#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module to provide common tools for analysis of Charge Injection acquisitions.

Created on Thu Sep 14 15:32:10 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop

from astropy.io import fits as fts
from pylab import plot,imshow,show
# END IMPORT


lineoffsets = dict(E=0, F=0, G=0, H=0)

def extract_injection_lines(ccdobj, Q, pattern, VSTART=0,
                            VEND=2066, suboffmean=False):
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
    
    if masked:
        stack_2d = np.ma.masked_array(stack_2d)

    for ii in (rowix):

        row = quaddata[npre:-1*npost, ii].copy()
        
        rowmedians[ii-VSTART] = np.nanmedian(row)
        rowstds[ii-VSTART] = np.nanstd(row)

        icycle = max((ii+lineoffset) / npercycle, 0)
        ix_in_cycle = (ii+lineoffset) - npercycle * icycle

        stack_2d[icycle, :, ix_in_cycle] = row.copy()
    
    if masked:
        stack_2d.mask[np.where(np.isnan(stack_2d.data))]=True

    if suboffmean:
        for icycle in range(nrep):
            stack_2d[icycle, :, :] -= np.nanmean(stack_2d[icycle, :, non:])
    
    
    stacked_2d = np.nanmean(stack_2d, axis=0)
    
    #mask2d = np.sum(stack_2d.mask,axis=0) # TESTS    
    #fts.writeto('stacked_2d.fits',stacked_2d.data.transpose(),overwrite=True) # TESTS
    #fts.writeto('mask2d.fits',mask2d.transpose().astype('float32'),overwrite=True) # TESTS

    avprof_alcol = np.nanmean(stacked_2d, axis=0)
    avprof_alrow = np.nanmean(stacked_2d[:, 0:non], axis=1)
    
    ixnan = np.where(np.isnan(stacked_2d))
    if len(ixnan[0])>0:
        stacked_2d.mask[ixnan] = True
    
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
    


def predict_inj_level(ID, IGs, id_timing, toi_ch, sectag):
    """ """

    IDL, IDH = ID
    IG1, IG2 = IGs
    id_wid, id_dly = id_timing

    inj_threshold = 7.5

    if id_wid < 10.:
        return 0.

    discrim = id_dly / float(toi_ch)
    if 1. < discrim < 2.:
        injsect = 'B'
    elif 2. < discrim < 3:
        injsect = 'T'
    else:
        injsect = 'None'

    if injsect != sectag:
        return np.nan

    if (IDL < IG1+inj_threshold):
        injlevel = 2000. + max(0, IG2-IG1)/0.5*2000.
    else:
        injlevel = np.nan

    return injlevel
