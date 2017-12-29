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

# END IMPORT



def extract_injection_lines(quaddata,pattern,VSTART=1,
            VEND=2066,suboffmean=False,lineoffset=0):
    """     
    quaddata: quadrant data, array
    pattern: non,noff,nrep (lines on, off, repeatitions)
    VSTART: VSTART
    VEND: VEND
    suboffmean: bool, subtract median of non-injected lines
    lineoffset: integer, to account for a shift between readout lines and charge
                   injection pattern.
    """
    
    npre = 51
    npost = 20
    
    non,noff,nrep = pattern
    
    npercycle = non+noff
    
    NX = quaddata.shape[0]
    
    nlines = VEND-VSTART+1
    
    rowmedians = np.zeros(nlines,dtype='float32')
    rowstds = np.zeros(nlines,dtype='float32')
    rowix = np.arange(VSTART,VEND+1)
    
    stack_2d = np.zeros((nrep,NX-(npre+npost),npercycle),dtype='float32') + np.nan
    
    for ii in (rowix-1):
        
        row = quaddata[npre:-1*npost,ii].copy()
        
        rowmedians[ii-VSTART+1] = np.nanmedian(row)
        rowstds[ii-VSTART+1] = np.nanmedian(row)
        
        icycle = max((ii+lineoffset+1) / npercycle,0)
        ix_in_cycle = (ii+lineoffset) - npercycle * icycle
        
        stack_2d[icycle,:,ix_in_cycle] = row.copy()
    
    if suboffmean:
        for icycle in range(nrep):
            stack_2d[icycle,:,:] -= np.nanmean(stack_2d[icycle,:,non:])    
    
    stacked_2d = np.nanmean(stack_2d,axis=0)
    

    avprof_alcol = np.nanmean(stacked_2d[:,:],axis=0)
    avprof_alrow = np.nanmean(stacked_2d[:,0:non],axis=1)
        
    avinjection = np.mean(stacked_2d[:,0:non])
    
    stats_injection = [np.nanmedian(stacked_2d[:,0:non]),
                       np.nanstd(stacked_2d[:,0:non]),
                       np.percentile(stacked_2d[:,0:non],5),
                       np.percentile(stacked_2d[:,0:non],95)]
    
    results = dict(avinjection=avinjection,avprof_alrow=avprof_alrow,
                   avprof_alcol=avprof_alcol,stats_injection=stats_injection)
    
    return results
