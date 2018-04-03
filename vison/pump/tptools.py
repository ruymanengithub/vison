#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""


Created on Fri Mar 16 14:38:51 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import copy
from scipy import  ndimage as nd
import pandas
from collections import  OrderedDict

from vison.datamodel import ccd as ccdmod

# END IMPORT

def get_injprofile_tpnorm(ccdobj,vstart,vend):
    """Produces a 2D Map of charge injection to be used in trap-pumping analysis,
    to obtain dipole maps."""
    
    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    NAXIS1 = ccdobj.NAXIS1
    NAXIS2 = ccdobj.NAXIS2
    
    ccdobj_dummy = copy.deepcopy(ccdobj)
    
    for Q in Quads:
        
        qimg = ccdobj.extract_region(Q,'img',vstart=vstart,vend=vend,canonical=True).copy()
        qimgmod = qimg.mean(axis=1).reshape(qimg.shape[0],1).repeat(qimg.shape[1],axis=1)
        
        qmod = np.ones((NAXIS1/2,NAXIS2/2),dtype='float32')
        
        qmod[prescan:-overscan,vstart:vend] = qimgmod.copy()
        
        ccdobj_dummy.set_quad(qmod,Q,canonical=True,extension=-1)
    
    injprofile_tpnorm = ccdobj_dummy.extensions[-1].data.copy()
    
    return injprofile_tpnorm

def gen_raw_dpmap_vtpump(ccdobj,Nrows=-1,vstart=0,vend=ccdmod.NrowsCCD,extension=-1):
    """ """
    
    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    
    for Q in Quads:
        qdata = ccdobj.get_quad(Q,canonical=True,extension=extension).copy()
        
        if Nrows == -1:
            avgrow = qdata[prescan:-overscan,vstart:vend].mean(axis=1)
            model = avgrow[:,None]
        else:
            model = nd.filters.median_filter(qdata[prescan:-overscan,vstart:vend],size=(1,Nrows),
                                             mode='reflect')
        
        dipmap = np.ones_like(qdata)
        dipmap[prescan:-overscan,vstart:vend] = qdata[prescan:-overscan,vstart:vend]/model
        
        ccdobj.set_quad(dipmap,Q,canonical=True,extension=extension)
        
    return ccdobj


def find_dipoles_vtpump(ccdobj,threshold,Q,vstart=0,vend=ccdmod.NrowsCCD,coosys='ccd',extension=-1):
    """ """
    
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    NAXIS2 = ccdobj.NAXIS2
    
    qrawmap = ccdobj.get_quad(Q,canonical=True,extension=extension)[prescan:-overscan,vstart:vend].copy()
    
    qraw_p1 = qrawmap[:,1:].copy()
    qraw_m1 = qrawmap[:,0:-1].copy()
    deltamap = qraw_p1-qraw_m1
    
    rows0,cols0 = np.where((np.abs(deltamap) > threshold*2.) & 
                           (np.abs(qraw_m1-1.)>threshold) &
                           (np.abs(qraw_p1-1.)>threshold) &
                           ((2.*np.abs(qraw_m1-1.)-np.abs(deltamap))<0.2*threshold) &
                           ((2.*np.abs(qraw_p1-1.)-np.abs(deltamap))<0.2*threshold))
    
    if len(rows0) == 0:
        return None
    
    if coosys == 'Q':
    
        X = rows0 + prescan 
        Y = cols0 + vstart 
    
    elif coosys == 'ccd':
        
        BB = ccdobj.QuadBound[Q]
        Xq = rows0 + prescan
        Yq = cols0 + vstart
        
        if Q in ['F','G']:
            BB = [BB[1]-1.,BB[0],BB[2],BB[3]]
        if Q in ['E','F']:
            BB = [BB[0],BB[1],BB[3]-1.,BB[2]]
        
        X = BB[0] + Xq * (BB[1]-BB[0])/np.abs(BB[1]-BB[0])
        Y = BB[2] + Yq * (BB[3]-BB[2])/np.abs(BB[3]-BB[2]) 
        
        
    S = np.zeros_like(X)
    S[np.where(deltamap[rows0,cols0]>0)] = 1
    S[np.where(deltamap[rows0,cols0]<0)] = 0
    A = deltamap[rows0,cols0]/2.
    
    indict = OrderedDict(X=X,Y=Y,S=S,A=A)
    df = pandas.DataFrame(indict,columns=['X','Y','S','A'])
         
    return df
    

def save_dipcat2D_as_ds9regs(df,Q,regfilename,clobber=True):
    """ """
    color = 'green'
    width = 2
    radius = 6.0
    
    hdr = ['# Region file format: DS9 version 4.1',
           'global color=%s dashlist=8 3 width=%i font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1' % \
            (color,width),
           'physical']
    body = []
    
    Nreg = len(df)
    for i in range(Nreg):
        body.append('circle(%.1f,%.1f,%.1f)' % (df['X'][i]+1.,df['Y'][i]+1.,radius))
    
    if clobber and os.path.exists(regfilename):
        os.system('rm %s' % regfilename)
    
    f = open(regfilename,'a')
    for line in hdr: print >> f, line
    for line in body: print >>f, line
    f.close()
    

    
    
    