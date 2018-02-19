#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary script to ccd.py

Created on Mon Feb 19 13:14:02 2018

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from astropy.io import fits as fts
import numpy as np
import os
from pdb import set_trace as stop
import sys
import datetime
import itertools
from collections import OrderedDict
import warnings

from scipy import ndimage as nd
from scipy import interpolate

from vison import __version__
# END IMPORT


def extract_region(ccdobj,Q,area='img',vstart=0,vend=2086,
                   Full=False,extension=-1):
    """ """
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    imgx = ccdobj.NAXIS1/2 - prescan - overscan
    
    Qdata = ccdobj.get_quad(Q,canonical=True,extension=extension)
    
    BB = [None,None,None,None]
    
    if area == 'pre':
        BB[0] = 0
        BB[1] = prescan
    elif area == 'img':
        BB[0] = prescan
        BB[1] = prescan + imgx
    elif area == 'ove':
        BB[0] = prescan + imgx
        BB[1] = None
    elif area == 'all':
        BB[0] = 0
        BB[1] = None
          
    BB[2] = vstart
    BB[3] = vend
    
    subregion = Qdata[BB[0]:BB[1],BB[2]:BB[3]].copy()
    
    if not Full:
        return subregion
    else:
        return subregion, BB


class Model2D():
    """Class for 2D models of images and images sections."""
    
    def __init__(self):
        """ """

def filter_img(img,filtsize=15,filtertype='median',Tests=False):
    """ """
    NX,NY = img.shape
    #xx,yy = np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
    
    
    if Tests:
        filtered = np.ones_like(img)
    else:
        if filtertype == 'median':    
            filtered = nd.median_filter(img,size=filtsize,mode='nearest') # 'constant',cval=0)
        elif filtertype == 'mean':
            filtered = nd.uniform_filter(img,size=filtsize,mode='nearest') # 'constant',cval=0)

    return filtered


def get_ilum_splines(img,sampling,splinemethod):
    """ """
    
    NX,NY = img.shape
    
    nX = NX/sampling
    nY = NY/sampling
    
    samplebin = 10
    
    sx = np.arange(samplebin/2,nX*samplebin+samplebin/2,samplebin)
    sy = np.arange(samplebin/2,nY*samplebin+samplebin/2,samplebin)
    
    sxx,syy = np.meshgrid(sx,sy,indexing='ij')
    
    zz = img[(sxx,syy)]
    
    xx,yy = np.mgrid[0:NX:NX*1j,0:NY:NY*1j]
    
    
    pilum = interpolate.griddata((sxx.flatten(),syy.flatten()),zz.flatten(),
                                (xx,yy), method=splinemethod,fill_value=np.nan)
    
    nans = np.isnan(pilum)
    pilum[nans] = img[nans].copy()
    
    
    ILUM = dict(img=img,polyfit=pilum,polycoefs=[])    
    
    return ILUM


def fit2Dpol_xyz(xx,yy,zz,degree=1):
    """ """
    from astropy.modeling import models, fitting
    
    p_init = models.Polynomial2D(degree=degree)
    fit_p = fitting.LinearLSQFitter()
    
    with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter (if changing fitter...)
        warnings.simplefilter('ignore')
        p = fit_p(p_init, xx, yy, zz)
    
    return p


def get_ilum_poly2D(img,sampling,pdegree=5):
    """ """
    
    NX,NY = img.shape
    
    nX = NX/sampling
    nY = NY/sampling
    
    x = np.arange(sampling/2,nX*sampling+sampling/2,sampling)
    y = np.arange(sampling/2,nY*sampling+sampling/2,sampling)
    
    xx,yy = np.meshgrid(x,y,indexing='ij')
    
    zz = img[(xx,yy)]
    
    p = fit2Dpol_xyz(xx,yy,zz,degree=pdegree)
    
    xp,yp= np.mgrid[:NX,:NY]
    pilum = p(xp, yp)
    
    ILUM = dict(img=img.coopy(),polyfit=pilum,polycoefs=p)
    
    return ILUM


def get_2Dmodel(ccdobj,Q,area='img',kind='spline',splinemethod='cubic',pdegree=2,
                doFilter=False,filtsize=1,filtertype='mean',
                vstart=0,vend=2086,extension=-1):
    """ """
    
    #prescan = ccdobj.prescan
    #overscan = ccdobj.overscan
    #imgx = ccdobj.NAXIS1/2 - prescan - overscan
    
    subregion, BB = ccdobj.extract_region(Q,area=area,vstart=vstart,vend=vend,
                                      Full=True,
                                          extension=extension)
    
    if doFilter:
        if filtsize > 1:
            subregion = filter_img(subregion,filtsize=filtsize,
                                   filtertype=filtertype)
    
    if kind == 'poly2D':
        ILUM = get_ilum_poly2D(subregion,sampling=filtsize,pdegree=pdegree)
    elif kind == 'spline':
        ILUM = get_ilum_splines(subregion,sampling=filtsize,
                                splinemethod=splinemethod)
    
    return ILUM
    
    
class Profile1D(object):
    """Class for 1D profiles of images and images sections."""
    
    def __init__(self,x,y):
        """ """
        assert len(x) == len(y)
        assert x.shape == y.shape
        
        self.data = OrderedDict()
        self.data['x'] = x.copy()
        self.data['y'] = y.copy()


def get_1Dprofile(ccdobj,Q,orient='hor',area='img',stacker='mean',vstart=0,
                  vend=2086,extension=-1):
    """ """
    
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    imgx = ccdobj.NAXIS1/2 - prescan - overscan
    
    subregion, BB = ccdobj.extract_region(Q,area=area,vstart=vstart,vend=vend,
                                      Full=True,
                                          extension=extension)
    
    if isinstance(subregion,np.ma.masked_array):
        stacker_dict = dict(mean=np.ma.mean,median=np.ma.median)
    else:
        stacker_dict = dict(mean=np.mean,median=np.median)
    
    if orient == 'hor': stackaxis = 1
    elif orient == 'ver': stackaxis = 0
    
    y = stacker_dict[stacker](subregion,axis=stackaxis)
    
    if orient == 'hor':
        if Q in ['E','H']:
            x = np.arange(BB[0],BB[1])
        elif Q in ['F','G']:
            x = np.arange(BB[1],BB[0],-1)+prescan+imgx+overscan
    elif orient == 'ver':
        if Q in ['G','H']:
            x = np.arange(BB[2],BB[3])
        elif Q in ['E','F']:
            x = np.arange(BB[3],BB[2],-1)+ccdobj.NAXIS2/2
       
    v1d = Profile1D(x,y)
    
    return v1d


    
    