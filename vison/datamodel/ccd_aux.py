#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary script to ccd.py

Created on Mon Feb 19 13:14:02 2018

:author: raf

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
from scipy import signal

from vison import __version__
# END IMPORT


def rebin(arr, new_shape, stat='mean'):
    """"Rebin 2D array arr to shape new_shape by averaging."""
    if stat == 'mean':
        fgroup = np.mean
    elif stat == 'median':
        fgroup = np.median

    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])

    return fgroup(fgroup(arr.reshape(shape), axis=-1), axis=1)


def extract_region(ccdobj, Q, area='img', vstart=0, vend=2086,
                   Full=False, canonical=True, extension=-1):
    """ """
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    imgx = ccdobj.NAXIS1/2 - prescan - overscan

    # the quadrant is "extracted" always in canonical orientation,
    # to ease extraction of sub-area (pre-scan, over-scan, image-area)
    Qdata = ccdobj.get_quad(Q, canonical=True, extension=extension)

    BB = [None, None, None, None]

    if area == 'pre':
        BB[0] = 0
        BB[1] = prescan
    elif area == 'img':
        BB[0] = prescan
        BB[1] = prescan + imgx
    elif area == 'ove':
        BB[0] = prescan + imgx
        BB[1] = prescan + imgx + overscan
    elif area == 'all':
        BB[0] = 0
        BB[1] = prescan + imgx + overscan

    BB[2] = vstart
    BB[3] = vend

    subregion = Qdata[BB[0]:BB[1], BB[2]:BB[3]].copy()

    # if not canonical we return area in "DS9" orientation for quadrant

    if not canonical:
        if Q == 'E':
            subregion = subregion[:, ::-1].copy()
        elif Q == 'F':
            subregion = subregion[::-1, ::-1].copy()
        elif Q == 'G':
            subregion = subregion[::-1, :].copy()
        elif Q == 'H':
            subregion = subregion[:, :].copy()

    if not Full:
        return subregion
    else:
        return subregion, BB


class Model2D():
    """Class for 2D models of images and images sections."""

    def __init__(self, img, corners=[]):
        """ """
        assert isinstance(img, np.ndarray)
        self.img = img.copy()
        self.binnedimg = None
        self.XXbin = None
        self.YYbin = None

        if len(corners) > 0:
            NX = corners[1]-corners[0]+1
            NY = corners[3]-corners[2]+1
            assert (NY, NX) == img.shape

        self.corners = corners
        self.imgmodel = np.zeros_like(img, dtype='float32')
        self.polycoeffs = dict()
        self.polyfunc = None

    def bin_img(self, boxsize, stat='median'):
        """ """

        ishape = self.img.shape
        _nshape = (ishape[0]/boxsize, ishape[1]/boxsize)
        trimshape = _nshape[0]*boxsize, _nshape[1]*boxsize
        _img = self.img[0:trimshape[0], 0:trimshape[1]].copy()

        rebinned = rebin(_img, _nshape, stat=stat)

        sx = np.arange(boxsize/2, rebinned.shape[0]*boxsize+boxsize/2, boxsize)
        sy = np.arange(boxsize/2, rebinned.shape[1]*boxsize+boxsize/2, boxsize)

        sxx, syy = np.meshgrid(sx, sy, indexing='ij')

        self.binnedimg = rebinned.copy()
        self.XXbin = sxx.copy()
        self.YYbin = syy.copy()

    def filter_img(self, filtsize=15, filtertype='median', Tests=False):
        """ """

        if Tests:
            filtered = np.ones_like(self.img)
        else:
            if filtertype == 'median':
                # 'constant',cval=0)
                filtered = nd.median_filter(
                    self.img, size=filtsize, mode='nearest')
                #filtered = signal.medfilt2d(self.img,kernel_size=filtsize)
            elif filtertype == 'mean':
                filtered = nd.uniform_filter(
                    self.img, size=filtsize, mode='nearest')  # 'constant',cval=0)

        self.img = filtered.copy()

        return None

    def get_model_splines(self, sampling=1, splinemethod='cubic', useBin=False):
        """ """

        #corners = self.corners

        NX, NY = self.img.shape

        if useBin:

            sxx = self.XXbin.copy()
            syy = self.YYbin.copy()
            zz = self.binnedimg.copy()

        else:

            nX = NX/sampling
            nY = NY/sampling

            #samplebin = 10

            sx = np.arange(sampling/2, nX*sampling+sampling/2, sampling)
            sy = np.arange(sampling/2, nY*sampling+sampling/2, sampling)

            sxx, syy = np.meshgrid(sx, sy, indexing='ij')

            zz = self.img[(sxx, syy)]

        
        if isinstance(zz, np.ma.masked_array):
            selix = np.where(~zz.mask)
            sxx = sxx[selix]
            syy = syy[selix]
            zz = zz[selix]
            
        xx, yy = np.mgrid[0:NX:NX*1j, 0:NY:NY*1j]


        pilum = interpolate.griddata((sxx.flatten(), syy.flatten()), zz.flatten(),
                                     (xx, yy), method=splinemethod, fill_value=np.nan)

        nans = np.isnan(pilum)
        pilum[nans] = self.img[nans].copy()
        
        self.imgmodel = pilum.copy()

        return None

    def fit2Dpol_xyz(self, xx, yy, zz, degree=1):
        """ """
        from astropy.modeling import models, fitting

        p_init = models.Polynomial2D(degree=degree)
        fit_p = fitting.LinearLSQFitter()

        with warnings.catch_warnings():
            # Ignore model linearity warning from the fitter (if changing fitter...)
            warnings.simplefilter('ignore')
            p = fit_p(p_init, xx, yy, zz)

        return p

    def get_model_poly2D(self, sampling=1, pdegree=5, useBin=False):
        """ """

        NX, NY = self.img.shape

        if useBin:

            xx = self.XXbin.copy()
            yy = self.YYbin.copy()
            zz = self.binnedimg.copy()

        else:

            nX = NX/sampling
            nY = NY/sampling

            x = np.arange(sampling/2, nX*sampling+sampling/2, sampling)
            y = np.arange(sampling/2, nY*sampling+sampling/2, sampling)

            xx, yy = np.meshgrid(x, y, indexing='ij')

            zz = self.img[(xx, yy)]

        if isinstance(zz, np.ma.masked_array):
            selix = np.where(~zz.mask)
            xx = xx[selix]
            yy = yy[selix]
            zz = zz[selix]

        p = self.fit2Dpol_xyz(xx, yy, zz, degree=pdegree)

        xp, yp = np.mgrid[:NX, :NY]
        pilum = p(xp, yp)

        self.imgmodel = pilum.copy()
        self.polycoeffs = dict(zip(p.param_names, p.parameters))
        self.polyfunc = p

        return None


def get_region2Dmodel(ccdobj, Q, area='img', kind='spline', splinemethod='cubic', pdegree=2,
                      doFilter=False, doBin=True, filtsize=1, binsize=1, filtertype='mean',
                      vstart=0, vend=2086, canonical=True, extension=-1):
    """ """

    #prescan = ccdobj.prescan
    #overscan = ccdobj.overscan
    #imgx = ccdobj.NAXIS1/2 - prescan - overscan

    subregion, BB = ccdobj.extract_region(Q, area=area, vstart=vstart, vend=vend,
                                          Full=True, canonical=canonical,
                                          extension=extension)

    regmodel = Model2D(subregion)

    if doFilter:
        if filtsize > 1:
            regmodel.filter_img(filtsize=filtsize, filtertype=filtertype,
                                Tests=False)

    if doBin:
        regmodel.bin_img(boxsize=binsize, stat=filtertype)

    if kind == 'poly2D':
        regmodel.get_model_poly2D(
            sampling=filtsize, pdegree=pdegree, useBin=doBin)
    elif kind == 'spline':
        regmodel.get_model_splines(sampling=filtsize,
                                   splinemethod=splinemethod, useBin=doBin)

    return regmodel


class Profile1D(object):
    """Class for 1D profiles of images and images sections."""

    def __init__(self, x, y):
        """ """
        assert len(x) == len(y)
        assert x.shape == y.shape
        

        self.data = OrderedDict()
        self.data['x'] = x.copy()
        self.data['y'] = y.copy()


def get_1Dprofile(ccdobj, Q, orient='hor', area='img', stacker='mean', vstart=0,
                  vend=2086, extension=-1):
    """ """

    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    imgx = ccdobj.NAXIS1/2 - prescan - overscan

    subregion, BB = ccdobj.extract_region(Q, area=area, vstart=vstart, 
                                          vend=vend,
                                          Full=True,
                                          extension=extension)

    if isinstance(subregion, np.ma.masked_array):
        stacker_dict = dict(mean=np.ma.mean, median=np.ma.median)
    else:
        stacker_dict = dict(mean=np.mean, median=np.median)

    if orient == 'hor':
        stackaxis = 1
    elif orient == 'ver':
        stackaxis = 0

    y = stacker_dict[stacker](subregion, axis=stackaxis)

    if orient == 'hor':
        if Q in ['E', 'H']:
            x = np.arange(BB[0], BB[1])
        elif Q in ['F', 'G']:
            x = np.arange(BB[1], BB[0], -1)+prescan+imgx+overscan
    elif orient == 'ver':
        if Q in ['G', 'H']:
            x = np.arange(BB[2], BB[3])
        elif Q in ['E', 'F']:
            x = np.arange(BB[3], BB[2], -1)+ccdobj.NAXIS2/2
    v1d = Profile1D(x, y)

    return v1d
