#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 11:55:12 2018

@author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
from collections import OrderedDict

from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp
# END IMPORT


def get_bgd_model(ccdobj, extension=-1, margins=None, 
        reg2Dmodkwargs=None):
    
    reg2Dmodkwargs_def = dict(area='img',
                        kind='poly2D', 
                        splinemethod='cubic',
                        pdegree=5,
                        doFilter=False, 
                        doBin=False,
                        filtsize=30, 
                        binsize=1,
                        filtertype='mean',
                        recoveredges=False,
                        vstart=0,
                        vend=2066,
                        canonical=True,
                        extension=extension)

    if reg2Dmodkwargs is not None:
        reg2Dmodkwargs_def.update(reg2Dmodkwargs)

    prescan = ccdobj.prescan

    bgd = ccdmod.CCD()
    bgd.add_extension(np.zeros_like(
        ccdobj.extensions[extension].data, dtype='float32'))

    for Quad in ccdobj.Quads:

        Qimgmodel = ccdobj.get_region2Dmodel(Quad, 
            **reg2Dmodkwargs_def).imgmodel.copy()

        # pre/overscans will not be modified unles margins is not None

        Qmodel = ccdobj.get_quad(
            Quad, canonical=True, extension=extension).copy()
        if margins is not None:
            Qmodel[:] = margins

        Qmodel[prescan:prescan + Qimgmodel.shape[0],
               0:Qimgmodel.shape[1]] = Qimgmodel.copy()

        bgd.set_quad(Qmodel, Quad, canonical=True, extension=-1)

    return bgd.extensions[-1].data.copy()


def get_Thresholding_DefectsMask(maskdata, thresholds):
    """ """
    #maskdata = ccdobj.extensions[extension].data
    mask = ((maskdata < thresholds[0]) | (maskdata > thresholds[1])).astype('int32')
    return mask


def set_extrascans(mask, val=0):
    """ """

    assert isinstance(mask,np.ndarray)

    soverscan = mask.shape[0]/2-ccdmod.prescan-ccdmod.NcolsCCD
    withpover = mask.shape[1]/2==(ccdmod.NrowsCCD+ccdmod.voverscan)

    mskccd = ccdmod.CCD(withpover=withpover,overscan=soverscan)
    mskccd.add_extension(data=mask)

    Qshape = mskccd.shape[0]/2, mskccd.shape[1]/2

    prescan = mskccd.prescan

    for Q in mskccd.Quads:

        Qmaskimg = mskccd.extract_region(Q, area='img', 
            vstart=0, vend=ccdmod.NrowsCCD,
            Full=False, canonical=True,
            extension=-1)
        
        Qmask = np.zeros(Qshape,dtype='int32')

        Qmask[prescan:prescan + Qmaskimg.shape[0],
               0:Qmaskimg.shape[1]] = Qmaskimg.copy()

        mskccd.set_quad(Qmask, Q, canonical=True, extension=-1)

    mask = mskccd.extensions[-1].data.copy()

    return mask

def mask_badcolumns(mask,colthreshold=200):
    """Flags entire column of pixels if N>colthreshold pixels in column are bad. """

    assert isinstance(mask,np.ndarray)

    soverscan = mask.shape[0]/2-ccdmod.prescan-ccdmod.NcolsCCD
    withpover = mask.shape[1]/2==(ccdmod.NrowsCCD+ccdmod.voverscan)

    mskccd = ccdmod.CCD(withpover=withpover,overscan=soverscan)
    mskccd.add_extension(data=mask)

    for Q in mskccd.Quads:

        quaddata = mskccd.get_quad(Q,canonical=False,extension=0)

        bX, bY = np.where(quaddata)
        uCol, uColcounts = np.unique(bX, return_counts=True)

        if len(uCol)>0:
            for ix,col in enumerate(uCol):
                if uColcounts[ix]>=colthreshold:
                    quaddata[col,:] = 1

        mskccd.set_quad(quaddata,Quadrant=Q,canonical=False,extension=0)

    maskout = mskccd.extensions[0].data.copy()

    return maskout