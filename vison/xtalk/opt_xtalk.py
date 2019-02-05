#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Module for analysis of Cross-Talk data obtained optical stimulation.

Created on Tue Jan 29 17:33:23 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import numpy as np
import copy
import os
from astropy.io import fits as fts

from vison.image import sextractor as sex
from vison.support import context
from vison.datamodel import ccd
from vison.support.files import cPickleRead
# END IMPORT

CCDs = [1,2,3]
Quads = context.Quads 

def exe_SEx(ccdobj,tag):
    """ """
    
    sexconfig = dict(MINAREA=2.,
             DET_THRESH=14.,
             MAG_ZERPOINT=20.,
             SATUR_LEVEL=65535.,
             SEEING_FWHM=1.,
             PIXEL_SCALE=1.,
             GAIN=1.
             )
    
    res = OrderedDict()
    
    for Q in Quads:
        
        
        res[Q]=OrderedDict()
        
        hdr = ccdobj.extensions[-1].header
        vstart = int(hdr['VSTART'])
        vend = int(hdr['VEND'])
        
        img = ccdobj.get_quad(Q,canonical=True,extension=-1)
        
        img = img[51:-20,vstart:vend]
        bgd = np.nanmedian(img)
        img -= bgd
        
        SExCatroot = '%s_%s' % (tag,Q)
        SExCat = sex.easy_run_SEx(img, SExCatroot, 
                     sexconfig=sexconfig,                                      
                     cleanafter=False)
        os.system('rm %s_BACK.fits' % SExCatroot)
        
        segm = fts.getdata('%s_SEGM.fits' % SExCatroot)
        ixsources = np.where(segm>0)
        ADUsources = img[ixsources]
        
        os.system('rm %s_SEGM.fits' % SExCatroot)
        
        res[Q]['ixsources'] = ixsources
        res[Q]['ADU'] = ADUsources.copy()
        res[Q]['SExCat'] = copy.deepcopy(SExCat)
        res[Q]['bgd'] = bgd
    
    return res


def get_rawct_plotdict(rawcrosstalks,CCDso,Qso):
    """ """
    
    
    CCDsok = 'CCD%i' % CCDso
    
    rct_plot = OrderedDict()
    
    emptyx = dict()
    emptyy = dict()
    for key in ['OPTICAL']:
        emptyx[key] = [None]# [0,2**16]
        emptyy[key] = [None] #[0,2**16]
    
    
    for CCD in CCDs:
        CCDk = 'CCD%i' % CCD
        rct_plot[CCDk] = OrderedDict()
        
        for Q in Quads:
            
            rct_plot[CCDk][Q] = OrderedDict()
            rct_plot[CCDk][Q]['x'] = OrderedDict()
            rct_plot[CCDk][Q]['y'] = OrderedDict()
            
            if CCD==CCDso and Q==Qso:
                rct_plot[CCDk][Q]['x'] = emptyx.copy()
                rct_plot[CCDk][Q]['y'] = emptyy.copy()
                continue
            else:
                _cdict = rawcrosstalks[CCDsok][Qso][CCDk][Q]
                
                rct_plot[CCDk][Q]['x']['OPTICAL'] = copy.deepcopy(_cdict['x'])
                rct_plot[CCDk][Q]['y']['OPTICAL'] = copy.deepcopy(_cdict['y'])
                
    
    rct_plot['labelkeys'] = ['OPTICAL']
    
    
    return rct_plot


def update_victim_rawcrosstalks(rawcrosstalks,imgv,sexdict,CCDv,Qv):
    """ """
    
    CCDvk = 'CCD%i' % CCDv
    
    sexvic = sexdict[CCDvk][Qv].copy()
    
    for CCDso in CCDs:
        for Qso in Quads:
            if CCDso==CCDv and Qso==Qv:
                continue
            CCDsok = 'CCD%i' % CCDso
            
            sexsou = sexdict[CCDsok][Qso].copy()
            
            mask = np.zeros_like(imgv)
            
            ixvic = sexvic['ixsources']
            ixsou = sexsou['ixsources']
            
            mask[ixsou] = 1
            mask[ixvic] -= 1
            
            imgsou = np.zeros_like(imgv)
            imgsou[ixsou] = sexsou['ADU']
                
            ixsel = np.where(mask>0)
                
            foox = imgsou[ixsel].copy().tolist()
            fooy = imgv[ixsel].copy().tolist()
            
            rawcrosstalks[CCDsok][Qso][CCDvk][Qv]['x'] += foox
            rawcrosstalks[CCDsok][Qso][CCDvk][Qv]['y'] += fooy
    
    return rawcrosstalks




def get_rawcrosstalk_mx(fitsdict):
    
    obsids= fitsdict['obsids']
    
    rawcrosstalks = OrderedDict()
    
    for CCDso in CCDs:
        CCDsk = 'CCD%i' % CCDso
        rawcrosstalks[CCDsk] = OrderedDict()
        for Qs in Quads:
            rawcrosstalks[CCDsk][Qs] = OrderedDict()
            for CCDv in CCDs:
                CCDvk = 'CCD%i' % CCDv
                rawcrosstalks[CCDsk][Qs][CCDvk] = OrderedDict()
                for Qv in Quads:
                    rawcrosstalks[CCDsk][Qs][CCDvk][Qv] = OrderedDict()
                    rawcrosstalks[CCDsk][Qs][CCDvk][Qv]['x'] = []
                    rawcrosstalks[CCDsk][Qs][CCDvk][Qv]['y'] = []
    
    
    for iobs,obsid in enumerate(obsids):
        
        print('\nGetting cross-talks from obsid %i/%i\n' % (iobs+1,len(obsids)))
        
        sexdict = cPickleRead(os.path.join(respath,'EUC_%i.pick' % obsid))
                
        for CCDv in CCDs:
            
            ifits = str(fitsdict['CCD%i' % CCDv][iobs])
            
            iccdobj = ccd.CCD(ifits)
            
            for Qv in Quads:
                
                hdr = iccdobj.extensions[-1].header
                vstart = int(hdr['VSTART'])
                vend = int(hdr['VEND'])
                
                imgv = iccdobj.get_quad(Qv,canonical=True,extension=-1)
                imgv = imgv[51:-20,vstart:vend]
                
                bgd = sexdict['CCD%i' % CCDv][Qv]['bgd']
                
                imgv -= bgd
                
                rawcrosstalks = update_victim_rawcrosstalks(rawcrosstalks,imgv,sexdict,CCDv,Qv)

    
    return rawcrosstalks


