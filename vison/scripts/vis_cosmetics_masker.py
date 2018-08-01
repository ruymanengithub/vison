#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to create cosmetics masks in VIS Ground Calibration Campaign.

Created on Wed Aug 1 11:02:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
from optparse import OptionParser
import sys
#from astropy.io import fits as fts
#import pyds9
#import itertools
import os
#from astropy.io import ascii
from collections import OrderedDict
#import numpy as np
#import copy
#import tempfile
#import subprocess
#import string as st
import glob
import inspect

from vison.datamodel import cdp
from vison.support import vistime
from vison.support import files
from vison.support import vjson
from vison.datamodel import ccd as ccdmod
from vison.image import cosmetics
# END IMPORT



def read_OBSID_list(ff):
    """ """
    f = open(ff)
    OBSID_list = f.readlines()
    f.close()
    return [int(item) for item in OBSID_list]

def get_FITS_list(datapath, OBSID_list, CCD):
    FITS_list = []
    for i, OBSID in enumerate(OBSID_list):

        res = glob.glob(os.path.join(datapath,'EUC_%i_*_ROE1_CCD%i.fits' % (OBSID,CCD)))
        FITS_list.append(res[0])
    return FITS_list

def pre_process(FITS_list, subOffset=False):
    """ """
    
    Quads = ['E', 'F', 'G', 'H']
    
    ccdobj_list = []
    
    for i, FITS in enumerate(FITS_list):
        ccdobj = ccdmod.CCD(FITS)
        
        if subOffset:
            for Q in Quads:
                ccdobj.sub_offset(Q, method='median', scan='ove',
                                  trimscan=[5,5], ignore_pover=True,
                                extension=-1)
        ccdobj_list.append(ccdobj)
        
    
    return ccdobj_list
    

def do_Mask(inputs, masktype, subbgd=True, normbybgd=False):
    """ """
    
    assert subbgd != normbybgd
    
    onTests = False
    
    CCD = inputs['CCD']
    sn_ccd = inputs['sn_ccd']
    tag = inputs['tag']
    sp_inputs = inputs['inputs_%s' % masktype]
    OBSID_list_f = sp_inputs['OBSID_list']
    datapath = sp_inputs['datapath']
    thresholds = sp_inputs['thresholds']
    
    outfilename = 'EUC_%s_%s_CCD_SN_%s.fits' % (masktype.upper(), tag, sn_ccd)
    
    OBSID_list = read_OBSID_list(OBSID_list_f)
    
    FITS_list = get_FITS_list(datapath, OBSID_list, CCD)
    
    # TESTS
    
    if not onTests:
    
        ccdobjList = pre_process(FITS_list, subOffset=True)
        
        # SUBTRACT Background
        
        if subbgd:
            
            for i, iccdobj in enumerate(ccdobjList):
                bgdmodel = cosmetics.get_bgd_model(iccdobj, extension=-1)
                iccdobj.sub_bias(bgdmodel, extension=-1)
        
        elif normbybgd:
            
            for i, iccdobj in enumerate(ccdobjList):
                bgdmodel = cosmetics.get_bgd_model(iccdobj, extension=-1)
                iccdobj.divide_by_flatfield(bgdmodel, extension=-1)
        
        files.cPickleDumpDictionary(dict(ccdobjs=ccdobjList),'%s_ccdobjs.pick' % masktype)
    else:
        ccdobjList = files.cPickleRead('%s_ccdobjs.pick' % masktype)['ccdobjs']
    
    # MEDIAN STACK
    if len(ccdobjList)>1:
        ccdpile = ccdmod.CCDPile(ccdobjList = ccdobjList, extension=-1)
        stacked = ccdpile.stack(method='median')
    else:
        stacked = ccdobjList[0].extensions[-1].data.copy()
    
    # THRESHOLDING
    
    mask = cosmetics.get_Thresholding_DefectsMask(stacked, thresholds)
    
    # SAVING to a CDP
    
    data = OrderedDict()
    data['MASK'] = mask.copy()
    data['labels'] = ['MASK']
    meta = OrderedDict()
    meta['MASKTYPE'] = masktype
    meta['SN']=sn_ccd
    meta['THRESH']=thresholds.__repr__()
    meta['SUBBGD'] = subbgd
    meta['NORMBYBGD'] = normbybgd
    meta['FUNCTION']=inspect.stack()[0][3]
    
    maskcdp = cdp.CCD_CDP(ID=vistime.get_time_tag(),
                      BLOCKID='Uknown',
                      CHAMBER='Unknown')
    
    maskcdp.ingest_inputs(data=data, meta=meta)
    
    maskcdp.savehardcopy(outfilename)
    
    return outfilename

def do_DarkMask(inputs):
    return do_Mask(inputs, 'dkmask', subbgd=True, normbybgd=False)

def do_FlatMask(inputs):
    return do_Mask(inputs, 'flmask', subbgd=False, normbybgd=True)

def do_MergeMasks(inputs):
    pass

def run_maskmaker(inputs):
    """ """
    
    if inputs['doDarkMask']:
        do_DarkMask(inputs)
    
    if inputs['doFlatMask']:
        do_FlatMask(inputs)
    
    if inputs['doMerge']:
        do_MergeMasks(inputs)
    
    


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-j", "--json", dest="json",
                      default='', help="json file with inputs")

    (options, args) = parser.parse_args()

    if options.json == '':
        parser.print_help()
        sys.exit()

    inputs = vjson.load_jsonfile(options.json, useyaml=True)
    
    run_maskmaker(inputs)

