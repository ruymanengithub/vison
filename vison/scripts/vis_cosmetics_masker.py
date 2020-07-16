#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to create cosmetics masks in VIS Ground Calibration Campaign.

:History:
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
import numpy as np
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

isthere = os.path.exists


def read_OBSID_list(ff):
    """ """
    f = open(ff)
    OBSID_list = f.readlines()
    f.close()
    return [int(item) for item in OBSID_list if item != '']


def get_FITS_list(datapath, OBSID_list, CCD):
    FITS_list = []
    for i, OBSID in enumerate(OBSID_list):
        res = glob.glob(os.path.join(datapath, 'EUC_%i_*_ROE1_CCD%i.fits' % (OBSID, CCD)))
        FITS_list.append(res[0])
    return FITS_list


def pre_process(FITS_list, subOffset=False, validrange=None):
    """ """

    Quads = ['E', 'F', 'G', 'H']

    ccdobj_list = []

    for i, FITS in enumerate(FITS_list):
        ccdobj = ccdmod.CCD(FITS)

        if subOffset:
            for Q in Quads:
                ccdobj.sub_offset(Q, method='median', scan='pre',
                                  trimscan=[25, 5], ignore_pover=True,
                                  extension=-1)

        if validrange is not None:
            img = ccdobj.extensions[-1].data.copy()
            mask = np.zeros_like(img)
            mask[np.where((img < validrange[0]) | (img > validrange[1]))] = 1
            ccdobj.get_mask(mask)

        ccdobj_list.append(ccdobj)

    return ccdobj_list


def do_Mask(inputs, masktype, subbgd=True, normbybgd=False, validrange=None,
            flagWholeColumns=False):
    """ """

    assert subbgd != normbybgd

    onTests = False

    CCD = inputs['CCD']
    sn_ccd = inputs['sn_ccd']
    tag = inputs['tag']
    outpath = inputs['outpath']
    sp_inputs = inputs['inputs_%s' % masktype]
    OBSID_list_raw = sp_inputs['OBSID_list']
    datapath = sp_inputs['datapath']
    thresholds = sp_inputs['thresholds']

    outfilename = os.path.join(
        outpath, 'EUC_%s_%s_CCD%i_SN_%s.fits' %
        (masktype.upper(), tag, CCD, sn_ccd))

    if isinstance(OBSID_list_raw, list):
        OBSID_list = OBSID_list_raw
    elif isinstance(OBSID_list_raw, str):
        OBSID_list = read_OBSID_list(OBSID_list_raw)

    FITS_list = get_FITS_list(datapath, OBSID_list, CCD)

    # TESTS

    tmp_pickf = os.path.join(outpath, '%s_%s_%s_ccdobjs.pick' % (masktype, tag, sn_ccd))

    if not onTests:

        ccdobjList = pre_process(FITS_list, subOffset=True, validrange=validrange)

        # SUBTRACT Background

        if subbgd:

            for i, iccdobj in enumerate(ccdobjList):
                bgdmodel = cosmetics.get_bgd_model(iccdobj, extension=-1, margins=0.)
                iccdobj.sub_bias(bgdmodel, extension=-1)

                if iccdobj.masked:
                    for iext in range(iccdobj.nextensions):
                        iccdobj.extensions[iext].data = iccdobj.extensions[iext].data.data.copy()
                    iccdobj.masked = False  # HACK

        elif normbybgd:

            for i, iccdobj in enumerate(ccdobjList):
                bgdmodel = cosmetics.get_bgd_model(iccdobj, extension=-1, margins=1.)
                iccdobj.divide_by_flatfield(bgdmodel, extension=-1)

        files.cPickleDumpDictionary(dict(ccdobjs=ccdobjList), tmp_pickf)
    else:
        ccdobjList = files.cPickleRead(tmp_pickf)['ccdobjs']

    # MEDIAN STACK
    if len(ccdobjList) > 1:
        ccdpile = ccdmod.CCDPile(ccdobjList=ccdobjList, extension=-1)
        stacked = ccdpile.stack(method='median')
    else:
        stacked = ccdobjList[0].extensions[-1].data.copy()

    if isinstance(stacked, np.ma.masked_array):
        stacked = stacked.data.copy()

    # THRESHOLDING

    mask = cosmetics.get_Thresholding_DefectsMask(stacked, thresholds)

    # SETTING PRE/OVERSCANS TO ZERO (those can't have cosmetic defects)

    mask = cosmetics.set_extrascans(mask, val=0)

    # DISCARDING WHOLE COLUMNS [OPTIONAL]

    if flagWholeColumns:
        mask = cosmetics.mask_badcolumns(mask,colthreshold=200)

    # SAVING to a CDP

    data = OrderedDict()
    data['MASK'] = mask.copy()
    data['labels'] = ['MASK']
    meta = OrderedDict()
    meta['MASKTYPE'] = masktype
    meta['SN'] = sn_ccd
    meta['THRESH'] = thresholds.__repr__()
    meta['SUBBGD'] = subbgd
    meta['NORMED'] = normbybgd
    meta['FUNCTION'] = inspect.stack()[0][3]

    maskcdp = cdp.CCD_CDP(ID=vistime.get_time_tag(),
                          BLOCKID='Uknown',
                          CHAMBER='Unknown')

    maskcdp.ingest_inputs(data=data, meta=meta)

    maskcdp.savehardcopy(outfilename)

    os.system('rm %s' % tmp_pickf)

    return outfilename


def do_DarkMask(inputs):
    return do_Mask(inputs, 'dkmask', subbgd=True, normbybgd=False,
                   validrange=[0.0, 1.E6],
                   flagWholeColumns=True)


def do_FlatMask(inputs):
    return do_Mask(inputs, 'flmask', 
        validrange=[0.0, 1.E6],
        subbgd=False, 
        normbybgd=True,
        flagWholeColumns=True)


def do_MergeMasks(inputs):

    CCD = inputs['CCD']
    sn_ccd = inputs['sn_ccd']
    tag = inputs['tag']
    outpath = inputs['outpath']

    dkmask_f = os.path.join(outpath, 'EUC_DKMASK_%s_CCD%i_SN_%s.fits' % (tag, CCD, sn_ccd))
    flmask_f = os.path.join(outpath, 'EUC_FLMASK_%s_CCD%i_SN_%s.fits' % (tag, CCD, sn_ccd))
    assert isthere(dkmask_f) and isthere(flmask_f)

    outfilename = os.path.join(outpath, 'EUC_MASK_%s_CCD%i_SN_%s.fits' % (tag, CCD, sn_ccd))

    dkmaskobj = ccdmod.CCD(dkmask_f)
    flmaskobj = ccdmod.CCD(flmask_f)

    dkmask = dkmaskobj.extensions[-1].data.copy()
    flmask = flmaskobj.extensions[-1].data.copy()

    mask = dkmask.astype('int32') | flmask.astype('int32')

    # SAVING to a CDP

    data = OrderedDict()
    data['MASK'] = mask.copy()
    data['labels'] = ['MASK']
    meta = OrderedDict()
    meta['MASKTYPE'] = 'COMBO'
    meta['SN'] = sn_ccd
    meta['DKMASK'] = os.path.split(dkmask_f)[-1]
    meta['FLMASK'] = os.path.split(flmask_f)[-1]
    meta['FUNCTION'] = inspect.stack()[0][3]

    maskcdp = cdp.CCD_CDP(ID=vistime.get_time_tag(),
                          BLOCKID='Uknown',
                          CHAMBER='Unknown')

    maskcdp.ingest_inputs(data=data, meta=meta)

    maskcdp.savehardcopy(outfilename)

    return outfilename


def run_maskmaker(inputs):
    """ """

    outpath = inputs['outpath']
    if not isthere(outpath):
        os.system('mkdir -p %s' % outpath)

    outputs = dict()

    if inputs['doDarkMask']:
        print('Generating Defects in Darkness Mask...')
        dkmask_f = do_DarkMask(inputs)
        outputs['DARK'] = dkmask_f

    if inputs['doFlatMask']:
        print('Generating Defects in PR Mask...')
        flmask_f = do_FlatMask(inputs)
        outputs['FLAT'] = flmask_f

    if inputs['doMerge']:
        print('Generating Cosmetics Mask (Dark | PR)...')
        mgmask_f = do_MergeMasks(inputs)
        outputs['MERGE'] = mgmask_f

    return outputs


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-j", "--json", dest="json",
                      default='', help="json file with inputs")

    (options, args) = parser.parse_args()

    if options.json == '':
        parser.print_help()
        sys.exit()

    inputs = vjson.load_jsonfile(options.json, useyaml=True)

    _ = run_maskmaker(inputs)
