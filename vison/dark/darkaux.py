#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary functions to BIAS01 and DARK01 test scripts.

Created on Tue Nov 14 13:25:01 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import copy
import numpy as np
from pdb import set_trace as stop
from collections import OrderedDict

from vison.datamodel.cdp import CDP as CDPClass
from vison import __version__
from vison.image import cosmetics
from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp
from vison.support.files import cPickleRead
# END IMPORT


def get_DarkDefectsMask_CDP(
        taskobj,
        darkccdobj,
        thresholds,
        subbgd=True,
        bgdmodel=None,
        extension=-1,
        sn_ccd='Unknown'):
    """ """

    if subbgd:
        if bgdmodel is None:

            reg2Dmodkwargs = dict(kind='poly2D', pdegree=5,
                doBin=True, binsize=200, recoveredges=True,
                filtertype='median', vstart=0, vend=2066)

            bgdmodel = cosmetics.get_bgd_model(darkccdobj, extension=extension,
                reg2Dmodkwargs=reg2Dmodkwargs)

        darkccdobj.sub_bias(bgdmodel, extension=extension)


    darkdata = darkccdobj.extensions[extension].data.copy()
    mask = cosmetics.get_Thresholding_DefectsMask(darkdata, thresholds)
    mask = cosmetics.mask_badcolumns(mask,colthreshold=200)

    data = OrderedDict(mask=copy.deepcopy(mask),
                        labels=['mask'])
    meta = OrderedDict(
                       MASKTYPE='DARK',
                       SN=sn_ccd,
                       THRESH=thresholds.__repr__(), 
                       SUBBGD=subbgd,
                       NORMED=False,
                       FUNCTION=cosmetics.get_Thresholding_DefectsMask.__name__)

    maskcdp = cdp.CCD_CDP(ID=taskobj.ID,
                      BLOCKID=taskobj.BLOCKID, 
                      CHAMBER=taskobj.CHAMBER)

    maskcdp.ingest_inputs(data=data, meta=meta)

    return maskcdp




def produce_MasterDark(outfits, infitsList=[], ccdpickList=[], mask=None, settings={}):
    """Produces a Master Dark out of a number of flat-illumination exposures."""

    assert (len(infitsList) == 0) != (len(ccdpickList) == 0)


    if len(ccdpickList) > 0:
        ccdobjList = [copy.deepcopy(cPickleRead(item)) for item in ccdpickList]
        ccdpile = ccdmod.CCDPile(ccdobjList=ccdobjList, extension=-1, withpover=True)
        NCOMB = len(ccdpickList)
    elif len(infitsList) > 0:
        ccdpile = ccdmod.CCDPile(
            infitsList=infitsList, extension=0, withpover=True)
        NCOMB = len(infitsList)

    stackimg, stackstd = ccdpile.stack(method='median', dostd=True)

    if isinstance(stackimg, np.ma.MaskedArray):
        stackmask = stackimg.mask.copy()
        stackimg = stackimg.data.copy()
        stackstd = stackstd.data.copy()
        if mask is not None and len(stackmask.shape)>0:
            mask = mask | stackmask
        if mask is not None and len(stackmask.shape)==0:
            pass
        else:
            mask = stackmask.copy()

    metabag = dict(CCDTempTop=np.nan,
                   CCDTempBot=np.nan,
                   CCDSerial='Unknown',
                   ID='0', BLOCKID='0', CHAMBER='None')
    metabag.update(settings)

    dkdata = dict(Dark=stackimg, eDark=stackstd)
    dkmeta = dict(NCOMB=NCOMB,
                  CCDTTOP=metabag['CCDTempTop'],
                  CCDTBOT=metabag['CCDTempBot'],
                  CCDSN=metabag['CCDSerial']
                  )

    mdk = DarkCDP(data=dkdata, meta=dkmeta, withpover=True,
                  ID=metabag['ID'], BLOCKID=metabag['BLOCKID'],
                  CHAMBER=metabag['CHAMBER'])

    if mask is not None and len(mask.shape) > 0:
        mdk.get_mask(mask)

    mdk.writeto(outfits, clobber=True, unsigned16bit=False)


class DarkCDP(ccdmod.CCD, CDPClass):
    """ """

    def __init__(
            self,
            fitsfile='',
            data=None,
            meta=None,
            withpover=True,
            ID=None,
            BLOCKID=None,
            CHAMBER=None):
        """ """

        self.BLOCKID = BLOCKID
        self.ID = ID
        self.CHAMBER = CHAMBER
        self.vison = __version__
        #self.Dark = None
        #self.eDark = None
        #self.Mask = None

        if data is None:
            data = dict()
        if meta is None:
            meta = dict()

        print('TODO: darkaux.DarkCDP needs improvemenents: masking')

        if fitsfile != '':
            super(DarkCDP, self).__init__(infits=fitsfile,
                                          getallextensions=True, 
                                          withpover=withpover)
            self.parse_fits()
        else:
            super(DarkCDP, self).__init__(infits=None, withpover=withpover)

            assert isinstance(data, dict)
            assert 'Dark' in list(data.keys())
            assert 'eDark' in list(data.keys())

            assert isinstance(meta, dict)
            assert 'NCOMB' in list(meta.keys())
            self.ncomb = meta['NCOMB']

            self.add_extension(data=None, header=None,
                               label=None, headerdict=meta)
            self.add_extension(data=data['Dark'].copy(), label='DARK')
            self.add_extension(data=data['eDark'].copy(), label='EDARK')

            if 'Mask' in list(data.keys()):
                self.add_extension(data=data['Mask'].copy(), label='MASK')

    def parse_fits(self,):
        """ """

        assert self.nextensions >= 3

        self.ncomb = self.extensions[0].header['NCOMB']
        assert self.extensions[1].label.upper() == 'DARK'
        #self.Dark = self.extensions[1].data
        assert self.extensions[2].label.upper() == 'EDARK'
        #self.eDark = self.extensions[2].data

        if self.nextensions > 3:
            assert self.extensions[3].label.upper() == 'MASK'
        #    self.Mask = self.extensions[3].data
        #else:
        #    self.Mask = None
