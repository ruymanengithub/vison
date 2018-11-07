# -*- coding: utf-8 -*-
"""

Flat-fielding Utilities.


Created on Fri Apr 22 16:13:22 2016

@author: raf
"""

# IMPORT STUFF

from multiprocessing import Pool
from pdb import set_trace as stop
import numpy as np
import datetime
import warnings
from scipy import interpolate

#import vison

#from matplotlib import pyplot as plt

from vison.datamodel import ccd as ccdmodule
from vison.datamodel.cdp import CDP as CDPClass
from vison.support.files import cPickleRead
from vison import __version__

from scipy import ndimage as nd
#from scipy import signal
from astropy.io import fits as fts
# END IMPORT

Quads = ['E', 'F', 'G', 'H']


def produce_SingleFlatfield(infits, outfits, settings=None, runonTests=False):
    """ """
    #runonTests = False

    insettings = dict(kind='spline', splinemethod='cubic',
                      doBin=True,binsize=50,filtertype='median')
                      #doFilter=True, filtsize=50, filtertype='mean')
    
    if settings is not None:
        insettings.update(settings)
    

    ccdin = cPickleRead(infits)
    inwithpover = ccdin.withpover
    NrowsCCD = ccdin.NrowsCCD

    ccdout = ccdmodule.CCD(withpover=True)
    oshape = ccdout.shape

    ccdout.add_extension(np.ones(oshape, dtype='float32'), label='FLAT')
    ccdout.add_extension(np.zeros(oshape, dtype='float32'), label='MODEL')

    for Q in Quads:
        Qregmodel = ccdin.get_region2Dmodel(Q, area='img', vstart=0, vend=NrowsCCD,
                                            canonical=True, extension=-1,
                                            **insettings)
        
        # Set 1st Extension: image

        Qimg = ccdin.get_quad(Q, canonical=True, extension=-1).copy()
        QFF = np.ones_like(Qimg)
        QFF[ccdin.prescan:-ccdin.overscan,
            0:NrowsCCD] = Qimg[ccdin.prescan:-ccdin.overscan, 0:NrowsCCD].copy()

        ccdout.set_quad(QFF, Q, canonical=True, extension=0)
        
        # Set 2nd Extension: model

        QMod = np.ones(Qimg.shape, dtype='float32')
        QMod[ccdout.prescan:-ccdout.overscan,
             0:NrowsCCD] = Qregmodel.imgmodel.copy()

        ccdout.set_quad(QMod, Q, canonical=True, extension=1)

        # divide image by model
        
    ccdout.divide_by_flatfield(ccdout.extensions[1].data, extension=0)

    ccdout.writeto(outfits, clobber=True)
    
    return None


def _produce_SingleFlatfield(args):
    produce_SingleFlatfield(*args)


def produce_IndivFlats(infitsList, outfitsList, settings, runonTests, processes=6):
    """ """

    assert len(infitsList) == len(outfitsList)

    arglist = []

    for ix in range(len(infitsList)):
        arglist.append((infitsList[ix], outfitsList[ix], settings, runonTests))

    # generate flats using multiprocessing
    pool = Pool(processes=processes)
    #_produce_SingleFlatfield(arglist[0]) # TESTS
    pool.map(_produce_SingleFlatfield, arglist)


def produce_MasterFlat(infitsList, outfits, mask=None, settings={}):
    """Produces a Master Flat out of a number of flat-illumination exposures.
    Takes the outputs from produce_IndivFlats."""

    ccdpile = ccdmodule.CCDPile(
        infitsList=infitsList, extension=0, withpover=True)

    stackimg, stackstd = ccdpile.stack(method='median', dostd=True)

    metabag = dict(WAVEL=-1, ID='0', BLOCKID='0', CHAMBER='None',
                   CCDSerial='Uknown')
    metabag.update(settings)

    fdata = dict(Flat=stackimg, eFlat=stackstd)
    fmeta = dict(NCOMB=len(infitsList),
                 WAVEL=metabag['WAVEL'],
                 CCDSN=metabag['CCDSerial'])

    mff = FlatField(data=fdata, meta=fmeta, withpover=True,
                    ID=metabag['ID'], BLOCKID=metabag['BLOCKID'],
                    CHAMBER=metabag['CHAMBER'])

    if mask is not None:
        mff.get_mask(mask)
    

    mff.writeto(outfits, clobber=True, unsigned16bit=False)


class FlatField(ccdmodule.CCD, CDPClass):
    """ """

    def __init__(self, fitsfile='', data=None, meta=None, withpover=True, ID=None, BLOCKID=None, CHAMBER=None):
        """ """
        
        if data is None:
            data = dict()
        if meta is None:
            meta = dict()
            

        self.BLOCKID = BLOCKID
        self.ID = ID
        self.CHAMBER = CHAMBER
        self.vison = __version__

        # super(FlatField,self).__init__(infits=None)
        print 'TODO: FlatFielding.FlatField needs improvemenents: masking'

        if fitsfile != '':
            super(FlatField, self).__init__(infits=fitsfile,
                                            getallextensions=True, withpover=withpover)
            self.parse_fits()
        else:
            super(FlatField, self).__init__(infits=None, withpover=withpover)

            assert isinstance(data, dict)
            assert 'Flat' in data.keys()
            assert 'eFlat' in data.keys()

            assert isinstance(meta, dict)
            assert 'NCOMB' in meta.keys()
            self.ncomb = meta['NCOMB']
            assert 'WAVEL' in meta.keys()
            self.wavelength = meta['WAVEL']

            self.add_extension(data=None, header=None,
                               label=None, headerdict=meta)
            self.add_extension(data=data['Flat'].copy(), label='FLAT')
            self.add_extension(data=data['eFlat'].copy(), label='EFLAT')
            

            if 'Mask' in data.keys():
                self.add_extension(data=data['Mask'].copy(), label='MASK')

    def parse_fits(self,):
        """ """

        assert self.nextensions >= 3

        self.ncomb = self.extensions[0].header['NCOMB']
        self.wavelength = self.extensions[0].header['WAVEL']
        
        assert self.extensions[1].label.upper() == 'FLAT'
        self.Flat = self.extensions[1].data
        
        assert self.extensions[2].label.upper() == 'EFLAT'
        self.eFlat = self.extensions[2].data

        if self.nextensions > 3:
            assert self.extensions[3].label.upper() == 'MASK'
            self.Mask = self.extensions[3].data
        else:
            self.Mask = None
