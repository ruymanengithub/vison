# -*- coding: utf-8 -*-
"""

Data model for Euclid-VIS CCDs (ground testing at MSSL)

Created on Fri Nov 13 17:42:36 2015

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from astropy.io import fits as fts
import numpy as np
import os
from pdb import set_trace as stop
import sys
import collections
import datetime
import itertools
import ccd_aux
import ccdsim
from vison import __version__
# END IMPORT

isthere = os.path.exists

NAXIS1 = 4238
NrowsCCD = 2066
NcolsCCD = 2048
NAXIS2 = (NrowsCCD+20)*2  # 4132
prescan = 51
overscan = 20
voverscan = 20
#imgarea = [2048,2066]
RON = 1.4
gain = 3.1  # e/adu


QuadBound = dict(E=[0, NAXIS1/2, NAXIS2/2, NAXIS2],
                 F=[NAXIS1/2, NAXIS1, NAXIS2/2, NAXIS2],
                 G=[NAXIS1/2, NAXIS1, 0, NAXIS2/2],
                 H=[0, NAXIS1/2, 0, NAXIS2/2])


def f_get_QuadBound(NAXIS1, NAXIS2):
    return dict(E=[0, NAXIS1/2, NAXIS2/2, NAXIS2],
                F=[NAXIS1/2, NAXIS1, NAXIS2/2, NAXIS2],
                G=[NAXIS1/2, NAXIS1, 0, NAXIS2/2],
                H=[0, NAXIS1/2, 0, NAXIS2/2])


Quads = ['E', 'F', 'G', 'H']


class Extension():
    """Extension Class"""

    def __init__(self, data, header=None, label=None, headerdict=None):
        """ """

        self.data = data

        if header is None:
            header = fts.Header()

        if headerdict is not None:
            for key in headerdict:
                header[key] = headerdict[key]

        self.header = header

        if label is not None:
            self.header['EXTNAME'] = label

        self.label = label


def cooconv_arrays_decorate(func):
    """ """
    def cooconv_wrapper(*args):

        cooargs = args[1:]
        arearraysvec = []
        for cooarg in cooargs:
            arearraysvec.append(isinstance(
                cooargs[0], (list, tuple, collections.Sequence, np.ndarray)))
        assert np.all(np.array(arearraysvec) == arearraysvec[0])
        arearrays = arearraysvec[0]

        if arearrays:
            lengths = [len(cooarg) for cooarg in cooargs]
            assert np.all(np.array(lengths) == lengths[0])
            retlist = []
            for i in range(lengths[0]):
                retlist.append(
                    func(args[0], *tuple([cooarg[i] for cooarg in cooargs])))
            return tuple([np.array(item) for item in zip(*retlist)])
        else:
            return func(*args)

    return cooconv_wrapper


class CCD(object):
    """Class of CCD objects. 
    Euclid Images as acquired by ELVIS software (Euclid LabView Imaging Software).


    The class has been extended to handle multi-extension images. This is useful
    to also "host" calibration data-products, such as Flat-Fields.


    A note on Coordinates Systems:
        - 'CCD': referenced to the first pixel readout from channel H. All 4 quadrants
        in a single array, their detection nodes in the 4 "corners" of the 
        rectangle. Same system as images are displayed on DS9. In clock-wise
        sense, quadrants are H (bottom-left), E (top-left), F (top-right),
        and G (bottom-right).
        - 'Quadrant-canonical': Quadrant coordinates system in which the first pixel 
        is the first pixel read out (closest pixel to the readout node), and 
        the last is the last readout. In this system, the serial pre-scan comes 
        before the image area, and this before the serial overscan. Parallel 
        overscan comes after image area in the parallel direction. In this 
        system, coordinates of pixels across quadrants, for a single readout, 
        correspond to the same point in time. Useful when doing cross-talk analysis, 
        for example.
        - 'Quadrant-relative': quadrant coordinates system with the same relative orientation
        as in the 'CCD' system, but referenced to the 'lower-left' pixel of the
        given quadrant in such system. In this system, the readout node is in a 
        different corner for each quadrant: lower-left for H, top-left for E,
        top-right for F and bottom-right for G.

    """

    NrowsCCD = NrowsCCD
    NcolsCCD = NcolsCCD
    Quads = Quads

    get_1Dprofile = ccd_aux.get_1Dprofile
    get_region2Dmodel = ccd_aux.get_region2Dmodel
    extract_region = ccd_aux.extract_region

    simadd_flatilim = ccdsim.simadd_flatilum
    simadd_points = ccdsim.simadd_points
    simadd_bias = ccdsim.simadd_bias
    simadd_ron = ccdsim.simadd_ron
    simadd_poisson = ccdsim.simadd_poisson
    sim_window = ccdsim.sim_window
    simadd_injection = ccdsim.simadd_injection

    #rebin = ccd_aux.rebin

    def __init__(self, infits=None, extensions=[-1], getallextensions=False, withpover=True):
        """ """

        self.extnames = []
        self.extensions = []

        if infits is not None:

            assert type(
                infits) is str, "'%s' can't be a name for a file!" % infits
            assert isthere(infits), 'infits:%s is just not there :-(' % infits

            self.loadfromFITS(infits, extensions, getallextensions)

        else:

            self.extensions = []
            self.extnames = []

        self.nextensions = len(self.extensions)

        self.NAXIS1 = NAXIS1
        if withpover:
            self.NAXIS2 = NAXIS2
        else:
            self.NAXIS2 = NAXIS2-voverscan*2

        self.shape = (self.NAXIS1, self.NAXIS2)
        self.wQ = self.NAXIS1/2
        self.hQ = self.NAXIS2/2

        for iext in range(self.nextensions):
            if self.extensions[iext].data is not None:
                assert self.shape == self.extensions[iext].data.shape

        self.chinjlines = 2
        self.prescan = prescan
        self.overscan = overscan
        self.withpover = withpover
        if self.withpover:
            self.voverscan = voverscan
        else:
            self.voverscan = 0

        self.gain = dict(E=3.1, F=3.1, G=3.1, H=3.1)
        self.rn = dict(E=4.5, F=4.5, G=4.5, H=4.5)

        self.QuadBound = f_get_QuadBound(self.NAXIS1, self.NAXIS2)

        self.wQphys = self.NAXIS1/2-self.prescan-self.overscan
        self.hQphys = self.NAXIS2/2-self.voverscan+self.chinjlines
        self.QuadBoundPhys = f_get_QuadBound(2*self.wQphys, 2*self.hQphys)

        self.masked = False

        self.historial = []

        self.Qmatrix = np.array([['H', 'G'],
                                 ['E', 'F']])

#        self.sectors = dict(LL='H',
#                        LR='G',
#                        UL='E',
#                        UR='F')

        self.xCCD2Physoffsets = dict(E=self.prescan,
                                     F=2*self.overscan+self.prescan,
                                     G=2*self.overscan+self.prescan,
                                     H=self.prescan)

        self.yCCD2Physoffsets = dict(E=2*self.voverscan-2*self.chinjlines,
                                     F=2*self.voverscan-2*self.chinjlines,
                                     G=0,
                                     H=0)

    @cooconv_arrays_decorate
    def cooconv_Qrel_2_CCD(self, x, y, Q):
        """Converts coordiates from Quadrant-relative to CCD."""
        BB = self.QuadBound[Q]
        X = x + BB[0]
        Y = y + BB[2]

        return X, Y

    @cooconv_arrays_decorate
    def cooconv_Qcan_2_CCD(self, x, y, Q):
        """Converts coordiates from Quadrant-canonical to CCD."""
        xr, yr = self.cooconv_Qcan_2_Qrel(x, y, Q)
        X, Y = self.cooconv_Qrel_2_CCD(xr, yr, Q)
        return X, Y

    @cooconv_arrays_decorate
    def cooconv_CCD_2_Qrel(self, x, y, Q):
        """Converts coordinates from CCD to Quadrant-relative."""

        BB = self.QuadBound[Q]
        X = x - BB[0]
        Y = y - BB[2]
        return X, Y

    # @cooconv_arrays_decorate
    def get_Q(self, x, y, w, h):
        """ """
        xix = int(x >= w/2)
        yix = int(y >= h/2)
        return self.Qmatrix[yix, xix]

    @cooconv_arrays_decorate
    def get_Q_of_CCDcoo(self, x, y):
        """ """
        w = self.NAXIS1
        h = self.NAXIS2
        return self.get_Q(x, y, w, h)

    @cooconv_arrays_decorate
    def get_Q_of_Physcoo(self, x, y):
        """ """
        #assert isinstance(x,type(y))
        #arearrays = not np.isscalar(x)
        # if arearrays:
        #    return map(self.get_Q_of_Physcoo,zip(x,y))

        w = self.NAXIS1-self.prescan*2-self.overscan*2
        h = self.NAXIS2-self.voverscan*2+self.chinjlines*2
        return self.get_Q(x, y, w, h)

    @cooconv_arrays_decorate
    def cooconv_CCD_2_Qcan(self, x, y, Q):
        """Converts coordiates from CCD to Quadrant-canonical."""
        xp, yp = self.cooconv_CCD_2_Qrel(x, y, Q)
        X, Y = self.cooconv_Qrel_2_Qcan(xp, yp, Q)
        return X, Y

    @cooconv_arrays_decorate
    def cooconv_Qrel_2_Qcan(self, x, y, Q):
        """Converts coordiates from Quadrant-relative to Quadrant-canonical."""
        BB = self.QuadBound[Q]
        xzero = 0.
        yzero = 0.
        xsign = 1.
        ysign = 1.
        if Q in ['F', 'G']:
            xzero = BB[1]-BB[0]-1.
            xsign = -1.
        if Q in ['E', 'F']:
            yzero = BB[3]-BB[2]-1.
            ysign = -1.
        return xsign*x + xzero, ysign*y + yzero

    @cooconv_arrays_decorate
    def cooconv_Qcan_2_Qrel(self, x, y, Q):
        """Converts coordiates from Quadrant-canonical to Quadrant-relative."""
        return self.cooconv_Qrel_2_Qcan(x, y, Q)

    @cooconv_arrays_decorate
    def cooconv_CCD_2_Phys(self, x, y):
        """ """
        Q = self.get_Q_of_CCDcoo(x, y)

        xcan, ycan = self.cooconv_CCD_2_Qcan(x, y, Q)

        isonSi = ((xcan >= self.prescan) & (xcan < (self.wQ-self.overscan)) &
                  (ycan >= 0) & (ycan < self.hQ-self.voverscan))

        xp, yp = np.nan, np.nan
        if isonSi:
            xp = x - self.xCCD2Physoffsets[Q]
            yp = y - self.yCCD2Physoffsets[Q]
        return xp, yp

    @cooconv_arrays_decorate
    def cooconv_Phys_2_CCD(self, x, y):
        """ """
        Q = self.get_Q_of_Physcoo(x, y)

        isonSi = ((x >= 0) & (x < self.wQphys*2) &
                  False == ((y >= self.hQphys) & (y <= self.hQphys + self.chinjlines-1)))
        xp, yp = np.nan, np.nan
        if isonSi:
            xp = x + self.xCCD2Physoffsets[Q]
            yp = y + self.yCCD2Physoffsets[Q]
        return xp, yp

    def cooconvert(self, x, y, insys, outsys, Q='U'):
        """Coordinates conversion between different systems."""

        conversion_dict = dict(CCD=dict(Qrel=self.cooconv_CCD_2_Qrel, Qcan=self.cooconv_CCD_2_Qcan),
                               Qrel=dict(CCD=self.cooconv_Qrel_2_CCD,
                                         Qcan=self.cooconv_Qrel_2_Qcan),
                               Qcan=dict(Qrel=self.cooconv_Qcan_2_Qrel, CCD=self.cooconv_Qcan_2_CCD))

        converter = conversion_dict[insys][outsys]

        return converter(x, y, Q)

    def dummyrebin(self, arr, new_shape, stat='median'):
        """ """
        return ccd_aux.rebin(arr, new_shape, stat)

    def loadfromFITS(self, fitsfile, extensions=[-1], getallextensions=False):
        """Loads contents of self from a FITS file."""

        hdulist = fts.open(fitsfile)

        nextensions = len(hdulist)

        if getallextensions:
            extensions = np.arange(nextensions)

        for iext in extensions:

            hdu = hdulist[iext]

            if hdu.data is not None:
                data = hdu.data.transpose().astype('float32').copy()
            else:
                data = None
            header = hdu.header

            if 'EXTNAME' in hdu.header:
                label = hdu.header['EXTNAME']
            else:
                label = None

            self.extensions.append(Extension(data, header, label))
            self.extnames.append(label)

        hdulist.close()

    def add_extension(self, data, header=None, label=None, headerdict=None):
        """Appends an extension to self (extensions are in a list)."""
        if data is not None:
            assert data.shape == self.shape

        self.extensions.append(Extension(data, header, label, headerdict))
        self.nextensions += 1

    def del_extension(self, extension):
        """Deletes an extension from self, by index."""

        self.extensions.pop(extension)
        self.nextensions -= 1

    def set_extension(self, data, header=None, label=None, headerdict=None, extension=-1):
        """Sets extension 'extension' in self."""
        assert data.shape == self.shape
        self.extensions[extension] = Extension(data, header, label, headerdict)

    def get_quad(self, Quadrant, canonical=False, extension=-1):
        """Returns a quadrant in canonical or non-canonical orientation.

        :param Quadrant: Quadrant, one of 'E', 'F', 'G', 'H'
        :type Quadrant: char

        :param canonical: 

        Canonical [True] = with readout-node at pixel index (0,0) regardless of quadrant. 
        This is the orientation which corresponds to the data-reading order (useful 
        for cross-talk measurements, for example). Non-Canonical [False] = with 
        readout-node at corner matching placement of quadrant on the CCD. 
        This is the orientation that would match the representation of the image on DS9.        

        :type canonical: bool

        :param extension: extension number. Default = -1 (last)
        :type extension: int

        """

        edges = self.QuadBound[Quadrant]
        Qdata = self.extensions[extension].data[edges[0]:edges[1], edges[2]:edges[3]]

        if canonical:
            if Quadrant == 'E':
                Qdata = Qdata[:, ::-1].copy()
            elif Quadrant == 'F':
                Qdata = Qdata[::-1, ::-1].copy()
            elif Quadrant == 'G':
                Qdata = Qdata[::-1, :].copy()
            elif Quadrant == 'H':
                Qdata = Qdata[:, :].copy()

        return Qdata

    def get_tile_coos(self, Quadrant, wpx, hpx):
        """ 
        Returns a dictionary with a tiling [coordinates of corners of tiles]
        of quadrant Q, with tiles of size wpx[width] x hpx[height].
           CAUTION: Returned coordinates are Q-relative.

        :param Quadrant: str, Quadrant, one of ['E','F','G','H']
        :param wpx: int, width [along NAXIS1] of tiles, in pixels.
        :param hpx: int, height [along NAXIS2] of tiles, in pixels.
        :return: tiles_dict = dict(
                          wpx='Width of tiles, integer', 
                          hpx='Height of tiles, integer', 
                          llpix='Lower left corner of tiles, list of tuples',
                          ccpix= 'Central pixel of tiles, list of tuples', 
                          Nsamps='Number of tiles, integer')

        """

        tiles_dict = dict()

        # Retrieve boundaries (cols and rows) of image sections: s-prescan, image, s-overscan, p-overscan
        _prestarth, _preendh, _imgstarth, _imgendh, _ovstarth, _ovendh = self.getsectioncollims(
            Quadrant)

        _imgstartv, _imgendv, _ovstartv, _ovendv = self.getsectionrowlims(
            Quadrant)

        # IMAGE area boundaries

        imglims = [_imgstarth, _imgendh,
                   _imgstartv, _imgendv]

        # sampling points (lower-left)

        xsamp = np.arange(imglims[0], imglims[1]-wpx-1, step=wpx)
        ysamp = np.arange(imglims[2], imglims[3]-hpx-1, step=hpx)

        # lower-left pixels of tiles
        llpix = list(itertools.product(xsamp, ysamp))
        # Number of samples-tiles
        Nsamps = len(llpix)
        # centre-pixels of tiles
        ccpix = [(llpix[i][0]+wpx/2., llpix[i][1]+hpx/2.)
                 for i in range(Nsamps)]
        # dictionary with results
        tiles_dict = dict(wpx=wpx, hpx=hpx, llpix=llpix,
                          ccpix=ccpix, Nsamps=Nsamps)

        return tiles_dict

    def get_tiles(self, Quadrant, tile_coos, extension=-1):
        """ """

        tiles = []
        Nsamps = tile_coos['Nsamps']
        wpx = tile_coos['wpx']
        hpx = tile_coos['hpx']
        for i in range(Nsamps):
            llx = tile_coos['llpix'][i][0]
            lly = tile_coos['llpix'][i][1]
            cc = [llx, llx+wpx+1, lly, lly+hpx+1]

            tiles.append(self.get_cutout(
                cc, Quadrant, canonical=False, extension=extension))
        return tiles

    def get_tiles_stats(self, Quad, tile_coos, statkey, extension=-1):
        """ """

        if isinstance(self.extensions[extension].data, np.ma.masked_array):
            stat_dict = dict(
                mean=np.ma.mean, median=np.ma.median, std=np.ma.std)
        else:
            stat_dict = dict(mean=np.mean, median=np.median, std=np.std)

        estimator = stat_dict[statkey]

        tiles = self.get_tiles(Quad, tile_coos, extension=extension)

        vals = np.array(map(estimator, tiles))

        return vals

    def get_cutout(self, corners, Quadrant, canonical=False, extension=-1):
        """Returns a cutout from the CCD image, either in 
        canonical or non-canonical orientation.


        :param corners: [x0,x1,y0,y1]
        :type corners: list (of int)

        :param Quadrant: Quadrant, one of 'E', 'F', 'G', 'H'
        :type Quadrant: char

        :param canonical: 
         Canonical [True] = with readout-node at pixel index (0,0) regardless of quadrant. 
         This is the orientation which corresponds to the data-readin order 
         (useful for cross-talk measurements, for example).
         Non-Canonical [False] = with readout-node at corner matching placement of  
         quadrant on the CCD. This is the orientation that would match the representation 
         of the image on DS9.        
        :type canonical: bool

        :param extension: extension number. Default = -1 (last)
        :type extension: int


        """
        Qdata = self.get_quad(Quadrant, canonical, extension)
        cutout = Qdata[corners[0]:corners[1], corners[2]:corners[3]].copy()
        return cutout

    def set_quad(self, inQdata, Quadrant, canonical=False, extension=-1):
        """ """
        edges = self.QuadBound[Quadrant]
        #Qdata = self.data[edges[0]:edges[1],edges[2]:edges[3]]

        if canonical:
            if Quadrant == 'E':
                Qdata = inQdata[:, ::-1].copy()
            elif Quadrant == 'F':
                Qdata = inQdata[::-1, ::-1].copy()
            elif Quadrant == 'G':
                Qdata = inQdata[::-1, :].copy()
            elif Quadrant == 'H':
                Qdata = inQdata[:, :].copy()
        else:
            Qdata = inQdata.copy()

        self.extensions[extension].data[edges[0]:edges[1], edges[2]:edges[3]] = Qdata.copy()

        return None

    def getsectioncollims(self, Q):
        """Returns limits of [HORIZONTAL] sections: prescan, image and overscan"""

        if Q in ['E', 'H']:

            prestart = 0
            preend = self.prescan-1
            ovstart = self.wQ-self.overscan
            ovend = ovstart + self.overscan - 1
            imgstart = preend+1
            imgend = ovstart-1

        elif Q in ['F', 'G']:

            ovstart = 0
            ovend = self.overscan-1
            prestart = self.wQ-self.prescan
            preend = prestart + self.prescan-1
            imgstart = ovend+1
            imgend = prestart-1

        return (prestart, preend, imgstart, imgend, ovstart, ovend)

    def getsectionrowlims(self, Q):
        """Returns limits of [VERTICAL] sections: image [and vertical overscan]"""

        vover = self.voverscan
        img = self.NrowsCCD

        if Q in ['E', 'F']:

            imgstart = vover
            imgend = imgstart + img - 1

            if self.voverscan != 0:
                ovstart = 0
                ovend = ovstart + vover - 1
            else:
                ovstart = None
                ovend = None

        elif Q in ['G', 'H']:

            imgstart = 0
            imgend = img - 1

            if self.voverscan != 0:
                ovstart = img
                ovend = img + vover - 1
            else:
                ovstart = None
                ovend = None

        return (imgstart, imgend, ovstart, ovend)

    def get_stats(self, Quadrant, sector='img', statkeys=['mean'], trimscan=[0, 0],
                  ignore_pover=True, extension=-1, VSTART=0, VEND=NrowsCCD+voverscan):
        """ """

        Qdata = self.get_quad(Quadrant, canonical=True, extension=extension)

        if isinstance(Qdata, np.ma.masked_array):
            stat_dict = dict(
                mean=np.ma.mean, median=np.ma.median, std=np.ma.std)
        else:
            stat_dict = dict(mean=np.mean, median=np.median, std=np.std)

        if ignore_pover:
            vlims = [VSTART, min(self.NrowsCCD, VEND)]
        else:
            vlims = [VSTART, VEND]

        if sector == 'pre':
            hlims = [0, self.prescan]
        elif sector == 'ove':
            hlims = [self.NAXIS1/2-self.overscan, self.NAXIS1/2]
        elif sector == 'img':
            hlims = [self.prescan, self.NAXIS1/2-self.overscan]
        elif sector == 'all':
            hlims = [0, None]

        hlims[0] += trimscan[0]
        hlims[1] -= trimscan[1]

        results = []

        for statkey in statkeys:
            results.append(stat_dict[statkey](
                Qdata[hlims[0]:hlims[1], vlims[0]:vlims[1]]))

        return results

    def sub_offset(self, Quad, method='row', scan='pre', trimscan=[3, 2],
                   ignore_pover=True, extension=-1):
        """ """

        if self.masked:
            median = np.ma.median
        else:
            median = np.median

        quaddata = self.get_quad(Quad, canonical=True, extension=extension)

        if scan == 'pre':
            hlims = [0, self.prescan]
        elif scan == 'ove':
            hlims = [self.NAXIS1/2-self.overscan, self.NAXIS1/2]
        else:
            sys.exit('ccd.sub_offset: scan=%s unkonwn' % scan)

        hlims[0] += trimscan[0]
        hlims[1] -= trimscan[1]

        if ignore_pover:
            vlims = [0, self.NrowsCCD]
        else:
            vlims = [0, None]

        if method == 'row':
            offsets = []
            for ix in range(self.NAXIS2/2):
                offset = median(quaddata[hlims[0]:hlims[1], ix])
                quaddata[:, ix] -= offset
                offsets.append(offset)

        elif method == 'median':

            offset = median(quaddata[hlims[0]:hlims[1], vlims[0]:vlims[1]])
            #if self.masked : offset = offset.data
            quaddata -= offset
            offsets = [offset]

        B = self.QuadBound[Quad]
        self.extensions[extension].data[B[0]:B[1], B[2]:B[3]
                                        ] = self.flip_tocanonical(quaddata, Quad).copy()

        params = dict(Quad=Quad, method=method, scan=scan, trimscan=trimscan,
                      ignore_pover=ignore_pover, extension=extension)
        self.add_to_hist('sub_offset', extension, params=params)

        return offsets

    def sub_bias(self, superbias, extension=-1):
        """Subtracts a superbias"""

        assert self.shape == superbias.shape
        self.extensions[extension].data -= superbias
        self.add_to_hist('sub_bias', extension, vison=__version__,
                         params=dict(superbias=superbias))

    def divide_by_flatfield(self, FF, extension=-1):
        """Divides by a Flat-field"""
        print 'TODO: ccd.CCD.divide_by_flatfield needs improvements: handling of masked values, overscans (V & H)'

        assert self.shape == FF.shape
        self.extensions[extension].data /= FF

        self.add_to_hist('divide_by_flatfield', extension, vison=__version__,
                         params=dict(FF=FF))

    def add_to_hist(self, action, extension=-1, vison=__version__, params=dict()):
        """ """
        tstamp = (datetime.datetime.now()).strftime('%d%m%yD%H%M%ST')
        hist = [dict(timestamp=tstamp, action=action, extension=extension,
                     vison=vison, params=params)]
        self.historial += hist

    def flip_tocanonical(self, array, Quad):

        if Quad == 'E':
            return array[:, ::-1].copy()
        elif Quad == 'F':
            return array[::-1, ::-1].copy()
        elif Quad == 'G':
            return array[::-1, :].copy()
        elif Quad == 'H':
            return array.copy()

    def do_Vscan_Mask(self, VSTART, VEND):
        """ """

        VscanMask = np.ones((self.NAXIS1, self.NAXIS2), dtype='bool')

        for Quad in self.QuadBound.keys():

            B = self.QuadBound[Quad]

            tmp = self.flip_tocanonical(VscanMask[B[0]:B[1], B[2]:B[3]], Quad)
            tmp[:, VSTART:VEND+1] = False
            VscanMask[B[0]:B[1], B[2]:B[3]] = self.flip_tocanonical(
                tmp, Quad).copy()

        return VscanMask

    def get_mask(self, mask):
        """ """
        assert self.shape == mask.shape

        for iext in range(self.nextensions):
            if self.extensions[iext].data is not None:
                masked = np.ma.masked_array(self.extensions[iext].data, mask)
                self.extensions[iext].data = masked.copy()

        self.masked = True

    def writeto(self, fitsf, clobber=False, unsigned16bit=False):
        """ """

        #prihdu = fts.PrimaryHDU()

        firstextension = self.extensions[0]

        if firstextension.data is not None:
            if self.masked:
                pridata = firstextension.data.data.transpose().copy()
            else:
                pridata = firstextension.data.transpose().copy()
        else:
            pridata = None

        prihdr = firstextension.header

        prihdu = fts.PrimaryHDU(data=pridata, header=prihdr)

        comments = ['FITS file generated by vison.datamodel.ccd',
                    'point of contact: Ruyman Azzollini (r.azzollini_at_ucl.ac.uk)', ]
        for comm in comments:
            prihdu.header.add_comment(comm)

        hdulist = fts.HDUList([prihdu])

        if self.nextensions > 1:

            for iext in range(1, self.nextensions):

                if self.masked:
                    idata = self.extensions[iext].data.data.T.copy()
                else:
                    idata = self.extensions[iext].data.T.copy()

                iheader = self.extensions[iext].header
                iname = self.extensions[iext].label

                ihdu = fts.ImageHDU(data=idata, header=iheader, name=iname)

                if unsigned16bit:
                    ihdu.scale('int16', '', bzero=32768)
                    ihdu.header.add_history(
                        'Scaled to unsigned 16bit integer!')

                hdulist.append(ihdu)

        hdulist.writeto(fitsf, overwrite=clobber)


class CCDPile(CCD):
    """Class to hold and operate (e.g. stack) on a bunch of CCD images.
    Each image (a single extension picked from each) becomes an extension in the pile.
    """

    def __init__(self, infitsList=[], ccdobjList=[], extension=-1, withpover=True):
        """ """
        
        self.masked = False
        self.extensions = []

        self.NAXIS1 = NAXIS1
        self.withpover = withpover
        if self.withpover:
            self.NAXIS2 = NAXIS2
        else:
            self.NAXIS2 = NAXIS2-40

        self.shape = (self.NAXIS1, self.NAXIS2)

        if len(infitsList) > 0:

            for i, infits in enumerate(infitsList):
                
                infits = str(infits)

                assert type(infits) is str,\
                                "'%s' can't be a name for a file!" % infits

                assert isthere(
                    infits), 'infits:%s is just not there :-(' % infits

                iccdobj = CCD(infits, extensions=[extension], withpover=self.withpover,
                              getallextensions=False)

                assert self.shape == iccdobj.shape

                self.extensions.append(iccdobj.extensions[extension])

        elif len(ccdobjList) > 0:
            

            for i, iccdobj in enumerate(ccdobjList):

                assert self.shape == iccdobj.shape

                self.extensions.append(iccdobj.extensions[extension])
                
            self.masked = np.any([isinstance(item.data,np.ma.core.MaskedArray) for item\
                                in self.extensions])

        else:

            self.extensions = []

        self.nextensions = len(self.extensions)

        self.prescan = prescan
        self.overscan = overscan
        if self.withpover:
            self.voverscan = voverscan
        else:
            self.voverscan = 0

        self.gain = dict(E=3.1, F=3.1, G=3.1, H=3.1)
        self.rn = dict(E=4.5, F=4.5, G=4.5, H=4.5)

        self.QuadBound = QuadBound

        

    def stack(self, method='median', dostd=False):
        """ """
        
        fstack_dict = dict(
                median=dict(
                        masked=np.ma.median,
                        unmasked=np.median),
                mean=dict(
                        masked=np.ma.mean,
                        unmasked=np.mean))
        if self.masked:
            fstack = fstack_dict[method]['masked']
            fstd = np.ma.std
            stackimg = np.ma.zeros(self.shape,dtype='float32')
            if dostd:
                stackstd = np.ma.zeros(stackimg.shape, dtype='float32')
            blankarray = np.ma.core.zeros
        else:
            fstack = fstack_dict[method]['unmasked']
            fstd = np.std
            stackimg = np.zeros(self.shape,dtype='float32')
            if dostd:
                stackstd = np.zeros_like(stackimg,dtype='float32')
            blankarray = np.zeros
        

        NAXIS1 = self.NAXIS1
        NAXIS2 = self.NAXIS2
        
        for i in range(NAXIS2):
            
            imgrow = blankarray((NAXIS1, self.nextensions), dtype='float32')
            
            for j in range(self.nextensions):
                imgrow[:, j] = self.extensions[j].data[:, i]
            
            stackimg[:, i] = fstack(imgrow, axis=1).copy()

            if dostd:
                stackstd[:, i] = fstd(imgrow, axis=1).copy()
        
        if dostd:
            return stackimg, stackstd
        else:
            return stackimg


def test_create_from_scratch():
    """ """

    NAXIS1, NAXIS2 = 4238, 4172

    img = np.ones((NAXIS1, NAXIS2), dtype='float32')
    eimg = np.ones((NAXIS1, NAXIS2), dtype='float32') * 0.1

    ccdout = CCD()

    fitsname = 'test_create_from_scratch.fits'

    ccdout.add_extension(data=img, label='IMAGE')
    ccdout.add_extension(data=eimg, label='UNCERTAINTY')
    ccdout.writeto(fitsname, clobber=True)

    ccdin = CCD(fitsname, getallextensions=True)

    print 'Number of extensions = %i' % ccdin.nextensions
    stop()


def test_load_ELVIS_fits():
    """ """

    fitsname = '/home/raf/WORK/EUCLID/REPOS/vison/vison/data/EUC_2112_231016D_135042T_ROE1_CCD1.fits'

    ccd = CCD(infits=fitsname, getallextensions=True)

    ccd.writeto('ccd_test_load_ELVIS_fits.fits', clobber=True)


if __name__ == '__main__':

    # test_create_from_scratch()
    test_load_ELVIS_fits()
