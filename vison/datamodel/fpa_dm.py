#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

FPA Data Model.
LE1 FITS files.

Created on Thu Aug  1 17:05:12 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from astropy.io import fits as fts
import os
import copy
import numpy as np
from collections import OrderedDict

from vison.fpa import fpa as fpamod
from vison.datamodel import ccd as ccdmod
# END IMPORT


Quads = ['E', 'F', 'G', 'H']
NSLICES = fpamod.NSLICES
NCOLS = fpamod.NCOLS

ext_ids = [None]
for j in range(1, NSLICES + 1):
    for i in [1, 2, 3]:

        ix = (j - 1) * 6 * 4 + (i - 1) * 4 + 1

        ext_ids.append((ix + 0, j, i, 'E', (0, 0)))
        ext_ids.append((ix + 1, j, i, 'F', (1, 0)))
        ext_ids.append((ix + 2, j, i, 'G', (1, 1)))
        ext_ids.append((ix + 3, j, i, 'H', (0, 1)))
    for i in [4, 5, 6]:

        ix = (j - 1) * 6 * 4 + (i - 1) * 4 + 1

        ext_ids.append((ix + 0, j, i, 'G', (0, 0)))
        ext_ids.append((ix + 1, j, i, 'H', (1, 0)))
        ext_ids.append((ix + 2, j, i, 'E', (1, 1)))
        ext_ids.append((ix + 3, j, i, 'F', (0, 1)))


PRIhdr_dict = OrderedDict()
PRIhdr_dict['BUNIT'] = 'ADU'
PRIhdr_dict['RA'] = None
PRIhdr_dict['DEC'] = None
PRIhdr_dict['DATE-OBS'] = None
PRIhdr_dict['SEQID'] = None

EXThdr_dict = OrderedDict()
EXThdr_dict['PCOUNT'] = 0
EXThdr_dict['PCOUNT'] = 1
EXThdr_dict['BZERO'] = 32768
EXThdr_dict['EXPTIME'] = 0
EXThdr_dict['CCDID'] = '0-0'
EXThdr_dict['QUADID'] = 'A'
EXThdr_dict['EXTNAME'] = '0-0.A'
EXThdr_dict['PRESCANX'] = 51
EXThdr_dict['OVRSCANX'] = 29
EXThdr_dict['GAIN'] = 3.5
EXThdr_dict['WCSAXES'] = 2
EXThdr_dict['INSTRUME'] = 'EUCLID-VIS'


class FPA_LE1(object):
    """Class for LE1 fits files built from system-level image data (whole FPA).

    This class relies on ccdobj instances to do analysis of the images (one CCD at a time).

    """

    Quads = ['E', 'F', 'G', 'H']

    def __init__(self, infits=None):
        """ """

        self.NEXTENSIONS = 145
        self.QNAXIS1 = 2128
        self.QNAXIS2 = 2086
        self.Qshape = (self.QNAXIS1, self.QNAXIS2)
        self.ext_ids = copy.deepcopy(ext_ids)
        self.extensions = []
        self.extnames = []
        if infits is not None:
            self.loadfromFITS(infits)
        self.fpamodel = fpamod.FPA()
        self.fillval = 0

    def add_extension(self, data, header, label=None, headerdict=None):
        """Adds an extension."""
        if data is not None:
            assert data.shape == self.Qshape

        self.extensions.append(ccdmod.Extension(data, header, label, headerdict))

    def set_extension(self, iext, data, header, label=None, headerdict=None):
        """Changes the contents of an extension."""
        if data is not None:
            assert data.shape == self.Qshape

        self.extensions[iext] = ccdmod.Extension(data, header, label, headerdict)


    def del_extension(self, ixextension):
        """Deletes an extension."""
        self.extensions.pop(ixextension)

    def initialise_as_blank(self, fillval=None):
        """Initialises object as a blank (filled with 'fillval') image of the FPA."""

        if fillval is None:
            fillval = self.fillval

        headerdict0 = PRIhdr_dict.copy()

        self.add_extension(data=None, header=None, label=None, headerdict=headerdict0)

        for iext in range(1, self.NEXTENSIONS):
            data = np.zeros((self.QNAXIS1, self.QNAXIS2), dtype='float32')+fillval
            header = None
            headerdict = EXThdr_dict.copy()

            _extid = self.ext_ids[iext]

            headerdict['CCDID'] = '%i-%i' % (_extid[1], _extid[2])
            headerdict['QUADID'] = _extid[3]
            headerdict['EXTNAME'] = '%s.%s' % (headerdict['CCDID'],headerdict['QUADID'])

            self.add_extension(data, header, label=None, headerdict=headerdict)

            self.extnames.append(None)

    def loadfromFITS(self, infits):
        """Loads contents of self from a FITS file."""

        hdulist = fts.open(infits)

        nextensions = len(hdulist)

        assert nextensions == self.NEXTENSIONS

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

            self.extensions.append(ccdmod.Extension(data, header, label))
            self.extnames.append(label)

        hdulist.close()

    def savetoFITS(self, outfits, clobber=True, unsigned16bit=False):
        """Dumps self to a FITS file."""

        prihdr = self.extensions[0].header

        prihdu = fts.PrimaryHDU(data=None,
                                header=prihdr)

        hdulist = fts.HDUList([prihdu])

        for iext in range(1, self.NEXTENSIONS):

            idata = self.extensions[iext].data.T.copy()

            iheader = self.extensions[iext].header
            iname = self.extensions[iext].label

            ihdu = fts.ImageHDU(data=idata, header=iheader, name=iname)

            if unsigned16bit:
                ihdu.scale('int16', '', bzero=32768)
                ihdu.header.add_history(
                    'Scaled to unsigned 16bit integer!')

            hdulist.append(ihdu)

        hdulist.writeto(outfits, overwrite=clobber)

    def get_CCDID_from_BLCCD(self, BLOCK, CCD):
        """

        Retrieves CCD ID (e.g. C_11) given BLOCK name and CCD in block.

        Parameters
        ----------

            :param BLOCK: block nickname (e.g. 'CURIE')
            :param CCD: 1, 2 or 3

        """

        return self.fpamodel.get_Ckey_from_BlockCCD(BLOCK, CCD)

    def get_extid(self, CCDID, Q):
        """
        Retrieves extension ID given CCDID and quadrant.
        
        Parameters
        ----------

            :param CCDID: e.g. 'C_11'
            :param Q: 'E', 'F', 'G' or 'H'

        """

        jC, iS = int(CCDID[2]), int(CCDID[3])

        for ix in range(1, self.NEXTENSIONS):
            _extid = self.ext_ids[ix]
            if _extid[1] == jC and _extid[2] == iS and _extid[3] == Q:
                return _extid

        return None

    def get_ccdobj(self, CCDID):
        """Returns a CCD Object given a CCDID"""

        ccdobj = ccdmod.CCD(withpover=True, overscan=29)

        blnk = np.zeros((ccdobj.NAXIS1, ccdobj.NAXIS2), dtype='float32')
        ccdobj.add_extension(blnk, header=None, label=None, headerdict=None)

        for Q in self.Quads:
            extid = self.get_extid(CCDID, Q)
            extix = extid[0]
            os_coo = extid[4]
            flip = os_coo

            Qdata = self.extensions[extix].data.copy()
            Qdata = self.fpamodel.flip_img(Qdata, flip)

            padQdata = np.zeros((ccdobj.NAXIS1 // 2, ccdobj.NAXIS2 // 2), dtype='float32')
            padQdata[:, :] = Qdata.copy()
            ccdobj.set_quad(padQdata, Q, canonical=True, extension=-1)

        return ccdobj

    def _padd_extra_soverscan(self, Qdata):
        """Padds extra serial overscan with self.fillval. 

        Used when building synthetic data of LE1 format using as input 
        images from the VGCC that only have 20 columns of serial 
        overscan (FPA images have 29)."""
        pQdata = np.zeros((self.QNAXIS1, self.QNAXIS2), dtype=Qdata.dtype)+self.fillval
        pQdata[0:Qdata.shape[0], 0:Qdata.shape[1]] = Qdata.copy()
        return pQdata

    def set_ccdobj(self, ccdobj, CCDID, inextension=-1):
        """Sets the image contents of input ccdobj to the image in self, at CCDID."""

        for Q in self.Quads:

            extid = self.get_extid(CCDID, Q)
            extix = extid[0]
            os_coo = extid[4]
            flip = os_coo

            Qdata = ccdobj.get_quad(Q, canonical=True, extension=inextension)

            if Qdata.shape[0] < self.QNAXIS1:
                Qdata = self._padd_extra_soverscan(Qdata)

            Qdata = self.fpamodel.flip_img(Qdata, flip)

            _extid = self.get_extid(CCDID,Q)

            extname = '%i-%i.%s' % (_extid[1],_extid[2],_extid[3])

            self.extensions[extix].data = Qdata.copy()
            self.extensions[extix].header['EXTNAME'] = extname


    def _core_funct_simul(self, ccdobj, CCDID=None, simputs=None):
        """ """
        raise NotImplementedError("child class implements abstract method")

    def simul(self, simputs=None, zerofirst=False):
        """Accessory function to simulate FPA images."""

        for jY in range(1, self.fpamodel.NSLICES + 1):
            for iX in range(1, self.fpamodel.NSLICES + 1):
                CCDID = 'C_%i%i' % (jY, iX)
                #print('Simulating CCD: %s' % CCDID)
                kccdobj = self.get_ccdobj(CCDID)

                if zerofirst:
                    kccdobj.extensions[-1].data *= 0

                kccdobj = self._core_funct_simul(kccdobj, CCDID, simputs)

                self.set_ccdobj(kccdobj, CCDID, inextension=-1)

    def apply_function_to_ccds(self, ccdfunction, **kwargs):
        """Applies a function to each CCD in self."""

        for jY in range(1, self.fpamodel.NSLICES + 1):
            for iX in range(1, self.fpamodel.NSLICES + 1):
                CCDID = 'C_%i%i' % (jY, iX)
                #print('Simulating CCD: %s' % CCDID)
                kccdobj = self.get_ccdobj(CCDID)

                kccdobj = ccdfunction(kccdobj, **kwargs)

                self.set_ccdobj(kccdobj, CCDID, inextension=-1)

    def get_as_FPAmosaic(self):
        """Returns a copy of self as an FPA mosaic"""

        QNAXIS1 = self.QNAXIS1
        QNAXIS2 = self.QNAXIS2
        NSLICES = self.fpamodel.NSLICES
        NCOLS = self.fpamodel.NCOLS

        img = np.zeros((QNAXIS1*2*NCOLS,QNAXIS2*2*NSLICES))

        for iext in range(1,self.NEXTENSIONS):
            
            extdata = self.extensions[iext].data.copy() # notice the y-inversion
            extname = self.extensions[iext].header['EXTNAME']
            islice = int(extname[0])
            icol = int(extname[2])
            Q = extname[4]

            llx = (icol-1)*2*QNAXIS1
            lly = (islice-1)*2*QNAXIS2

            if icol<=3:
                xdeltaQ = dict(E=0,F=1,G=1,H=0)
                ydeltaQ = dict(E=0,F=0,G=1,H=1)
            else:
                xdeltaQ = dict(E=1,F=0,G=0,H=1)
                ydeltaQ = dict(E=1,F=1,G=0,H=0)

            llx += xdeltaQ[Q]*QNAXIS1
            lly += ydeltaQ[Q]*QNAXIS2

            img[llx:llx+QNAXIS1,lly:lly+QNAXIS2] = extdata.copy()

        return img





def test1():
    """ """
    dpath = 'IWS_DATA_FORMATS/samples_SC456_AUG19/'

    LE1fits = os.path.join(
        dpath,
        'EUC_SIM_VIS-SCIENCE-52929-3_0509A4C85402-0115313_20181012T074449.852486Z_SC456-C7a_T2.fits')

    le1 = FPA_LE1(LE1fits)

    nextensions = le1.NEXTENSIONS

    ccd11 = le1.get_ccdobj('C_11')

    stop()


def test2():
    """ """

    class ThisFPA_LE1(FPA_LE1):

        def _core_funct_simul(self, ccdobj, CCDID=None, simputs=None):

            iCCDID  = int(CCDID[2:])
            
            for Q in self.Quads:

                Qdata = ccdobj.get_quad(Q,canonical=True)
                NX,NY = Qdata.shape
                pNY = NY-ccdobj.voverscan
                yvector = np.expand_dims(np.arange(pNY)/float(pNY) * 100.,0)
                Qdata[ccdobj.prescan:-ccdobj.overscan,0:-ccdobj.voverscan] = yvector.copy()
                Qdata[ccdobj.prescan:-ccdobj.overscan,0:-ccdobj.voverscan] += iCCDID * 500

                ccdobj.set_quad(Qdata,Q,canonical=True)

            return ccdobj

    le1 = ThisFPA_LE1()
    le1.initialise_as_blank(fillval=0)

    le1.simul(zerofirst=True)

    img = le1.get_as_FPAmosaic()

    fts.writeto('imgFPA.fits', img.T, overwrite=True)

    

if __name__ == '__main__':
    test2()
