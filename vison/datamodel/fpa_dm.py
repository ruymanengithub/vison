#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

FPA Data Model(s).

Created on Thu Aug  1 17:05:12 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from astropy.io import fits as fts
import os
import copy
import numpy as np

from vison.fpa import fpa as fpamod
from vison.datamodel import ccd as ccdmod
# END IMPORT


Quads = ['E','F','G','H']
NSLICES = 6 
NCOLS = 6

ext_ids = [None]
for j in range(1,NSLICES+1):    
    for i in [1, 2, 3]:
        
        ix = (j-1)*6*4 + (i-1)*4+1
        
        ext_ids.append((ix+0, j, i, 'E', (0,0)))
        ext_ids.append((ix+1, j, i, 'F', (1,0)))
        ext_ids.append((ix+2, j, i, 'H', (0,1)))
        ext_ids.append((ix+3, j, i, 'G', (1,1)))
    for i in [4,5,6]:
        
        ix = (j-1)*6*4 + (i-1)*4+1
        
        ext_ids.append((ix+0, j, i, 'G', (0,0)))
        ext_ids.append((ix+1, j, i, 'H', (1,0)))
        ext_ids.append((ix+2, j, i, 'F', (0,1)))
        ext_ids.append((ix+3, j, i, 'E', (1,1)))


class FPA_LE1(object):
    """ """
    
    Quads = ['E','F','G','H']
    
    def __init__(self,infits):
        """ """
        
        self.NEXTENSIONS = 145
        self.ext_ids = copy.deepcopy(ext_ids)
        self.extensions = []
        self.extnames = []
        self.loadfromFITS(infits)
        self.fpamodel = fpamod.FPA()
        
        
    def loadfromFITS(self,infits):
        """ """
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
    
    def savetoFITS(self,outfits, clobber=True, unsigned16bit=False):
        """ """
        
        prihdr = self.extensionsp[0].header

        prihdu = fts.PrimaryHDU(data=None, 
                                header=prihdr)
        
        hdulist = fts.HDUList([prihdu])
        
        
        for iext in range(1, self.nextensions):


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
        
            
    def get_CCDID(self, BLOCK, CCD):
        """
        BLOCK: block nickname (e.g. 'CURIE')
        CCDk: 'CCD1', 'CCD2' or 'CCD3'
        
        """
        
        return self.fpamodel.get_Ckey_from_BlockCCD(BLOCK,CCD)
        
    
    def get_extid(self,CCDID, Q):
        """ 
        CCDID: e.g. 'C_11'
        Q: 'E', 'F', 'G' or 'H'
        
        """
        
        jC, iS = int(CCDID[2]), int(CCDID[3])
        
        for ix in range(1,self.NEXTENSIONS):
            _extid = self.ext_ids[ix]
            if _extid[1]==jC and _extid[2]==iS and _extid[3]==Q:
                return _extid
        
        return None
        

    def get_ccdobj(self, CCDID):
        """Returns a CCD Object given a CCDID."""
        
                
        ccdobj = ccdmod.CCD(withpover=True)    
        blnk = np.zeros((ccdobj.NAXIS1,ccdobj.NAXIS2),dtype='float32')        
        ccdobj.add_extension(blnk, header=None, label=None,headerdict=None)
        
        for Q in self.Quads:
            extid = self.get_extid(CCDID, Q)
            extix = extid[0]
            os_coo = extid[4]
            flip = os_coo
            
            Qdata = self.extensions[extix].data.copy()            
            Qdata = self.fpamodel.flip_img(Qdata, flip)
            
            padQdata = np.zeros((ccdobj.NAXIS1/2,ccdobj.NAXIS2/2),dtype='float32')
            padQdata[:,0:ccdobj.NrowsCCD] = Qdata.copy()
            ccdobj.set_quad(padQdata, Q, canonical=True, extension=-1)
        
                
        return ccdobj
        


def test1():
    """ """
    dpath = 'IWS_DATA_FORMATS/samples_SC456_AUG19/'
    
    LE1fits = os.path.join(dpath,'EUC_SIM_VIS-SCIENCE-52929-3_0509A4C85402-0115313_20181012T074449.852486Z_SC456-C7a_T2.fits')
    
    le1 = FPA_LE1(LE1fits)
    
    
    nextensions = le1.NEXTENSIONS
    
    ccd11= le1.get_ccdobj('C_11')
    
    stop()
    
    
if __name__ == '__main__':
    test1()
