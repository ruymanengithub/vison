# -*- coding: utf-8 -*-
"""

Data model for Euclid-VIS CCDs (ground testing at MSSL)

Created on Fri Nov 13 17:42:36 2015

:Author: Ruyman Azzollini
"""

# IMPORT STUFF
from astropy.io import fits as fts
import numpy as np
import os
from pdb import set_trace as stop
import sys
# END IMPORT

isthere = os.path.exists

NAXIS1 = 4238
NAXIS2 = 4132
prescan = 51
overscan = 20
imgarea = [2048,2066]

QuadBound = dict(E=[0,NAXIS1/2,NAXIS2/2,NAXIS2],\
               F=[NAXIS1/2,NAXIS1,NAXIS2/2,NAXIS2],\
               G=[NAXIS1/2,NAXIS1,0,NAXIS2/2],\
               H=[0,NAXIS1/2,0,NAXIS2/2])


class Extension():
    """Extension Class"""
    
    def __init__(self,data,header=None,label=None,headerdict=None):
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



class CCD(object):
    """Class of CCD objects. 
    Euclid Images as acquired by ELVIS software (Euclid LabView Imaging Software).
    
    
    The class has been extended to handle multi-extension images. This is useful
    to also "host" calibration data-products, such as Flat-Fields.
    
    """
    
    def __init__(self,infits=None,extensions=[-1],getallextensions=False):
        """ """
        
        self.extnames = []
        self.extensions = []
        
        
        if infits is not None:
        
            assert type(infits) is str, "infits can't be a name for a file!"
            assert isthere(infits), 'infits is just not there :-('
            
            
            self.loadfromFITS(infits,extensions,getallextensions)

        else:
            
            self.extensions = []
            self.extnames = []
        
        self.nextensions = len(self.extensions)
        
        self.NAXIS1 = NAXIS1
        self.NAXIS2 = NAXIS2        
        self.shape = (NAXIS1,NAXIS2)
        
        for iext in range(self.nextensions):
            if self.extensions[iext].data is not None:
                assert self.shape == self.extensions[iext].data.shape                
        
        self.prescan = prescan
        self.overscan = overscan
        self.gain = dict(E=3.1,F=3.1,G=3.1,H=3.1)
        self.rn = dict(E=4.5,F=4.5,G=4.5,H=4.5)
        
        self.QuadBound = QuadBound 
        
        self.masked = False
    
    
    def loadfromFITS(self,fitsfile,extensions=[-1],getallextensions=False):
        
        hdulist = fts.open(fitsfile)
        
        
        nextensions = len(hdulist)
            
        if getallextensions:
            extensions = np.arange(nextensions)
            
        for iext in extensions:
                
            hdu = hdulist[iext]
            
            
            if hdu.data is not None:
                data = hdu.data.transpose().astype('float32').copy()
            else: data = None
            header = hdu.header
                
            if 'EXTNAME' in hdu.header:
                label = hdu.header['EXTNAME']
            else:
                label = None
                
            self.extensions.append(Extension(data,header,label))
            self.extnames.append(label)
            
        hdulist.close()

    
    def add_extension(self,data,header=None,label=None,headerdict=None):
        """ """    
        self.extensions.append(Extension(data,header,label,headerdict))
        self.nextensions += 1
        
    
    def get_quad(self,Quadrant,canonical=False,extension=-1):
        """Returns a quadrant in canonical or non-canonical orientation.
        
        :param Quadrant: Quadrant, one of 'E', 'F', 'G', 'H'
        :type Quadrant: char
        
        :param canonical: 
        
        Canonical [True] = with readout-node at pixel index (0,0) regardless of quadrant. This is the orientation which corresponds to the data-reading order (useful for cross-talk measurements, for example).
        Non-Canonical [False] = with readout-node at corner matching placement of quadrant on the CCD. This is the orientation that would match the representation of the image on DS9.        

        :type canonical: bool
        
        :param extension: extension number. Default = -1 (last)
        :type extension: int
        
        """
        
        edges = self.QuadBound[Quadrant]        
        Qdata = self.extensions[extension].data[edges[0]:edges[1],edges[2]:edges[3]]
        
        if canonical:
            if Quadrant == 'E': Qdata = Qdata[:,::-1].copy()
            elif Quadrant == 'F': Qdata = Qdata[::-1,::-1].copy()
            elif Quadrant == 'G': Qdata = Qdata[::-1,:].copy()
            elif Quadrant == 'H': Qdata = Qdata[:,:].copy()
        
        
        return Qdata
        
    def get_cutout(self,corners,Quadrant,canonical=False,extension=-1):
        """Returns a cutout from the CCD image, either in 
        canonical or non-canonical orientation.
        
        
        :param corners: [x0,x1,y0,y1]
        :type corners: list (of int)
        
        :param Quadrant: Quadrant, one of 'E', 'F', 'G', 'H'
        :type Quadrant: char
        
        :param canonical: 
        
         Canonical [True] = with readout-node at pixel index (0,0) regardless of quadrant. This is the orientation which corresponds to the data-readin order (useful for cross-talk measurements, for example).
         Non-Canonical [False] = with readout-node at corner matching placement of  quadrant on the CCD. This is the orientation that would match the representation of the image on DS9.        
                      
        :type canonical: bool
        
        :param extension: extension number. Default = -1 (last)
        :type extension: int
        
        
        """       
        Qdata = self.get_quad(Quadrant,canonical,extension)        
        section = Qdata[corners[0]:corners[1],corners[2]:corners[3]].copy()
        return section
        
        
    def set_quad(self,inQdata,Quadrant,canonical=False,extension=-1):
        """ """
        
        edges = self.QuadBound[Quadrant]        
        Qdata = self.extensions[extension].data[edges[0]:edges[1],edges[2]:edges[3]]
        
        if canonical:
            if Quadrant == 'E': Qdata = inQdata[:,::-1].copy()
            elif Quadrant == 'F': Qdata = inQdata[::-1,::-1].copy()
            elif Quadrant == 'G': Qdata = inQdata[::-1,:].copy()
            elif Quadrant == 'H': Qdata = inQdata[:,:].copy()
        else:
            Qdata = inQdata.copy()
        
        self.extensions[extension].data = Qdata.copy()
        
        return None


    def getsectioncollims(self,QUAD):
        """Returns limits of sections: prescan, image and overscan"""
        
        semiNAXIS1 = self.NAXIS1/2
        
        if QUAD in ['E','H']:
                
            prestart = 0
            preend = self.prescan-1
            ovstart = semiNAXIS1-self.overscan 
            ovend = ovstart + self.overscan - 1
            imgstart = preend+1
            imgend = ovstart-1
            
        elif QUAD in ['F','G']:
            
            ovstart = 0
            ovend = self.overscan-1
            prestart = semiNAXIS1-self.prescan
            preend = prestart + self.prescan-1
            imgstart = ovend+1
            imgend = prestart-1
                
        return (prestart,preend,imgstart,imgend,ovstart,ovend)


   
    def meas_bias(self,Quadrant,sector='both',detail=False,extension=-1):
        """ """
        
        Qdata = self.get_quad(Quadrant,canonical=True,extension=extension)
        
        prescan = Qdata[4:self.prescan-1,3:-3]
        ovscan = Qdata[self.prescan+2048+4:-3,3:-3]
        
        #stop()
        
        if sector == 'both':
            virpixels = np.concatenate((prescan,ovscan))
        elif sector == 'pre':
            virpixels = prescan.copy()
        elif sector == 'over':
            virpixels = ovscan.copy()
        
        bias = np.mean(virpixels)
        stdbias = np.std(virpixels)
        
                
        if not detail:
            return bias,stdbias
        else:
            biasodd = np.mean(virpixels[0::2,:])
            stdodd = np.std(virpixels[0::2,:])
            biaseven = np.mean(virpixels[1::2,:])
            stdeven = np.std(virpixels[1::2,:])
            #if Quadrant=='F': stop()
            return bias,stdbias,biasodd,stdodd,biaseven,stdeven
                
    
    def sub_offset(self,Quad,method='row',scan='pre',trimscan=[3,2],extension=-1):
        """ """
        
        if self.masked:
            median= np.ma.median
        else:
            median = np.median
        
        quaddata = self.get_quad(Quad,canonical=True,extension=extension)
        
        
        if scan == 'pre':
            lims = [0,self.prescan]
        elif scan == 'ove':
            lims = [self.NAXIS1/2-self.overscan,self.NAXIS1/2]
        else:
            sys.exit('ccd.sub_offset: scan=%s unkonwn' % scan)
            
        lims[0] += trimscan[0]
        lims[1] -= trimscan[1]

        if method == 'row':
            offsets = []
            for ix in range(self.NAXIS2/2):
                offset = median(quaddata[lims[0]:lims[1],ix])
                if self.masked : offset = offset.data
                quaddata[:,ix] -= offset
                offsets.append(offset)
            
        elif method == 'median':
            
            offset = median(quaddata[lims[0]:lims[1],:])
            #if self.masked : offset = offset.data
            quaddata -= offset
            offsets = [offset]
        
        B = self.QuadBound[Quad]
        self.extension[extension].data[B[0]:B[1],B[2]:B[3]] = self.flip_tocanonical(quaddata,Quad).copy()
        

        return offsets
        
    
    def sub_bias(self,superbias,extension=-1):
        """Subtracts a superbias"""
        
        assert self.shape == superbias.shape
        self.extensions[extension].data -= superbias
        
    
    def flip_tocanonical(self,array,Quad):
        
        if Quad == 'E': return array[:,::-1].copy()
        elif Quad == 'F': return array[::-1,::-1].copy()
        elif Quad == 'G': return array[::-1,:].copy()
        elif Quad == 'H': return array.copy()
    
    def do_Vscan_Mask(self,VSTART,VEND):
        
        VscanMask = np.ones((self.NAXIS1,self.NAXIS2),dtype='bool')
        
        for Quad in self.QuadBound.keys():
            
            B = self.QuadBound[Quad]
            
            tmp = self.flip_tocanonical(VscanMask[B[0]:B[1],B[2]:B[3]],Quad)
            tmp[:,VSTART:VEND+1] = False
            VscanMask[B[0]:B[1],B[2]:B[3]] = self.flip_tocanonical(tmp,Quad).copy()
        
        
        return VscanMask
    
    def get_mask(self,mask):
        """ """
        assert self.shape == mask.shape
        
        for iext in range(self.nextensions):
            
            masked = np.ma.masked_array(self.extensions[iext].data,mask)
            self.extensions[iext].data = masked.copy()
        
        self.masked = True
        
    
    def writeto(self,fitsf,clobber=False):
        """ """
        
        #prihdu = fts.PrimaryHDU()
        
        firstextension = self.extensions[0]
          
        if firstextension.data is not None:
            if self.masked: pridata = firstextension.data.data.transpose().copy()
            else: pridata = firstextension.data.transpose().copy()
        else:
            pridata = None
            
        prihdr = firstextension.header
                
        prihdu = fts.PrimaryHDU(data=pridata,header=prihdr)

          
        comments = ['FITS file generated by vison.datamodel.ccd',
        'point of contact: Ruyman Azzollini (r.azzollini_at_ucl.ac.uk)',]
        for comm in comments:
            prihdu.header.add_comment(comm)
            
        
    
        hdulist = fts.HDUList([prihdu])
        
        if self.nextensions > 1:
            
            for iext in range(1,self.nextensions):
                
                if self.masked: idata = self.extensions[iext].data.data.transpose().copy()
                else: idata = self.extensions[iext].data.transpose().copy()
                
                iheader = self.extensions[iext].header
                iname = self.extensions[iext].label
                
                ihdu = fts.ImageHDU(data=idata,header=iheader,name=iname)
            
                hdulist.append(ihdu)
        
        hdulist.writeto(fitsf,overwrite=clobber)
        
def test_create_from_scratch():
    """ """
    
    NAXIS1,NAXIS2 = 4238,4132
    
    img = np.ones((NAXIS1,NAXIS2),dtype='float32')
    eimg = np.ones((NAXIS1,NAXIS2),dtype='float32') * 0.1
    
    ccdout = CCD()
    
    fitsname = 'test_create_from_scratch.fits'
    
    ccdout.add_extension(data=img,label='IMAGE')
    ccdout.add_extension(data=eimg,label='UNCERTAINTY')
    ccdout.writeto(fitsname,clobber=True)
    
    ccdin = CCD(fitsname,getallextensions=True)
    
    print 'Number of extensions = %i' % ccdin.nextensions
    stop()

def test_load_ELVIS_fits():
    """ """
    
    fitsname = '/home/raf/WORK/EUCLID/REPOS/vison/vison/data/EUC_2112_231016D_135042T_ROE1_CCD1.fits'
    
    ccd = CCD(infits=fitsname,getallextensions=True)
    
    ccd.writeto('ccd_test_load_ELVIS_fits.fits',clobber=True)
    
    
if __name__ == '__main__':
    
    #test_create_from_scratch()
    test_load_ELVIS_fits()
    