# -*- coding: utf-8 -*-
"""

data model for CCDs

Created on Fri Nov 13 17:42:36 2015

@author: raf
"""

# IMPORT STUFF
from astropy.io import fits as fts
import numpy as np
import os
from pdb import set_trace as stop
import sys
# END IMPORT

isthere = os.path.exists

#class Quadrant():
#    
#    def __init__(self,data):
#                
#        self.data = data.copy()
#        shape = self.data.shape
#        
#        self.NAXIS1 = 2119
#        self.NAXIS2 = 2066
#        
#        assert self.NAXIS1 == shape[0]
#        assert self.NAXIS2 == shape[1]
#        
#        self.prescan = 51
#        self.overscan = 20
#        self.gain = 3.1
#        self.readnoise = 4.5

NAXIS1 = 4238
NAXIS2 = 4132
prescan = 51
overscan = 20
imgarea = [2048,2066]

QuadBound = dict(E=[0,NAXIS1/2,NAXIS2/2,NAXIS2],\
               F=[NAXIS1/2,NAXIS1,NAXIS2/2,NAXIS2],\
               G=[NAXIS1/2,NAXIS1,0,NAXIS2/2],\
               H=[0,NAXIS1/2,0,NAXIS2/2])

class CCD():
    
    def __init__(self,infits=None,extension=-1):
        """ """
        
        if infits is not None:
        
            assert type(infits) is str, 'infits is not a name of a file'
            assert isthere(infits), 'infits is just not there'
            
            f = fts.open(infits)
            
            #assert len(f) == 2, 'len of fits object should be 2'
            
            self.data = f[extension].data.transpose() # Beware!!
            self.data = self.data.astype('float32').copy()
            self.header = f[extension].header
        else:
            self.data = []
            self.header = ''
        
        self.NAXIS1 = 4238
        self.NAXIS2 = 4132
        
        if infits is not None:
            assert self.NAXIS1 == self.data.shape[0]
            assert self.NAXIS2 == self.data.shape[1]
        
        self.prescan = prescan
        self.overscan = overscan
        self.gain = dict(E=3.1,F=3.1,G=3.1,H=3.1)
        self.rn = dict(E=4.5,F=4.5,G=4.5,H=4.5)
        
        self.QuadBound = QuadBound 
        
        self.masked = False
        
    def get_quad(self,Quadrant,canonical=False):
        """ """
        
        edges = self.QuadBound[Quadrant]        
        Qdata = self.data[edges[0]:edges[1],edges[2]:edges[3]]
        
        if canonical:
            if Quadrant == 'E': Qdata = Qdata[:,::-1].copy()
            elif Quadrant == 'F': Qdata = Qdata[::-1,::-1].copy()
            elif Quadrant == 'G': Qdata = Qdata[::-1,:].copy()
            elif Quadrant == 'H': Qdata = Qdata[:,:].copy()
        
        
        return Qdata

    def set_quad(self,inQdata,Quadrant,canonical=False):
        """ """
        
        edges = self.QuadBound[Quadrant]        
        Qdata = self.data[edges[0]:edges[1],edges[2]:edges[3]]
        
        if canonical:
            if Quadrant == 'E': Qdata = inQdata[:,::-1].copy()
            elif Quadrant == 'F': Qdata = inQdata[::-1,::-1].copy()
            elif Quadrant == 'G': Qdata = inQdata[::-1,:].copy()
            elif Quadrant == 'H': Qdata = inQdata[:,:].copy()
        else:
            Qdata = inQdata.copy()
        
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


   
    def meas_bias(self,Quadrant,sector='both',detail=False):
        """ """
        
        Qdata = self.get_quad(Quadrant,canonical=True)
        
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
                
    
    def sub_offset(self,Quad,method='row',scan='pre',trimscan=[3,2]):
        """ """
        
        if self.masked:
            median= np.ma.median
        else:
            median = np.median
        
        quaddata = self.get_quad(Quad,canonical=True)
        
        #if Quad == 'E': quaddata = quaddata[:,::-1].copy()
        #elif Quad == 'F': quaddata = quaddata[::-1,::-1].copy()
        #elif Quad == 'G': quaddata = quaddata[::-1,:].copy()
        #elif Quad == 'H': quaddata = quaddata.copy()
        
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
        self.data[B[0]:B[1],B[2]:B[3]] = self.flip_tocanonical(quaddata,Quad).copy()
        
        #if Quad == 'E': self.data[B[0]:B[1],B[2]:B[3]] = quaddata[:,::-1].copy()
        #elif Quad == 'F': self.data[B[0]:B[1],B[2]:B[3]] = quaddata[::-1,::-1].copy()
        #elif Quad == 'G': self.data[B[0]:B[1],B[2]:B[3]] = quaddata[::-1,:].copy()
        #elif Quad == 'H': self.data[B[0]:B[1],B[2]:B[3]] = quaddata.copy()
        
        
        return offsets
        
    
    def sub_bias(self,superbias):
        """Subtracts a superbias"""
        
        assert self.data.shape == superbias.shape
        self.data -= superbias
        
    
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
        assert self.data.shape == mask.shape
        masked = np.ma.masked_array(self.data,mask)
        self.data = masked.copy()
        #self.mask = mask.copy()
        self.masked = True
        
    
    def writeto(self,fitsf,clobber=False):
        """ """
        
        if self.masked: data = self.data.data.copy()
        else: data = self.data.copy()
        
        hdu = fts.PrimaryHDU(self.data.transpose())
    
        #hdu.header['OBSID']=OBSID
        #hdu.header['CCD']=CCD   
        comments = ['FITS generated by vissim.datamodel.ccd',
        'point of contact: Ruyman Azzollini (r.azzollini@ucl.ac.uk)',]
        for comm in comments:
            hdu.header.add_comment(comm)
    
        hdulist = fts.HDUList([hdu])
        hdulist.writeto(fitsf,clobber=clobber)
        