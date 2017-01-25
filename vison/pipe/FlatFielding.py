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

#from matplotlib import pyplot as plt

from vissim.datamodel import ccd as ccdmodule
#from vissim.support.files import cPickleDumpDictionary, cPickleRead
#from vissim.datamodel.HKtools import parseHKfiles
from vissim.simulator.simulator import __version__

from scipy import ndimage as nd
#from scipy import signal
from astropy.io import fits as fts
# END IMPORT

Quads = ['E','F','G','H']


def get_ilum_splines(img,filtsize=25,filtertype='median',Tests=False):
    """ """
    
    NX,NY = img.shape
    #xx,yy = np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
    
    if Tests:
        filtered = np.ones_like(img)
    else:
        
        if filtertype == 'median':    
            filtered = nd.median_filter(img,size=filtsize,mode='nearest') # ,cval=0)
        elif filtertype == 'mean':
            filtered = nd.uniform_filter(img,size=filtsize,mode='nearest') #,cval=0)
        #filtered = signal.medfilt(img,filtsize)
        #filtered = np.ones_like(img)
    
    
    nX = NX/filtsize
    nY = NY/filtsize
    
    samplebin = 10
    
    sx = np.arange(samplebin/2,nX*samplebin+samplebin/2,samplebin)
    sy = np.arange(samplebin/2,nY*samplebin+samplebin/2,samplebin)
    
    sxx,syy = np.meshgrid(sx,sy,indexing='ij')
    
    zz = filtered[(sxx,syy)]
    #zz = img[(xx,yy)]
    
    xx,yy = np.mgrid[0:NX:NX*1j,0:NY:NY*1j]
    
    
    pilum = interpolate.griddata((sxx.flatten(),syy.flatten()),zz.flatten(),
                                (xx,yy), method='cubic',fill_value=np.nan)
    
    nans = np.isnan(pilum)
    pilum[nans] = filtered[nans].copy()
    
    
    ILUM = dict(filtered=filtered,polyfit=pilum,polycoefs=[])
    
    
    return ILUM


def fit2D(xx,yy,zz,degree=1):
    """ """
    from astropy.modeling import models, fitting
    
    p_init = models.Polynomial2D(degree=degree)
    fit_p = fitting.LinearLSQFitter()
    
    with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter (if changing fitter...)
        warnings.simplefilter('ignore')
        p = fit_p(p_init, xx, yy, zz)
    
    return p

def get_ilum(img,pdegree=5,filtsize=15,filtertype='median',Tests=False):
    """ """
    #import warnings
    
    
    NX,NY = img.shape
    #xx,yy = np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
    
    
    if Tests:
        filtered = np.ones_like(img)
    else:
        if filtertype == 'median':    
            filtered = nd.median_filter(img,size=filtsize,mode='nearest') # 'constant',cval=0)
        elif filtertype == 'mean':
            filtered = nd.uniform_filter(img,size=filtsize,mode='nearest') # 'constant',cval=0)
        #filtered = signal.medfilt(img,filtsize)
        #filtered = np.ones_like(img)
    
    nX = NX/filtsize
    nY = NY/filtsize
    
    x = np.arange(filtsize/2,nX*filtsize+filtsize/2,filtsize)
    y = np.arange(filtsize/2,nY*filtsize+filtsize/2,filtsize)
    
    xx,yy = np.meshgrid(x,y,indexing='ij')
    
    zz = filtered[(xx,yy)]
    #zz = img[(xx,yy)]
    
    p = fit2D(xx,yy,zz,degree=pdegree)
    
    xp,yp= np.mgrid[:NX,:NY]
    pilum = p(xp, yp)
    
    ILUM = dict(filtered=filtered,polyfit=pilum,polycoefs=p)
    
    return ILUM


def _produce_SingleFlatfield(args):
    produce_SingleFlatfield(*args)


def produce_SingleFlatfield(infits,outfits,settings={},runonTests=False):
    """ """
    #runonTests = False
    
    insettings = dict(bmethod='row',bscan='pre',ilumtype='polynomial')
    insettings.update(settings)
    
    print infits
    
    #ffits = glob(os.path.join(datapath,'EUC_%i_*_CCD%i.fits' % (OBSID,CCD)))[0]       
    ccd = ccdmodule.CCD(infits)
    NX,NY = ccd.NAXIS1,ccd.NAXIS2
    
    # Cosmetic Masking
    
    if 'cosmetics' in insettings.keys():
        cosmeticsmask = insettings['cosmetics'].copy()
        ccd.get_mask(cosmeticsmask)
    
    # Master Bias
    
    if 'superbias' in insettings.keys():
        superbias = insettings['superbias'].copy()
        ccd.sub_bias(superbias)
    
    cube = np.zeros((NX,NY,3),dtype='float32')
    
    bmethod = insettings['bmethod']
    bscan = insettings['bscan']
    filtertype = insettings['filtertype']
    ilumtype = insettings['ilumtype']
    
    for Quad in Quads:
        
        B = ccdmodule.QuadBound[Quad]
        
        #print 'OBSID=%i, Q=%s' % (OBSID,Quad)
        
        ccd.sub_offset(Quad,method=bmethod,scan=bscan)
        
        img = ccd.get_quad(Quad,canonical=True).copy()
        
        stripedimg = img[ccdmodule.prescan:ccdmodule.imgarea[0]+ccdmodule.prescan,:]
        
        if ilumtype == 'polynomial':
            ILUM = get_ilum(stripedimg,pdegree=5,filtsize=15,filtertype=filtertype,Tests=runonTests)
        elif ilumtype == 'spline':
            ILUM = get_ilum_splines(stripedimg,filtsize=25,filtertype=filtertype,Tests=runonTests)
            
        #divbyfiltered = np.zeros_like(img)
        filtered = np.zeros_like(img)
        #divbypoly = np.zeros_like(img)
        poly = np.zeros_like(img)
        #divfiltbypoly = np.zeros_like(img)
        
        #divbyfiltered[51:2048+51,:] = stripedimg / ILUM['filtered']
        
        filtered[ccdmodule.prescan:ccdmodule.imgarea[0]+ccdmodule.prescan,:] = ILUM['filtered'].copy()

        #divbypoly[51:2048+51,:] = stripedimg / ILUM['polyfit']
        poly[ccdmodule.prescan:ccdmodule.imgarea[0]+ccdmodule.prescan,:] = ILUM['polyfit'].copy()
        #divfiltbypoly[51:2048+51,:] = ILUM['filtered']/ILUM['polyfit']
        
        
        if Quad == 'E': 
            cube[B[0]:B[1],B[2]:B[3],0] = img[:,::-1].copy()
            cube[B[0]:B[1],B[2]:B[3],1] = filtered[:,::-1].copy()
            cube[B[0]:B[1],B[2]:B[3],2] = poly[:,::-1].copy()
        elif Quad == 'F': 
            cube[B[0]:B[1],B[2]:B[3],0] = img[::-1,::-1].copy()
            cube[B[0]:B[1],B[2]:B[3],1] = filtered[::-1,::-1].copy()
            cube[B[0]:B[1],B[2]:B[3],2] = poly[::-1,::-1].copy()
        elif Quad == 'G': 
            cube[B[0]:B[1],B[2]:B[3],0] = img[::-1,:].copy()
            cube[B[0]:B[1],B[2]:B[3],1] = filtered[::-1,:].copy()
            cube[B[0]:B[1],B[2]:B[3],2] = poly[::-1,:].copy()
        elif Quad == 'H': 
            cube[B[0]:B[1],B[2]:B[3],0] = img.copy()
            cube[B[0]:B[1],B[2]:B[3],1] = filtered.copy()
            cube[B[0]:B[1],B[2]:B[3],2] = poly.copy()
    
    
    hdu = fts.PrimaryHDU()
    hdulist = fts.HDUList([hdu])
    
    hdulist.append(fts.ImageHDU(cube.transpose()))
    #hdulist.append(fts.ImageHDU(divbyfiltered.transpose()))
    #hdulist.append(fts.ImageHDU(filtered.transpose()))
    #hdulist.append(fts.ImageHDU(divbypoly.transpose()))
    #hdulist.append(fts.ImageHDU(poly.transpose()))
    #hdulist.append(fts.ImageHDU(divfiltbypoly.transpose()))
    hdulist.writeto(outfits,clobber=True)


def produce_IndivFlats(infits,outfits,settings,runonTests,processes=6):
    """ """
    
    assert len(infits) == len(outfits)
        
    arglist = []
    
    for ix in range(len(infits)): 
        arglist.append((infits[ix],outfits[ix],settings,runonTests))

    # TESTS
    
    #produce_SingleFlatfield(*arglist[0])
    #stop()
    
    #generate flats using multiprocessing
    pool = Pool(processes=processes)    
    #_produce_SingleFlatfield(arglist[0]) # TESTS    
    pool.map(_produce_SingleFlatfield, arglist)
    
    
def produce_MasterFlat(infits,outfits,mask=None,settings={}):
    """Produces a Master Flat out of a number of flat-illumination exposures.
    Takes the outputs from produce_IndivFlats."""
    
    insettings = dict(ilum='fit')
    insettings.update(settings)
    
    NAXIS1 = ccdmodule.NAXIS1/2
    NAXIS2 = ccdmodule.NAXIS2/2
    
    mflat = np.zeros((NAXIS1*2,NAXIS2*2))
    eflat = np.zeros((NAXIS1*2,NAXIS2*2)) 
    
    nin = len(infits)
    
    for Quad in Quads:
        
        B = ccdmodule.QuadBound[Quad]
        x0 = B[0]
        x1 = B[1]
        y0 = B[2]
        y1 = B[3]
    
        cubeflat = np.zeros((NAXIS1,NAXIS2,len(infits)))
        
        
        
        for ix in range(nin):
        
            print 'Loading img %i of %i..., Q=%s' % (ix+1,nin,Quad)
            inhdulist = fts.open(infits[ix])
        
            #for Quad in Quads:
            #    cubeflat = np.zeros((2119,2066,len(list_ffits)))
            img = inhdulist[1].data.transpose()[:,:,0].copy() # offset-subtracted exposure
            filtered = inhdulist[1].data.transpose()[:,:,1].copy() # medianfiltered bias-sub exposure
            polyfit = inhdulist[1].data.transpose()[:,:,2].copy() # polynomial fit to  bias-sub exposure
            #iflat = np.zeros_like(img)
            #iflat[51:2048+51,:] = img[51:2048+51,:]/filtered[51:2048+51,:]
        
            inhdulist.close()
            
            
            if insettings['ilum'] == 'fit': subdivimg = img[x0:x1,y0:y1] / polyfit[x0:x1,y0:y1]
            elif insettings['ilum'] == 'filtered': subdivimg = img[x0:x1,y0:y1] / filtered[x0:x1,y0:y1]    
            
            if mask is not None:
               subdivimg = np.ma.masked_array(subdivimg,mask[x0:x1,y0:y1])
               subdivimg[np.where(subdivimg.mask)] = 0
            
            subdivimg[np.where(np.isinf(subdivimg))] = 0
            
            cubeflat[:,:,ix] = subdivimg

            #cubeflat[:,:,iOBS] = np.zeros((NAXIS1*2,NAXIS2*2))
    
        
        qflat = np.zeros((NAXIS1,NAXIS2))
        eqflat = np.zeros((NAXIS1,NAXIS2))
    
        for ix in range(NAXIS1):
            qflat[ix,:] = np.nanmedian(cubeflat[ix,:,:],axis=1)
            eqflat[ix,:] = np.nanstd(cubeflat[ix,:,:],axis=1)/np.sqrt(nin)
        
        
        
        mflat[x0:x1,y0:y1] = qflat.copy()
        eflat[x0:x1,y0:y1] = eqflat.copy()
    
    
    hdu = fts.PrimaryHDU()
    
    if 'header' in insettings:
        for key in insettings['header'].keys():
            hdu.header[key] = insettings['header'][key]
    
    hdu.header['NCOMB'] = nin 
    hdu.header['ILUM'] = settings['ilum']
    hdu.header.add_history('Master FLAT')
    hdu.header.add_history('Extension 1: FLAT')
    hdu.header.add_history('Extension 2: eFLAT')
    hdu.header.add_history('CCD273 EM1A Characterization Campaign')
    hdu.header.add_history('Created by VISsim (version=%s) at %s' % (__version__, datetime.datetime.isoformat(datetime.datetime.now())))
    hdu.header.add_history('Further Info: Ruyman Azzollini (r.azzollini at mssl.ucl.ac.uk)')
    
    hdulist = fts.HDUList([hdu])
    
    mhdu = fts.ImageHDU(mflat.transpose())
    mhdu.header['EXTNAME'] = 'Flat'    
    hdulist.append(mhdu)
    ehdu = fts.ImageHDU(eflat.transpose())
    ehdu.header['EXTNAME'] = 'Uncertainty'
    hdulist.append(ehdu)
    #hdulist.append(fts.ImageHDU(filtmflat.transpose()))
        
    hdulist.writeto(outfits,clobber=True)
