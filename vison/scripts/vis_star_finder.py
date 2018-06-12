#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to find point sources in VIS Ground Calibration Campaign.
Used to 'prime' the position tables of point-source objects.

Created on Tue Jun 12 16:09:31 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
from optparse import OptionParser
import sys
from astropy.io import fits as fts
import pyds9

from vison.image.ds9reg import save_spots_as_ds9regs
from vison.image import sextractor as sex

# END IMPORT


def run_SEx(img, catroot):
    
    vSEx = sex.VSExtractor(img=img.copy())
    vSEx.internals['params'] = ['NUMBER','EXT_NUMBER','X_IMAGE','Y_IMAGE',
    'A_IMAGE','B_IMAGE','THETA_IMAGE','ELONGATION','FWHM_IMAGE','MAG_AUTO']
    vSEx.internals.update(
            dict(MINAREA=3,
                 DET_THRESH=13.,
                 MAG_ZERPOINT=20.,
                 SATUR_LEVEL=65535.,
                 SEEING_FWHM=1.2,
                 PIXEL_SCALE=1.,
                 GAIN=1.
                    )
            )
    
    SExCatFile = vSEx.run_sex(catroot,
                          checks=['BACKGROUND','SEGMENTATION'],                         
                          cleanafter=False)
    SExCat = sex.aw.utils.ldac.get_table_from_ldac(SExCatFile)
    
    return SExCat



def run_starfinder(FITS):
    
    
    with fts.open(FITS,lazy_load_hdus=False) as hdulist:
        
        EXTNAME = hdulist[1].header['EXTNAME']
        CCDkey = EXTNAME[EXTNAME.find('_')+1:]
        
        img = hdulist[1].data.copy()
        NAXIS2, NAXIS1 = img.shape
        assert NAXIS1 == 4238
        assert NAXIS2 in (4132, 4172)
        
        if NAXIS2/2 == 2086:
            withpover = True
        else:
            withpover = False
        
        SExCatroot ='StarFinder_%s' % CCDkey
        
        
        SExCat = run_SEx(img, SExCatroot)
        
        SEG_FITS = '%s_SEGM.fits' % SExCatroot
        BGD_FITS = '%s_BACK.fits' % SExCatroot

        
        SExCat.write('%s.dat' % SExCatroot, format='ascii.fixed_width')
        
        hdulist.close()
    
    # Creating ds9 regions for all detections
    
    ds9regsf = '%s_ds9.reg' % SExCatroot
    regsdata = dict(
        X = SExCat['X_IMAGE'].data.copy()-1.,
        Y = SExCat['Y_IMAGE'].data.copy()-1.,
        THETA = SExCat['THETA_IMAGE'].data.copy(),
        A = 2.*SExCat['A_IMAGE'].data.copy(),
        B = 2.*SExCat['B_IMAGE'].data.copy()
        )
    save_spots_as_ds9regs(regsdata,regfilename=ds9regsf,regtype='ellipse',
                clobber=True)
    
    # Display image and (filtered) detections as regions
    
    d = pyds9.DS9()
    d.set("frame 1")
    d.set("file %s" % FITS)
    d.set("zoom to fit")
    d.set("scale mode zscale")
    d.set('regions load %s' % ds9regsf)
    
    d.set("frame 2")
    d.set("file %s" % SEG_FITS)
    d.set("zoom to fit")
    d.set("scale log zscale")
    
    d.set("frame 3")
    d.set("file %s" % BGD_FITS)
    d.set("zoom to fit")
    d.set("scale mode zscale")
    
    stop()


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-F","--FITS",dest="FITS",default='',help="input FITS file.")
    
    (options, args) = parser.parse_args()
    
    if options.FITS == '':
        parser.print_help()
        sys.exit()
    
    run_starfinder(options.FITS)
    