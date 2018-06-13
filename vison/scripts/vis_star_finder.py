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
import itertools
import os
from astropy.io import ascii
from collections import OrderedDict
import numpy as np
import copy
import tempfile
import subprocess
import string as st

from vison.point.startracker import StarTracker
from vison.image.ds9reg import save_spots_as_ds9regs
from vison.image import sextractor as sex

# END IMPORT

Quads = ['E','F', 'G', 'H']
Starnames = ['ALPHA', 'BRAVO', 'CHARLIE', 'DELTA', 'ECHO']


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


def write_ID_chart(filename, Quads, Starnames):
    """ """
    
    StarKeys = [item[0] for item in Starnames]
    IDs = ['%s%s' % item for item in list(*[itertools.product(*[Quads,StarKeys])])]
    
    hdr = '# ID NUMBER'

    with open(filename, 'wa') as f: 
        print >> f, hdr
        for item in IDs: print >> f, item
    f.close()
    
    
def load_ID_chart(IDfile):    
    data = ascii.read(IDfile)
    iddict = OrderedDict()
    iddict['ID'] = data['ID'].data.copy()
    iddict['NUMBER'] = data['NUMBER'].data.copy()
    return iddict

def get_pattern_table(IDs_dict, SExCat):
    
    StarKeys = [item[0] for item in Starnames]
    IDs = ['%s%s' % item for item in list(*[itertools.product(*[Quads,StarKeys])])]
    IDs = np.array(IDs)
    
    X = []
    Y = []
    
    for ID in IDs:
        iNUMBER = IDs_dict['NUMBER'][np.where(IDs_dict['ID'] == ID)][0]
        ixsel = np.where(SExCat['NUMBER'] == iNUMBER)
        X.append(SExCat[ixsel]['X_IMAGE'].data[0])
        Y.append(SExCat[ixsel]['Y_IMAGE'].data[0])
    
    X = np.array(X)
    Y = np.array(Y)
    
    pattern_table = OrderedDict()
    pattern_table['ID'] = IDs.copy()
    pattern_table['X'] = X.copy()
    pattern_table['Y'] = Y.copy()
    
    return pattern_table

def run_starfinder(FITS, tag=''):
    
    if tag != '':
        _tag = '_%s' % tag
    else:
        _tag = ''
    
#    debug = False
#    
#    if debug:
#        
#        SExCatroot ='StarFinder_CCD2'
#        SExCat_dat_name = '%s.dat' % SExCatroot
#        SExCat = ascii.read(SExCat_dat_name,format='fixed_width')
#        CCDkey = 'CCD2'
#        withpover = True
#        # Display image and (filtered) detections as regions
#        
#        Pattfile = 'Pattern_%s%s.txt' % (CCDkey,_tag)
#                            
#        d = pyds9.DS9()
#        
#    else:
    
    with fts.open(FITS,lazy_load_hdus=False) as hdulist:
        
        EXTNAME = hdulist[1].header['EXTNAME']
        CCDkey = EXTNAME[EXTNAME.find('_')+1:]
        
        Pattfile = 'Pattern_%s%s.txt' % (CCDkey,_tag)
        
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
        
        SExCat_dat_name = '%s.dat' % SExCatroot
        
        SExCat.write(SExCat_dat_name, format='ascii.fixed_width')
        
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
    
        
    
    # USER IDENTIFIES POINT SOURCES
    
    fahrtig = False
    
    while not fahrtig:
        
        IDfile = '%s_IDENTIFY.cat' % SExCatroot
        
        if not os.path.exists(IDfile):
            write_ID_chart(IDfile, Quads, Starnames)
        
        proc1 = subprocess.Popen(["nedit",IDfile])
        proc2 = subprocess.Popen(["nedit",SExCat_dat_name])
        
        print 'Fill in the identification chart: %s\n' % IDfile
        print 'Use the SExCat %s and the image & segmentation map on DS9\n' % SExCat_dat_name
        
        holdon = raw_input('Press any key when finished!\n ')
        
        IDs_dict = load_ID_chart(IDfile)
        
        pattern_ccd_table = get_pattern_table(IDs_dict, SExCat)
        
        # Show identifications on ds9
        
        d.set("frame 1")
        d.set("regions delete all")
        
        _X = pattern_ccd_table['X'].copy()-1.
        IDregs = dict(
            X = _X,
            Y = pattern_ccd_table['Y'].copy()-1.,
            THETA = np.zeros_like(_X),
            A = np.ones_like(_X)*15.,
            B = np.ones_like(_X)*15.
            )
        idds9regsf = tempfile.NamedTemporaryFile(mode='w+a',suffix='_ds9regs.txt',
                                           prefix='vis_star_finder',delete=False)
        
        save_spots_as_ds9regs(IDregs,regfile=idds9regsf,regtype='ellipse',
                    clobber=True)
        idds9regsf.close()
        
        d.set('regions load %s' % idds9regsf.name)
        
        os.system('rm %s' % idds9regsf.name)
        
        whileender = raw_input('Are you satisfied with identifications? Y/n ')
        if whileender.upper() == 'Y':
            fahrtig = True
        
        proc1.terminate()
        proc2.terminate()
    
    CCD = int(st.replace(CCDkey,'CCD',''))
    stracker = StarTracker(CCD,withpover=withpover)
    
    Xc = pattern_ccd_table['X'].copy()
    Yc = pattern_ccd_table['Y'].copy()
    
    Xp, Yp = stracker.convert_Phys_2_CCD(Xc,Yc)
    
    pattern_phys_table = OrderedDict()
    pattern_phys_table['ID'] = pattern_ccd_table['ID'].copy()
    pattern_phys_table['X'] = Xp.copy()
    pattern_phys_table['Y'] = Yp.copy()
    
    
    stracker.load_Patt_fromdict(pattern_phys_table)
    
    if os.path.exists(Pattfile):
        hold = raw_input('Beware, you are about to overwrite "%s"... any key to proceed. ' % Pattfile)
    
    stracker.save_Patt_tofile(Pattfile, overwrite=True)
    
    

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-F","--FITS",dest="FITS",default='',help="input FITS file.")
    parser.add_option("-t","--tag",dest="tag",default='',help="output pattern file tag")
    
    (options, args) = parser.parse_args()
    
    if options.FITS == '':
        parser.print_help()
        sys.exit()
    
    run_starfinder(options.FITS, options.tag)
    