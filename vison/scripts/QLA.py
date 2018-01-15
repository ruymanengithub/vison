#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python
# -*- coding: utf-8 -*-
"""


QLA script.
Aimed at quick inspection of data from Characterization and Calibration Campaigns
of Euclid-VIS.

:History:
Created on Wed Mar 16 11:33:21 2016

@author: Ruyman Azzollini

"""

from optparse import OptionParser
import sys
import os
from glob import glob
import string as st
import numpy as np

#from vissim.datamodel.HKtools import parseHKfiles,allHK_keys
#from vison.support.latex import LaTeX
from vison.datamodel import QLAtools as QLAt
from vison.pipe import lib as pilib
from pdb import set_trace as stop

isthere = os.path.exists


if __name__ == '__main__':
    """ 
    # TODO:
    #  find FITS files in a folder
    #  load FITS file into a CCD object
    #  obtain metrics on the data:
    #    image of 4 quadrants
    #    across-rows and across-columns plots
    #    statistics of image, pre and over scan regions: mea, med, std, p25, p75
    #  parse header
    #  do plots
    #  assemble all plots into a pdf file per image
    #  merge all pdfs into a single pdf file
     
    """
    
    parser = OptionParser()
    parser.add_option("-p","--path",dest="path",default='',help="path where to look for FITS files. Search is not recursive.")    
    parser.add_option("-r","--roe",dest="roe",default=None,help="ROE to select")
    parser.add_option("-c","--ccd",dest="ccd",default=None,help="CCD to select")
    #parser.add_option("-e","--elvis",dest="elvis",default='5.7.02',help="ELVIS version")
        
    (options, args) = parser.parse_args()
    
    
    if options.path == '':
        parser.print_help()
        sys.exit()
    
    path = options.path
    roe = int(options.roe)
    ccd = int(options.ccd)
    #elvis = options.elvis
  
    if not os.path.exists(path):
        sys.exit('QLA.py: %s does not exist' % path)  
    
    FITSlist = glob(os.path.join(path,'EUC_*_ROE%i_CCD%i.fits' % (roe,ccd)))
    #FITSlist = FITSlist[0:-1]
    
    nOBSIDs = len(FITSlist)
    firstFITS = os.path.split(FITSlist[0])[-1]
    datetag = st.split(firstFITS,'_')[2]
    datetag = datetag[0:datetag.index('D')]
    
    outpath = 'FITSmonitor_%s' % datetag
    if not isthere(outpath): os.system('mkdir %s' % outpath)
    
    
    PDFlist = []
    
    #ixstart = [ix for ix in range(nOBSIDs) if 'EUC_799' in FITSlist[ix]][0] # TESTS
    #PDFlist = glob(os.path.join(outpath,'EUC_*.pdf'))    
    
    #for FITSfile in FITSlist[ixstart:]:
    
    OBSIDs = []
    
    for FITSfile in FITSlist:  
        
        bareFITS = os.path.split(FITSfile)[-1]
        OBSID = int(st.split(bareFITS,'_')[1])
        OBSIDs.append(OBSID)
        
        FITSpdf = QLAt.reportFITS(FITSfile,outpath=outpath)
        PDFlist.append(FITSpdf)
        
    
    order = np.argsort(OBSIDs)
    PDFlist = np.array(PDFlist)[order]
    
    pdfFile = 'FITSmonitor_%s_ROE%s_CCD%s.pdf' % (datetag,roe,ccd)
    
    
    tempcommand = 'ghostscript -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s %s'
    PDFlistLine = st.join(PDFlist,' ')
    os.system(tempcommand % (pdfFile,PDFlistLine))
    
    for PDFfile in PDFlist: os.system('rm %s' % PDFfile)
    
    os.system('mv %s %s/' % (pdfFile,outpath))
