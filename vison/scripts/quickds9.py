#! $HOME/SOFTWARE/anaconda/envs/VISSIM/bin/ python
# -*- coding: utf-8 -*-
"""

Wrap-up of ds9 to quickly load a number of images, for inspection.

:History:
Created on Thu Mar 17 13:18:10 2016

@author: Ruyman Azzollini
"""

from optparse import OptionParser
import sys
import os
from glob import glob
import string as st
from pdb import set_trace as stop

from vison.support.ds9 import ds9class


isthere = os.path.exists

if __name__ == '__main__':
    """ 
     
    """
    usage='usage: %prog [options] arg1 [arg2]'
    parser = OptionParser(usage)
    parser.add_option("-p","--path",dest="path",default='',help="path where to look for FITS files. Search is not recursive.")    
    parser.add_option("-r","--roe",dest="roe",default=0,help="ROE to select")
    parser.add_option("-c","--ccd",dest="ccd",default=0,help="CCD to select")
    #parser.add_option("-s","--start",dest="start",default=None,help="Start OBSID")
    #parser.add_option("-e","--end",dest="end",default=None,help="End OBSID")
    #parser.add_option("-e","--elvis",dest="elvis",default='5.7.02',help="ELVIS version")
        
    (options, args) = parser.parse_args()
    
    if ((len(args) <1) or (len(args)>2)): sys.exit("incorrect number of arguments")
    else:
        obsid0 = int(args[0])
        if len(args) == 2:
            obsid1 = int(args[1])
        else:
            obsid1 = obsid0
    
    
    path = options.path
    roe = int(options.roe)
    ccd = int(options.ccd)
    #elvis = options.elvis
    
    
    FITSlist = []
    for iobs in range(obsid0,obsid1+1):
        tantFITS = os.path.join(path,'EUC_%i_*D_*T_ROE%i_CCD%i.fits' % (iobs,roe,ccd))
        try:
            FITS = glob(tantFITS)[0]
            FITSlist.append(FITS)
        except:
            sys.exit("Image %s not found!" % tantFITS)
    
    
    ds9obj = ds9class()
    if ds9obj.isOpen():
        sys.exit('Close all OPEN DS9 instances first!')
    ds9obj.launch()
    
    nFITS = len(FITSlist)
    for ifits,FITS in enumerate(FITSlist):
        
        
        ds9obj.xpaset('file %s' % FITS)
        if ifits ==0: ds9obj.xpaset('file %s' % FITS)
        ds9obj.xpaset('zoom to fit')
        ds9obj.xpaset('scale histequ zscale')
        if ifits < (nFITS-1): 
            ds9obj.xpaset('frame new')
    
