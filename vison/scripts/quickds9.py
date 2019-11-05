#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python
# -*- coding: utf-8 -*-
"""

Wrap-up of ds9 to quickly load a number of images, for inspection.

:History:
Created on Thu Mar 17 13:18:10 2016

@author: Ruyman Azzollini
"""
from pdb import set_trace as stop
from optparse import OptionParser
import sys
import os
from glob import glob

#from vison.support.ds9 import ds9class
try:
    import pyds9
except ImportError:
    print 'pyds9 not installed in the system. Will not be possible to liaise with SAO-ds9.'


isthere = os.path.exists

if __name__ == '__main__':

    usage = 'usage: %prog [options] arg1 [arg2]\n' +\
            'arg1: starting OBSID\n' +\
            '[arg2]: optional, end OBSID. If not provided, arg2==arg1.'
    parser = OptionParser(usage)
    parser.add_option("-p", "--path", dest="path", default='',
                      help="path where to look for FITS files. Search is not recursive.")
    parser.add_option("-r", "--roe", dest="roe",
                      default=0, help="ROE to select")
    parser.add_option("-c", "--ccd", dest="ccd",
                      default=0, help="CCD to select")
    parser.add_option("-t", "--template", dest="temp",
                      default='EUC_%i_*D*T_ROE%i_CCD%i.fits',
                      help="Image Name Template")
    parser.add_option("-d", "--DS9", dest="ds9target", default="*",
                      help="Specify DS9 target (pyds9)?")
    #parser.add_option("-s","--start",dest="start",default=None,help="Start OBSID")
    #parser.add_option("-e","--end",dest="end",default=None,help="End OBSID")
    #parser.add_option("-e","--elvis",dest="elvis",default='5.7.02',help="ELVIS version")

    (options, args) = parser.parse_args()

    if ((len(args) < 1) or (len(args) > 2)):
        sys.exit("incorrect number of arguments")
    else:
        obsid0 = int(args[0])
        if len(args) == 2:
            obsid1 = int(args[1])
        else:
            obsid1 = obsid0

    path = options.path
    roe = int(options.roe)
    ccd = int(options.ccd)
    template = options.temp
    ds9target = options.ds9target
    #elvis = options.elvis

    FITSlist = []
    for iobs in range(obsid0, obsid1 + 1):
        tantFITS = os.path.join(
            path, template % (iobs, roe, ccd))

        try:
            FITS = glob(tantFITS)[0]
            FITSlist.append(FITS)
        except IndexError as IOError:
            print("Image %s not found!" % tantFITS)

    try:
        d = pyds9.DS9(target=ds9target)
    except NameError:
        print "pyds9 not installed, can't open DS9!"
        sys.exit()

    nFITS = len(FITSlist)
    for ifits, FITS in enumerate(FITSlist):

        d.set("frame %i" % (ifits + 1,))
        d.set('file %s' % FITS)
        d.set('zoom to fit')
        d.set('scale histequ zscale')
