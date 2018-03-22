#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Master Script to measure and report cross-talk levels among 12 ROE channels.
Takes as input a data-set composed of 3x12 CCD images, corresponding to 
injecting a "ladder" of signal on each of the 12 channels, using the ROE-TAB.


Created on Thu Mar 22 16:17:39 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import sys, os
from pdb import set_trace as stop
from astropy.io import ascii
import datetime
import numpy as np
from optparse import OptionParser

from vison.xtalk import xtalk
from vison.support.files import cPickleDumpDictionary,cPickleRead
from vison.support import logger as lg
# END IMPORT

#HARDWIRED
colstart = 1
colend = 1600
thresholdinj = 1.E2 # ADU
#END HARDWIRED


def run_xtalk(incat,inpath='',respath='',label='',doCompute=False):
    """ """
    
    ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')
    
    if label !='': _label = '_%s' % label
    else: _label = ''
    logfile = os.path.join(respath,'analysis_Xtak_%s%s.log' % (ttag,_label))
    
    outfile = os.path.join(respath,'Xtalks_%s%s.pick' % (ttag,_label))

    
    indata = ascii.read(incat)
    CHANNELS = indata['CHANNELS'].data.copy()
    OBSIDs = indata['OBSID'].data.copy()
    
    CCDs = [1,2,3]
    Quads = ['E','F','G','H']
    
    if doCompute: 
        
       if os.path.exists(logfile): os.system('rm %s' % logfile)
        
       log = lg.setUpLogger(logfile)
       log.info(['Starting CROSS-TALK ANALYSIS on %s' % ttag])
    
       Xtalks = dict()
    
       for CCDref in CCDs:
       
           Xtalks['CCD%i' % CCDref] = dict()
       
           for Qref in Quads:
                          
               source = '%i%s' % (CCDref,Qref)
               if source not in CHANNELS:
                   continue
               
               OBSID = OBSIDs[np.where(CHANNELS == source)]
                
               Xtalks['CCD%i' % CCDref][Qref] = xtalk.processXtalk_single(CCDref,Qref,OBSID,\
               thresholdinj,colstart=colstart,colend=colend,savefigs=True,log=log,
               datapath=inpath,respath=respath)               
                   
       cPickleDumpDictionary(Xtalks,outfile)
    
    else:
        
        Xtalks = cPickleRead(outfile)
        
    stop()


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-p","--path",dest="path",default='',help="Data Path.")
    parser.add_option("-r","--respath",dest="respath",default='',help="Results Path.")
    parser.add_option("-i","--incat",dest="incat",default='',help="Inputs Catalog-like file.")
    parser.add_option("-l","--label",dest="label",default='',help="label [figures, file-names].")
    parser.add_option("-C","--compute",dest="doCompute",action='store_true',default=False,help="Compute cross-talk matrix OR [default] only report results.")
        
    (options, args) = parser.parse_args()
        
    path = options.path
    respath = options.respath
    incat = options.incat
    label = options.label
    doCompute = options.doCompute
    
    if incat == '':
        parser.print_help()
        sys.exit()
    
    if respath != '' and not os.path.exists(respath):
        os.system('mkdir %s' % respath)
    

    header = '\n\n#########################################################\n'+\
             '#                                                       #\n'+\
             '#                                                       #\n'+\
             '#       running XTALK analysis                          #\n'+\
             '#                                                       #\n'+\
             '#                                                       #\n'+\
             '#########################################################\n'
    
    print header
    
    
    
    run_xtalk(incat,path,respath,label=label,doCompute=doCompute)
    