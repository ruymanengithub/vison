#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python

# -*- coding: utf-8 -*-
"""
:TODO:
     find HK files in a folder
     parse HK files
     plot HK parameters vs. time
     assemble all plots into a pdf file
     
     DEBUG, calls nonexistent class LaTeX

Script to produce HK reports out of HK files in a folder.
Aimed at quick inspection of data from Characterization and Calibration Campaigns
of Euclid-VIS.




:History:
Created on Tue Mar 15 10:35:43 2016

@author: Ruyman Azzollini (MSSL)
"""

from optparse import OptionParser
import sys
import os
from glob import glob
import string as st
import numpy as np

from vison.support import context
#from vison.support.latex import LaTeX
#from vison.pipe import lib as pilib
from vison.datamodel.HKtools import doHKSinglePlot, parseHKfiles, allHK_keys, HKlims
#from vison.support.latex import LaTeX
from vison.support import report as repmod

from pdb import set_trace as stop

isthere = os.path.exists

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path", default='',
                      help="path where to look for HK files. Search is not recursive.")
    parser.add_option("-r", "--roe", dest="roe",
                      default=None, help="ROE to select")
    parser.add_option("-e", "--elvis", dest="elvis",
                      default=context.elvis, help="ELVIS version")
    parser.add_option("-o", "--obsid", dest="OBSID",
                      default=None, help="Starting OBSID")
    parser.add_option("-n", "--Nobs", dest="Nobs", default=-
                      1, help="Number of OBSIDs to plot")

    (options, args) = parser.parse_args()

    if options.path == '':
        parser.print_help()
        sys.exit()
    
    debug = False # REDTAG
    
    path = options.path
    roe = int(options.roe)
    elvis = options.elvis
    OBSID = options.OBSID
    Nobs = int(options.Nobs)

    if not os.path.exists(path):
        sys.exit('HKmonitory.py: %s does not exist' % path)

    HKlist = glob(os.path.join(path, 'HK_*_ROE%i.txt' % roe))
    HKlist.sort()

    if OBSID is not None:

        allOBSIDS = [st.split(item, '_')[3] for item in HKlist]
        ixstart = allOBSIDS.index(OBSID)

        if Nobs > 0:
            ixend = ixstart + Nobs
            HKlist = HKlist[ixstart:ixend]
        else:
            HKlist = HKlist[ixstart:]

    nOBSIDs = len(HKlist)
    firstHK = os.path.split(HKlist[0])[-1]
    datetag = st.split(firstHK, '_')[2]
    datetag = datetag[0:datetag.index('D')]

    outpath = 'HKmonitor_%s' % datetag
    if not isthere(outpath):
        os.system('mkdir %s' % outpath)
    
    if not debug:
        
        obsids, dtobjs, tdeltasec, HK_keys, HKdata = parseHKfiles(
            HKlist, elvis=elvis)

        fOBSID = obsids[0]
        lOBSID = obsids[-1]
    
        figlist = []
    
        HKkeys = allHK_keys[elvis]
        if 'TimeStamp' in HKkeys:
            HKkeys.pop(HKkeys.index('TimeStamp'))
        

        for ik,key in enumerate(HKkeys):
    
            figname = os.path.join(outpath, '%s_%s.eps' % (datetag, key,))
            
            HKvec = HKdata[:,0,ik+1].copy()
            
            iHKlims = HKlims[elvis]['S'][key]
            
            doHKSinglePlot(dtobjs, HKvec,
                           HKkey=key, ylabel=key, HKlims=iHKlims[1:], 
                           filename=figname)
            
            figlist.append(figname)    

    reportroot = 'HKmonitor_%s_%i-%i_OBSIDs_ROE%i' % (
        datetag, fOBSID, lOBSID, roe)
     
    report = repmod.Report(TestName='HKMonitor_%s' % datetag,
                           Model='XX')
    
    report.add_FigsTable(figlist,Ncols=2,figswidth='8',
                         caption='HK plots, built from HK-OBSID files only.'+\
                         ' CAUTION: HK collected inter-readouts is not shown here.')
    report.doreport(reportroot,cleanafter=True)
    os.system('mv %s.pdf %s/' % (reportroot, outpath))
