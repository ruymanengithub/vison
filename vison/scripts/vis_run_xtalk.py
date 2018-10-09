#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Master Script to measure and report cross-talk levels among 12 ROE channels.
Takes as input a data-set composed of 3x12 CCD images, corresponding to 
injecting a "ladder" of signal on each of the 12 channels, using the ROE-TAB.


Created on Thu Mar 22 16:17:39 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import sys
import os
from pdb import set_trace as stop
from astropy.io import ascii
import datetime
import numpy as np
from optparse import OptionParser
from collections import OrderedDict
import copy

from vison.datamodel import cdp
from vison.xtalk import xtalk
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import logger as lg
from vison.support import vjson
# END IMPORT

# HARDWIRED
rowstart = 1
rowend = 2086
colstart = 1
colend = 1600
thresholdinj = 1.E2  # ADU
# END HARDWIRED

meta_defaults = dict(Date="21.02.1980", ROE="Unknown",
                     ROE_FW="Uknown",
                     Injector="Uknown",
                     Injector_FW="Uknown",
                     Label=None)


def run_xtalk(incat, inpath='', respath='', metafile='', doCompute=False):
    """ """

    run_ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')

    meta = copy.deepcopy(meta_defaults)
    if metafile != '':
        inmeta = vjson.load_jsonfile(metafile)["inputs"]
        meta.update(inmeta)

    datetag = meta["AcqDate"]

    label = meta['Label']
    if label != '':
        _label = '_%s' % label
    else:
        _label = ''

    logfile = os.path.join(
        respath, 'analysis_Xtalk_%s%s.log' % (datetag, _label))
    outrootname = 'Xtalks_%s%s' % (datetag,_label)
    outpickfile = os.path.join(respath, '%s.pick' % outrootname)
    outexcelfile = os.path.join(respath, '%s.xlsx' % outrootname)

    indata = ascii.read(incat)
    CHANNELS = indata['CHANNELS'].data.copy()
    OBSIDs = indata['OBSID'].data.copy()

    CCDs = [1, 2, 3]
    Quads = ['E', 'F', 'G', 'H']
    
    
    if doCompute:

        if os.path.exists(logfile):
            os.system('rm %s' % logfile)

        log = lg.setUpLogger(logfile)
        log.info(['Starting CROSS-TALK ANALYSIS on %s' % run_ttag])

        Xtalks = dict(meta=meta)
        Xtalks['figs'] = dict()
        Xtalks['figs']['keys'] = []

        for CCDref in CCDs:

            Xtalks['CCD%i' % CCDref] = dict()

            for Qref in Quads:

                source = '%i%s' % (CCDref, Qref)
                if source not in CHANNELS:
                    continue

                OBSID = OBSIDs[np.where(CHANNELS == source)]

                Xtalks['CCD%i' % CCDref][Qref], fignames = xtalk.processXtalk_single(CCDref, Qref, OBSID,
                    thresholdinj, rowstart=rowstart, rowend=rowend, colstart=colstart, colend=colend,
                    savefigs=True, log=log, datapath=inpath, respath=respath)
                
                for CCD in CCDs:
                    Xtalks['figs']['CCD%i_R%i%s' % (CCD,CCDref,Qref)] = fignames['CCD%i' % CCD]
                    Xtalks['figs']['keys'].append('CCD%i_R%i%s' % (CCD,CCDref,Qref))
                

        cPickleDumpDictionary(Xtalks, outpickfile)

    else:

        Xtalks = cPickleRead(outpickfile)

    savefigs = True
    if savefigs:
        figname_R = os.path.join(
            respath, 'XTALK_summary_%s%s_RATIO.png' % (datetag, _label))
        figname_A = os.path.join(
            respath, 'XTALK_summary_%s%s_ADU.png' % (datetag, _label))
    else:
        figname_R = ''
        figname_A = ''

    suptitle = 'XTALK %s %s' % (label, datetag)
    xtalk.PlotSummaryFig(Xtalks, suptitle, figname_R, scale='RATIO')
    xtalk.PlotSummaryFig(Xtalks, suptitle, figname_A, scale='ADU')
    

    Xtalks['figs'].update(dict(
    RATIO=figname_R, ADU=figname_A))
    Xtalks['figs']['jump']=40
    Xtalks['figs']['keys'] = ['RATIO','ADU'] + Xtalks['figs']['keys']
    
    Xtalks['meta']['Analysis_Date'] = run_ttag

    report = xtalk.ReportXL_Xtalk(Xtalks)

    report.fill_Header()
    report.fill_Summary()
    report.fill_Xtalk_Ratios()
    report.fill_Xtalk_ADUs()
    report.fill_Figures()

    report.save(outexcelfile)
    


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path",
                      default='', help="Data Path.")
    parser.add_option("-r", "--respath", dest="respath",
                      default='', help="Results Path.")
    parser.add_option("-i", "--incat", dest="incat",
                      default='', help="Inputs Catalog-like file.")
    parser.add_option("-m", "--meta", dest="metafile",
                      default='', help="File with Meta-Information.")
    #parser.add_option("-l","--label",dest="label",default='',help="label [figures, file-names].")
    #parser.add_option("-t","--ttag",dest="ttag",default='',help="Time-Tag [figures, file-names].")
    parser.add_option("-C", "--compute", dest="doCompute", action='store_true',
                      default=False, help="Compute cross-talk matrix OR [default] only report results.")

    (options, args) = parser.parse_args()

    path = options.path
    respath = options.respath
    incat = options.incat
    metafile = options.metafile
    #label = options.label
    #ttag = options.ttag
    doCompute = options.doCompute

    if incat == '':
        parser.print_help()
        sys.exit()

    if respath != '' and not os.path.exists(respath):
        os.system('mkdir %s' % respath)

    header = '\n\n#########################################################\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#       running XTALK analysis                          #\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#########################################################\n'

    print header

    run_xtalk(incat, path, respath, metafile=metafile, doCompute=doCompute)
