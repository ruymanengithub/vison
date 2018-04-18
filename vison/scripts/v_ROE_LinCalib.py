#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 15:32:11 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from itertools import islice
import numpy as np
from pdb import set_trace as stop
from astropy.io import ascii
from scipy import ndimage, stats
from sklearn.cluster import KMeans
import os
import copy
from optparse import OptionParser
from collections import OrderedDict
import sys
import datetime

import NL_lib
from vison.datamodel.ccd import CCD as CCDClass
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import vjson
from vison.support.excel import ReportXL

from matplotlib import pyplot as plt

# END IMPORT


def find_adu_levels(qdata, Nlevels):
    """ """
    nsamples = qdata.size

    kmeans = KMeans(n_clusters=Nlevels, verbose=0,
                    n_jobs=-1).fit(qdata.flatten().reshape(nsamples, 1))

    raw_levels = kmeans.cluster_centers_.flatten()
    adu_levels = np.sort(raw_levels)

    return adu_levels


class ReportXL_ROENL(ReportXL):

    def __init__(self, datadict):
        """ """
        super(ReportXL_ROENL, self).__init__(datadict)

        self.wb.create_sheet('Header', 0)
        self.wb.create_sheet('Data', 1)
        self.wb.create_sheet('ROE_NONLINpol', 2)
        self.wb.create_sheet('ROE_NONLIN', 3)
        self.wb.create_sheet('Figures', 4)


def run_ROE_LinCalib(inputsfile, incatfile, datapath='', respath='', doExtractFits=True):
    """ """

    run_ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')
    degROE = 9

    inputs = vjson.load_jsonfile(inputsfile)['inputs']
    Date = inputs['Date']
    ROE = inputs['ROE']
    Injector_CalFile = inputs['Injector_CalFile']
    RTlevels = np.array(inputs['RTlevels'])

    if not os.path.exists(respath):
        os.system('mkdir %s' % respath)

    Nlevels = len(RTlevels)

    InjectorCal = cPickleRead(Injector_CalFile)

    # Master Inputs File

    indata = ascii.read(incatfile)

    CHANNELS = indata['CHANNEL'].data.copy()
    FitsList = indata['FITS'].data.copy()

    # Loop over CHANNELS to be calibrated

    CHANNELS = [CHANNELS[0]]  # TESTS

    # Initialisations

    outexcelfile = os.path.join(respath, 'ROE_NLCAL_%s_%s.xlsx' % (ROE, Date))

    adus_pickf = os.path.join(respath, '%s_%s_NL_ADUs.pick' % (ROE, Date))
    adu_levels_dict = OrderedDict()

    figs = dict()

    reportdict = OrderedDict()
    reportdict['meta'] = OrderedDict(degROE=degROE,
                                     Analysis_Date=run_ttag)
    reportdict['meta'].update(inputs)

    reportdict['CHANNELS'] = CHANNELS

    reportdict['data'] = OrderedDict()
    reportdict['data']['RT_DN'] = RTlevels
    reportdict['data']['RT_V'] = []
    reportdict['data']['ADU'] = []

    reportdict['ROE_NONLIN_pol'] = []

    reportdict['ROE_NONLIN'] = OrderedDict()
    reportdict['ROE_NONLIN']['MaxNL'] = []
    reportdict['ROE_NONLIN']['MaxNL_ADU'] = []

    # Main Loop

    injCHANs = InjectorCal['CHANNELS']

    for ix, CHAN in enumerate(CHANNELS):

        Q = CHAN[1]

        ixRT = injCHANs.index(CHAN)

        RT_pol_V2DN = InjectorCal['RT_pol_V2DN'][ixRT]

        # EXTRACTING LEVELS FROM IMAGE (ADUs)

        if doExtractFits:

            fitsf = os.path.join(datapath, FitsList[ix])
            ccdobj = CCDClass(fitsf)
            qdata = ccdobj.get_quad(Q, canonical=True)

            adu_levels = find_adu_levels(qdata, Nlevels)
            adu_levels_dict[CHAN] = adu_levels

            cPickleDumpDictionary(adu_levels_dict, adus_pickf)

        else:

            adu_levels = cPickleRead(adus_pickf)[CHAN]

        # NON-LINEARITY OF ROE

        ixgood = np.where(adu_levels < 2.**16-2.)

        adu_levels -= adu_levels[0]  # BIAS subtraction

        Vfitlevels = NL_lib.f_pol_byorigin(RTlevels, *RT_pol_V2DN)

        reportdict['data']['RT_V'].append(Vfitlevels)
        reportdict['data']['ADU'].append(adu_levels)

        R_pol_NL = NL_lib.find_NL_pol(Vfitlevels[ixgood], adu_levels[ixgood],
                                      deg=degROE)

        R_pol_LIN = NL_lib.fit_pol_byorigin(Vfitlevels[ixgood], adu_levels[ixgood],
                                            deg=1)

        R_data_LIN = NL_lib.f_pol_byorigin(Vfitlevels, *R_pol_LIN)

        R_data_NL = (adu_levels - R_data_LIN)/R_data_LIN * 100.

        xR = np.linspace(Vfitlevels[0], Vfitlevels[-1], 1000)

        R_bestfit_NL = NL_lib.f_pol_byorigin(xR, *R_pol_NL) * 100.

        reportdict['ROE_NONLIN_pol'].append(R_pol_NL)

        ixmaxNL = np.argmax(R_bestfit_NL)
        MaxNL = R_bestfit_NL[ixmaxNL]
        MaxNL_ADU = xR[ixmaxNL]

        reportdict['ROE_NONLIN']['MaxNL'].append(MaxNL)
        reportdict['ROE_NONLIN']['MaxNL_ADU'].append(MaxNL_ADU)

        pldatasetRNL = dict(
            data=dict(x=Vfitlevels, y=R_data_NL, marker='o', color='k'),
            bestfit=dict(x=xR, y=R_bestfit_NL, ls='--', color='r')
        )

        figs[CHAN] = os.path.join(
            respath, '%s_%s_%s_ROE_NL.png' % (CHAN, ROE, Date))

        NL_lib.myplot_NL(pldatasetRNL, xlabel='V', ylabel='NL[pc]', ylim=[-10., 10.],
                         title='ROE NonLin. CHAN=%s' % CHAN,
                         figname=figs[CHAN])

    reportdict['figs'] = figs

    report = ReportXL_ROENL(reportdict)

    report.fill_Header()
    report.fill_Data()
    report.fill_NONLINpols()
    report.fill_NONLINMaxima()
    report.fill_Figures()

    report.save(outexcelfile)


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path",
                      default='', help="Data Path.")
    parser.add_option("-r", "--respath", dest="respath",
                      default='', help="Results Path.")
    parser.add_option("-c", "--incat", dest="incat",
                      default='', help="Inputs Catalog-like file.")
    parser.add_option("-i", "--inputs", dest="inputs",
                      default='', help="Master inputs file.")
    parser.add_option("-F", "--FITS", dest="doExtractFits",
                      action="store_true", default=False, help="Extract ADUs from FITS?")

    (options, args) = parser.parse_args()

    datapath = options.path
    respath = options.respath
    incat = options.incat
    inputsfile = options.inputs
    doExtractFits = options.doExtractFits

    if incat == '' or inputsfile == '':
        parser.print_help()
        sys.exit()

    if respath != '' and not os.path.exists(respath):
        os.system('mkdir %s' % respath)

    header = '\n\n#########################################################\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#       running ROE NON-LINEARITY ANALYSIS              #\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#########################################################\n'

    print header

    run_ROE_LinCalib(inputsfile, incat, datapath=datapath, respath=respath,
                     doExtractFits=doExtractFits)