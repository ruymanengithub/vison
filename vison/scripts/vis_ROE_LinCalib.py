#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Non-Linearity Calibration of ROE (on bench).

Created on Thu Mar 15 15:32:11 2018

:author: Ruyman Azzollini

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
import pandas as pd

from vison.roe_fft import NL_lib
from vison.datamodel.ccd import CCD as CCDClass
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import vjson
#from vison.support.excel import ReportXL
from vison.datamodel import cdp
from vison.support import utils
from vison import __version__ as vison_version

from pylab import plot, show
#import matplotlib
# matplotlib.use('Agg')
from matplotlib import pyplot as plt

# END IMPORT


def find_adu_levels(qdata, Nlevels, debug=False):
    """ """
    nsamples = qdata.size

    kmeans = KMeans(n_clusters=Nlevels, verbose=0,
                    n_jobs=-1).fit(qdata.flatten().reshape(nsamples, 1))

    raw_levels = kmeans.cluster_centers_.flatten()
    adu_levels = np.sort(raw_levels)

    if debug:
        import pandas as pd
        import seaborn as sns
        anomaly = qdata - raw_levels[kmeans.labels_].reshape(qdata.shape)
        ixsel = (np.random.choice(anomaly.size, 20000),)
        fooy = anomaly.flatten()[ixsel].copy()
        foox = raw_levels[kmeans.labels_][ixsel].copy()
        df = pd.DataFrame.from_dict(dict(x=foox, y=fooy))
        #from astropy.io import fits as fts
        sns.jointplot(x="x", y="y", data=df, kind="kde")
        show()
        stop()

    return adu_levels


def run_ROE_LinCalib(inputsfile, incatfile, datapath='', respath='', doExtractFits=True,
                     dopolyRT=False, debug=False):
    """ """

    run_ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')
    degROE = 6

    inputs = vjson.load_jsonfile(inputsfile, useyaml=True)['inputs']
    Date = inputs['AcqDate']
    ROE = inputs['ROE']
    Injector_CalFile = inputs['Injector_CalFile']
    RTlevels = np.array(inputs['RTlevels'])
    rawtag = inputs['tag']
    if len(rawtag) == 0:
        tag = rawtag
    else:
        tag = '_%s' % rawtag

    if not os.path.exists(respath):
        os.system('mkdir %s' % respath)

    Nlevels = len(RTlevels)

    InjectorCal = cPickleRead(Injector_CalFile)

    # Master Inputs File

    indata = ascii.read(incatfile)

    CHANNELS = indata['CHANNEL'].data.copy()
    FitsList = indata['FITS'].data.copy()

    # Loop over CHANNELS to be calibrated

    if debug:
        ixchansel = np.where(CHANNELS == '1F')
        CHANNELS = CHANNELS[ixchansel]  # TESTS
        FitsList = FitsList[ixchansel]

    # Initialisations

    function, module = utils.get_function_module()
    CDP_header = OrderedDict()
    CDP_header['function'] = function
    CDP_header['module'] = module
    CDP_header['DATE'] = run_ttag
    CDP_header['vison version'] = vison_version
    CDP_header['tag'] = tag

    outfile = os.path.join('ROE_NLCAL_%s_%s%s' % (ROE, Date, tag))

    adu_levels_dict = OrderedDict()

    meta = OrderedDict(
        degROE=degROE)
    meta.update(inputs)
    meta['CHANNELS'] = CHANNELS.tolist()

    figs = OrderedDict()

    data = OrderedDict()

    data['ROE_NONLIN'] = OrderedDict()
    data['ROE_NONLIN']['CHANNEL'] = CHANNELS.copy()
    data['ROE_NONLIN']['MaxNL_pc'] = []
    data['ROE_NONLIN']['MaxNL_mV'] = []

    data['ROE_NONLIN_POLS'] = OrderedDict()
    data['ROE_NONLIN_POLS']['CHANNEL'] = CHANNELS.copy()
    for i in range(degROE + 1):
        data['ROE_NONLIN_POLS']['P%i' % i] = []

    for CHAN in CHANNELS:
        data[CHAN] = OrderedDict()
        data[CHAN]['RT_DN'] = RTlevels
        data[CHAN]['RT_mV'] = np.zeros_like(RTlevels, dtype='float32') + np.nan
        data[CHAN]['ADC_ADUS'] = np.zeros_like(RTlevels, dtype='float32') + np.nan
        data[CHAN]['fit_LS'] = np.zeros_like(RTlevels, dtype='float32') + np.nan
        data[CHAN]['NLpc_LS'] = np.zeros_like(RTlevels, dtype='float32') + np.nan

    # Main Loop

    degRT = InjectorCal['meta']['degRT']
    injCHANs = np.array(InjectorCal['meta']['CHANNELS'])

    for ix, CHAN in enumerate(CHANNELS):

        print('Processing Channel %s...' % CHAN)

        adus_pickf = os.path.join(respath, '%s_%s_%s_NL_ADUs.pick' % (ROE, Date, CHAN))

        Q = CHAN[1]

        ixRT = np.where(injCHANs == CHAN)

        if dopolyRT:

            RT_pol_DN2mV = np.zeros(degRT + 1)
            for ip in range(degRT + 1):
                RT_pol_DN2mV[ip] = InjectorCal['data']['RT_ADU_2_mV']['P%i' % ip].as_matrix()[
                    ixRT][0]
            mVfitlevels = np.polyval(RT_pol_DN2mV, RTlevels)

        else:
            assert np.all(RTlevels == InjectorCal['data'][CHAN]['RT_DN'].as_matrix())
            mVfitlevels = InjectorCal['data'][CHAN]['RT_mV'].as_matrix()
            if debug:
                #mVfitlevels = np.array([0.,120.,240.,350.,460.,580.,680.,780.,880.,1000.,1120.,1220.,1360.,1440.,1540.,1660.,1750.])
                # mVfitlevels = InjectorCal['data'][CHAN]['RT_DN'].as_matrix() # TESTS
                pass

        # EXTRACTING LEVELS FROM IMAGE (ADUs)

        if doExtractFits:

            fitsf = os.path.join(datapath, FitsList[ix])
            ccdobj = CCDClass(fitsf)
            qdata = ccdobj.get_quad(Q, canonical=True)

            adu_levels = find_adu_levels(qdata, Nlevels, debug=debug)
            adu_levels_dict[CHAN] = adu_levels

            cPickleDumpDictionary(adu_levels_dict, adus_pickf)

        else:

            adu_levels = cPickleRead(adus_pickf)[CHAN]

        # NON-LINEARITY OF ROE

        ixgood = np.where((adu_levels < 0.9 * 2.**16) & (mVfitlevels > 10.))

        # adu_levels -= adu_levels[0]  # BIAS subtraction

        data[CHAN]['RT_mV'] = mVfitlevels
        data[CHAN]['ADC_ADUS'] = adu_levels

        R_pol_NL, bag_NL = NL_lib.find_NL_pol(
            mVfitlevels[ixgood], adu_levels[ixgood], deg=degROE, Full=True, debug=debug)

        R_data_NL = bag_NL[4]

        xR = np.linspace(mVfitlevels[ixgood][0], mVfitlevels[ixgood][-1], 1000)

        R_bestfit_NL = np.polyval(R_pol_NL, xR) * 100.

        data[CHAN]['fit_LS'][ixgood] = bag_NL[2].copy()
        data[CHAN]['NLpc_LS'][ixgood] = R_data_NL * 100.

        for ip in range(degROE + 1):
            data['ROE_NONLIN_POLS']['P%i' % ip].append(R_pol_NL[ip])

        ixmaxNL = np.argmax(np.abs(R_bestfit_NL))
        MaxNL = R_bestfit_NL[ixmaxNL]
        MaxNL_mV = xR[ixmaxNL]

        data['ROE_NONLIN']['MaxNL_pc'].append(MaxNL)
        data['ROE_NONLIN']['MaxNL_mV'].append(MaxNL_mV)

        dyNL = bag_NL[1] - bag_NL[2]

        pldatasetANL = dict(
            data=dict(x=adu_levels[ixgood], y=dyNL, marker='o', ls='-', color='k'),
        )

        figs['%s_ANL' % CHAN] = os.path.join(
            respath, '%s_%s_%s_ROE_ABSNL%s.png' % (CHAN, ROE, Date, tag))

        if debug:
            figname1 = ''
        else:
            figname1 = figs['%s_ANL' % CHAN]

        NL_lib.myplot_NL(pldatasetANL, xlabel='ADU', ylabel='y-yLin [ADU]',  # ylim=[-10., 10.],
                         title='ROE Abs. NonLin. CHAN=%s' % CHAN,
                         figname=figname1)

        pldatasetRNL = dict(
            data=dict(x=mVfitlevels[ixgood], y=R_data_NL * 100., marker='o', color='k'),
            bestfit=dict(x=xR, y=R_bestfit_NL, ls='--', color='r')
        )

        figs['%s_NL' % CHAN] = os.path.join(
            respath, '%s_%s_%s_ROE_NL%s.png' % (CHAN, ROE, Date, tag))

        if debug:
            figname2 = ''
        else:
            figname2 = figs['%s_NL' % CHAN]

        NL_lib.myplot_NL(pldatasetRNL, xlabel='mV', ylabel='NL[percent]',  # ylim=[-10., 10.],
                         title='ROE NonLin. CHAN=%s' % CHAN,
                         figname=figname2)

    figs['keys'] = figs.keys()
    figs['jump'] = 26

    if debug:
        stop('Press c to finish... and potentially overwrite results')

    dddf = OrderedDict()
    for key in data.keys():
        dddf[key] = pd.DataFrame.from_dict(data[key])

    cdpRNL = cdp.Tables_CDP()
    cdpRNL.rootname = outfile
    cdpRNL.path = respath

    cdpRNL.ingest_inputs(
        data=dddf.copy(),
        meta=meta.copy(),
        header=CDP_header.copy(),
        figs=figs.copy()
    )

    cdpRNL.init_wb_and_fillAll(header_title='ROE NON LINEARITY Calibration Report')
    cdpRNL.savetopickle()
    cdpRNL.savehardcopy()


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path",
                      default='', help="Data Path.")
    parser.add_option("-d", "--debug", dest="debug",
                      action="store_true", default=False,
                      help="Debug mode")
    parser.add_option("-r", "--respath", dest="respath",
                      default='', help="Results Path.")
    parser.add_option("-c", "--incat", dest="incat",
                      default='', help="Inputs Catalog-like file.")
    parser.add_option("-i", "--inputs", dest="inputs",
                      default='', help="Master inputs file.")
    parser.add_option("-F", "--FITS", dest="doExtractFits",
                      action="store_true", default=False, help="Extract ADUs from FITS?")
    parser.add_option("-P", "--polynomial", dest="dopolyRT",
                      default=False, action="store_true",
                      help="Use Polynomial fit of RT-DN vs. RT-mV?")

    (options, args) = parser.parse_args()

    datapath = options.path
    respath = options.respath
    debug = options.debug
    incat = options.incat
    inputsfile = options.inputs
    doExtractFits = options.doExtractFits
    dopolyRT = options.dopolyRT

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

    print(header)

    run_ROE_LinCalib(inputsfile, incat, datapath=datapath, respath=respath,
                     doExtractFits=doExtractFits,
                     dopolyRT=dopolyRT,
                     debug=debug)
