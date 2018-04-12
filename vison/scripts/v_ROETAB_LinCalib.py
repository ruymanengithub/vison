#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Linearity Calibration of ROE-TAB.

Created on Tue Mar 27 14:42:00 2018

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
import sys
from optparse import OptionParser
from collections import OrderedDict
import pandas
import datetime

from matplotlib import pyplot as plt

import NL_lib
import Wave_bayes as WaveBay
from vison.datamodel.ccd import CCD as CCDClass
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import vjson
from vison.support.excel import ReportXL
from vison import __version__
# END IMPORT


def load_WF(WFf, chkNsamp=None, chkSampInter=None):
    """ """
    Nhead = 17
    with open(WFf) as myfile:
        hdr = list(islice(myfile, Nhead))
    Nsamp = int(hdr[0].split(',')[1])
    SampInter = float(hdr[1].split(',')[1])

    if chkNsamp is not None:
        assert Nsamp == chkNsamp
    if chkSampInter is not None:
        assert np.isclose(chkSampInter, SampInter)

    WF = ascii.read(WFf)
    timex = WF['col4'].data.copy()
    Vy = WF['col5'].data.copy() * (-1.)

    return (timex, Vy)


def plot_waveform(WF, disc_voltages=[], figname='', chan='Unknown'):
    """ """

    timex, Vy = WF

    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.plot(timex, Vy)

    for idislev in disc_voltages:
        ax.axhline(y=idislev, ls='--', color='r')

    ax.set_xlabel('time [s]')
    ax.set_ylabel('Voltage')
    ax.set_title('Linearity Calib. ROE-TAB. CHAN=%s' % chan)

    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()


def filter_Voltage(rV, filt_kernel):
    """ """
    fV = ndimage.filters.median_filter(rV, filt_kernel)
    return fV


def find_discrete_voltages_inwaveform(rV, levels, filtered=None):
    """ """

    if filtered is not None:
        assert len(rV) == len(filtered)
        iV = copy.deepcopy(filtered)
    else:
        iV = copy.deepcopy(rV)

    Nsamp = len(iV)
    Nclust = len(levels)

    kmeans = KMeans(n_clusters=Nclust, verbose=0,
                    n_jobs=-1).fit(iV.reshape(Nsamp, 1))

    labels = kmeans.labels_.copy()
    ulabels = np.unique(labels)

    discrete_levels = np.zeros(len(ulabels), dtype='float32')

    for i in range(len(ulabels)):
        discrete_levels[i] = np.mean(stats.sigmaclip(
            rV[np.where(labels == ulabels[i])], 3, 3).clipped)

    #discrete_levels = kmeans.cluster_centers_.flatten()

    sorted_levels = np.sort(discrete_levels)

    return sorted_levels


class ReportXL_RTLIN(ReportXL):
    """ """

    def __init__(self, datadict):
        """ """
        super(ReportXL_RTLIN, self).__init__(datadict)

        self.wb.create_sheet('Header', 0)
        self.wb.create_sheet('Data', 1)
        self.wb.create_sheet('RTADU_to_V', 2)
        self.wb.create_sheet('RT_NONLINpol', 3)
        self.wb.create_sheet('RT_NONLIN', 4)
        self.wb.create_sheet('Figures', 5)

    def fill_Header(self):
        """ """
        headerdict = self.data['meta']
        headerdict.update(dict(vison=__version__))
        headkeys = ['Date', 'Injector', 'Injector_FW', 'degRT', 'RTlevels', 'SampInter',
                    'Analysis_Date', 'vison']

        self.dict_to_sheet(headerdict, 'Header', keys=headkeys,
                           title='ROE-TAB Non-Linearity Calibration Report')

    def fill_Data(self):
        """ """
        CHANNELS = self.data['CHANNELS']
        RT_DN = self.data['data']['RT_DN']
        RT_V = self.data['data']['RT_V']
        RT_eV = self.data['data']['RT_eV']

        indict = OrderedDict()

        indict['DN'] = RT_DN

        for i, CHAN in enumerate(CHANNELS):
            indict['V_%s' % CHAN] = RT_V[i]
            indict['eV_%s' % CHAN] = RT_eV[i]

        df = pandas.DataFrame.from_dict(indict)
        self.df_to_sheet(df, 'Data', index=False, header=True)

    def fill_RTADU2V(self):
        """ """

        CHANNELS = self.data['CHANNELS']
        RT_ADU_2_V = self.data['RT_ADU_2_V']
        Npol = len(RT_ADU_2_V[0])

        indict = OrderedDict()
        indict['CHANNEL'] = CHANNELS

        for jCHAN, CHAN in enumerate(CHANNELS):

            for i in range(Npol):
                indict['p_%i' % (i+1)] = RT_ADU_2_V[jCHAN][i]

        df = pandas.DataFrame.from_dict(indict)
        self.df_to_sheet(df, 'RTADU_to_V', index=False, header=True)
        self.adjust_columns_width('RTADU_to_V', minwidth=10)

    def fill_RTNONLINpols(self):
        """ """

        CHANNELS = self.data['CHANNELS']
        RT_NONLIN_pol = self.data['RT_NONLIN_pol']
        Npol = len(RT_NONLIN_pol[0])

        indict = OrderedDict()
        indict['CHANNEL'] = CHANNELS

        for jCHAN, CHAN in enumerate(CHANNELS):

            for i in range(Npol):
                indict['p_%i' % (i+1)] = RT_NONLIN_pol[jCHAN][i]

        df = pandas.DataFrame.from_dict(indict)
        self.df_to_sheet(df, 'RT_NONLINpol', index=False, header=True)

        self.adjust_columns_width('RT_NONLINpol', minwidth=10)

    def fill_RTNONLINMaxima(self):
        """ """

        CHANNELS = self.data['CHANNELS']
        RT_NONLIN = self.data['RT_NONLIN']

        indict = OrderedDict()
        indict['CHANNEL'] = CHANNELS

        for jCHAN, CHAN in enumerate(CHANNELS):
            indict['MaxNL'] = RT_NONLIN['MaxNL'][jCHAN]
            indict['MaxNL_DN'] = RT_NONLIN['MaxNL_DN'][jCHAN]

        df = pandas.DataFrame.from_dict(indict)
        self.df_to_sheet(df, 'RT_NONLIN', index=False, header=True)

        self.adjust_columns_width('RT_NONLIN', minwidth=10)

    def fill_Figures(self):
        """ """
        figsdict = self.data['figs']
        CHANNELS = self.data['CHANNELS']

        rowcounter = 1
        rowjump = 26

        for iCHAN, CHAN in enumerate(CHANNELS):

            self.add_image(figsdict['WF'][CHAN], 'Figures', 'A%i' % rowcounter)
            self.add_image(figsdict['NL_RT'][CHAN],
                           'Figures', 'M%i' % rowcounter)

            rowcounter += rowjump


def run_ROETAB_LinCalib(inputsfile, incatfile, datapath='', respath='', doBayes=True):
    """ """

    run_ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')
    degRT = 7

    inputs = vjson.load_jsonfile(inputsfile)['inputs']
    SampInter = inputs['SampInter']
    RTlevels = np.array(inputs['RTlevels'])
    Injector = inputs['Injector']
    Date = inputs['Date']

    outexcelfile = os.path.join(
        respath, 'ROETAB_NLCAL_%s_%s.xlsx' % (Injector, Date))

    pixT0 = 14.245E-6  # s
    pixTx0 = int(np.round(pixT0 / SampInter))
    toi0 = pixT0/3.
    # SampInter = 100.E-9 # s
    filt_kernel = int(np.round(toi0/SampInter)/4.)  # waveform filtering kernel

    #print 'filt_kernel = %i' % filt_kernel

    Nlevels = len(RTlevels)

    # Master Inputs File

    indata = ascii.read(incatfile)

    CHANNELS = indata['CHANNEL'].data.copy()
    #FitsList = indata['FITS'].data.copy()
    WFList = indata['WAVEFORM'].data.copy()

    # Loop over CHANNELS to be calibrated

    figs = dict()
    figs['WF'] = dict()
    figs['NL_RT'] = dict()

    CHANNELS = [CHANNELS[0]]  # TESTS

    # Initialisations

    reportdict = OrderedDict()
    reportdict['meta'] = OrderedDict(degRT=degRT,
                                     Analysis_Date=run_ttag)
    reportdict['meta'].update(inputs)

    reportdict['CHANNELS'] = CHANNELS

    reportdict['data'] = OrderedDict()
    reportdict['data']['RT_DN'] = RTlevels
    reportdict['data']['RT_V'] = []
    reportdict['data']['RT_eV'] = []
    reportdict['RT_ADU_2_V'] = []
    reportdict['RT_NONLIN_pol'] = []

    reportdict['RT_NONLIN'] = OrderedDict()
    reportdict['RT_NONLIN']['MaxNL'] = []
    reportdict['RT_NONLIN']['MaxNL_DN'] = []

    for ix, CHAN in enumerate(CHANNELS):

        WFf = os.path.join(datapath, WFList[ix])

        timex, rV = load_WF(WFf, chkNsamp=1.E5, chkSampInter=SampInter)

        fV = filter_Voltage(rV, filt_kernel)

        # EXTRACTING INJECTED VOLTAGE LEVELS

        v_levels = find_discrete_voltages_inwaveform(rV, RTlevels, filtered=fV)
        ev_levels = np.ones_like(v_levels)*v_levels*0.1

        if doBayes:

            bayespick = os.path.join(
                respath, '%s_WFfit_%s_%s_bayes.pick' % (CHAN, Injector, Date))

            # ... DOING BAYESIAN ANALYSIS OF WAVEFORM

            Vrange = [min(v_levels), max(v_levels)]

            bayresults = WaveBay.forwardModel(fV, Vrange, pixTx0, Nlevels, burn=1000, run=2000, cores=8,
                                              doPlot=True, figkey='%s_%s_%s' % (CHAN, Injector, Date), figspath=respath)

            cPickleDumpDictionary(bayresults, bayespick)

            ai = []
            eai = []
            for i in range(1, Nlevels+1):
                ai.append(bayresults['a%i' % i])
                eai.append(bayresults['ea%i' % i])

            ixsort = np.argsort(ai)
            v_levels = np.array(ai)[ixsort]
            ev_levels = np.array(eai)[ixsort]

        relerr_v_levels = ev_levels/v_levels

        reportdict['data']['RT_V'].append(v_levels)
        reportdict['data']['RT_eV'].append(ev_levels)

        pldataset1 = dict(filtered=dict(x=timex[0:5000],
                                        y=fV[0:5000],
                                        ls='-',
                                        color='b'
                                        ))

        figs['WF'][CHAN] = os.path.join(
            respath, '%s_%s_%s_WAVEFORM.png' % (CHAN, Injector, Date))

        NL_lib.myplot_NL(pldataset1, xlabel='t[s]', ylabel='V', ylim=[],
                         title='ROE-TAB Injection. CHAN=%s' % CHAN,
                         figname=figs['WF'][CHAN])

        #  V-DN calibration and NON-LINEARITY of ROE-TAB

        RT_pol_DN2V = NL_lib.fit_pol_byorigin(RTlevels, v_levels,
                                              deg=degRT, sigma=relerr_v_levels)

        RT_pol_NL = NL_lib.find_NL_pol(RTlevels, v_levels,
                                       deg=degRT, sigma=relerr_v_levels)

        RT_pol_LIN = NL_lib.fit_pol_byorigin(RTlevels, v_levels, deg=1,
                                             sigma=relerr_v_levels)

        RT_V_LIN = NL_lib.f_pol_byorigin(RTlevels, *RT_pol_LIN)
        v_RT_data_NL = (v_levels - RT_V_LIN)/RT_V_LIN * 100.

        xRT = np.linspace(RTlevels[0], RTlevels[-1], 1000)

        v_RT_bestfit_NL = NL_lib.f_pol_byorigin(xRT, *RT_pol_NL) * 100.

        pldatasetRTNL = dict(
            data=dict(x=RTlevels, y=v_RT_data_NL, marker='o', color='k'),
            bestfit=dict(x=xRT, y=v_RT_bestfit_NL, ls='--', color='r')
        )

        figs['NL_RT'][CHAN] = os.path.join(
            respath, '%s_%s_%s_NL.png' % (CHAN, Injector, Date))

        NL_lib.myplot_NL(pldatasetRTNL, xlabel='', ylabel='NL[pc]', ylim=[-200., 100.],
                         title='ROE-TAB NonLin. CHAN=%s' % CHAN,
                         figname=figs['NL_RT'][CHAN])

        reportdict['RT_ADU_2_V'].append(RT_pol_DN2V)
        reportdict['RT_NONLIN_pol'].append(RT_pol_NL)

        ixmaxNL = np.argmax(v_RT_bestfit_NL)
        MaxNL = v_RT_bestfit_NL[ixmaxNL]
        MaxNL_DN = xRT[ixmaxNL]

        reportdict['RT_NONLIN']['MaxNL'].append(MaxNL)
        reportdict['RT_NONLIN']['MaxNL_DN'].append(MaxNL_DN)

    reportdict['figs'] = figs

    report = ReportXL_RTLIN(reportdict)

    report.fill_Header()
    report.fill_Data()
    report.fill_RTADU2V()
    report.fill_RTNONLINpols()
    report.fill_RTNONLINMaxima()
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
    parser.add_option("-B", "--Bayes", dest="doBayes",
                      action="store_true", default=False, help="Do Bayesian Analysis?")

    (options, args) = parser.parse_args()

    datapath = options.path
    respath = options.respath
    incat = options.incat
    inputsfile = options.inputs
    doBayes = options.doBayes

    if incat == '' or inputsfile == '':
        parser.print_help()
        sys.exit()

    if respath != '' and not os.path.exists(respath):
        os.system('mkdir %s' % respath)

    header = '\n\n#########################################################\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#       running ROE-TAB NON-LINEARITY ANALYSIS          #\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#########################################################\n'

    print header

    run_ROETAB_LinCalib(inputsfile, incat, datapath=datapath, respath=respath,
                        doBayes=doBayes)
