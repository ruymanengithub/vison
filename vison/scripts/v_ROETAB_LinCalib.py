#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Linearity Calibration of ROE-TAB.

Created on Tue Mar 27 14:42:00 2018
Modified on Fri Sep 14 10:53:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
from itertools import islice
import numpy as np
from astropy.io import ascii
from scipy import ndimage, stats
from sklearn.cluster import KMeans
import os
import copy
import sys
from optparse import OptionParser
from collections import OrderedDict
import pandas as pd
import datetime

from pylab import plot,show
from matplotlib import pyplot as plt

from vison.datamodel import cdp
from vison.support import utils
from vison.roe_fft import NL_lib
from vison.roe_fft import Wave_bayes as WaveBay
from vison.datamodel.ccd import CCD as CCDClass
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import vjson
from vison.support.excel import ReportXL
from vison import __version__ as vison_version
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
        ax.axhline(y=idislev, ls='--', c='r')

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
    fV = ndimage.filters.uniform_filter(rV, filt_kernel)
    return fV


def find_discrete_voltages_inwaveform(rV, levels, filtered=None, debug=False):
    """ """

    if filtered is not None:
        assert len(rV) == len(filtered)
        iV = copy.deepcopy(filtered)
    else:
        iV = copy.deepcopy(rV)

    Nsamp = len(iV)
    Nclust = len(levels)+1

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
    
    if debug:
        fktimex = np.arange(len(filtered))
        Wf = (fktimex[0:50000],filtered[0:50000])
        plot_waveform(Wf, disc_voltages=sorted_levels, 
                      figname='', chan='Unknown')

    return sorted_levels


def run_ROETAB_LinCalib(inputsfile, incatfile, datapath='', respath='', doBayes=False):
    """ """

    run_ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')
    degRT = 7

    inputs = vjson.load_jsonfile(inputsfile)['inputs']
    SampInter = inputs['SampInter']
    RTlevels = np.array(inputs['RTlevels'])
    Injector = inputs['Injector']
    Date = inputs['Date']
    outfileroot = 'ROETAB_NLCAL_%s_%s' % (Injector, Date)
    #outexcelfile = os.path.join(respath, '%s.xlsx' % outfileroot)
    #outpickfile = os.path.join(respath, '%s.pick' % outfileroot)

    pixT0 = 14.245E-6  # s
    pixTx0 = int(np.round(pixT0 / SampInter))
    toi0 = pixT0/3.
    # SampInter = 100.E-9 # s
    filt_kernel = max([10,int(np.round(toi0/SampInter)/32.)])  # waveform filtering kernel

    #print 'filt_kernel = %i' % filt_kernel

    Nlevels = len(RTlevels)

    # Master Inputs File

    indata = ascii.read(incatfile)

    CHANNELS = indata['CHANNEL'].data.copy()
    #FitsList = indata['FITS'].data.copy()
    WFList = indata['WAVEFORM'].data.copy()

    # Loop over CHANNELS to be calibrated

    figs = OrderedDict()

    CHANNELS = np.array([CHANNELS[0]])  # TESTS

    # Initialisations
    
    function, module = utils.get_function_module()
    CDP_header = OrderedDict()
    CDP_header['function'] = function
    CDP_header['module'] = module
    CDP_header['DATE'] = run_ttag
    CDP_header['vison version'] = vison_version
    
    meta = OrderedDict(degRT=degRT)
    meta.update(inputs)
    meta['CHANNELS'] = CHANNELS.tolist().__repr__()

    data = OrderedDict()    
    
    for CHAN in CHANNELS:
        data[CHAN] = OrderedDict()
        data[CHAN]['RT_DN'] = RTlevels
        data[CHAN]['RT_mV'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
        data[CHAN]['RT_emV'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
        data[CHAN]['RT_NL_pc'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
    
    data['RT_ADU_2_mV'] = OrderedDict()
    data['RT_ADU_2_mV']['CHANNEL'] = CHANNELS.copy()
    for i in range(degRT+1):
        data['RT_ADU_2_mV']['P%i' % i] = np.zeros_like(CHANNELS,dtype='float32')+np.nan
    
    
    data['RT_NONLIN_pol'] = OrderedDict()
    data['RT_NONLIN_pol']['CHANNEL'] = CHANNELS.copy()
    for i in range(degRT+1):
        data['RT_NONLIN_pol']['P%i' % i] = np.zeros_like(CHANNELS,dtype='float32')+np.nan
    
    
    data['RT_NONLIN'] = OrderedDict()
    data['RT_NONLIN']['CHANNEL'] = CHANNELS.copy()
    data['RT_NONLIN']['MaxNL_pc'] = np.zeros_like(CHANNELS,dtype='float32')+np.nan
    data['RT_NONLIN']['MaxNL_DN'] = np.zeros_like(CHANNELS,dtype='float32')+np.nan
    
    
    for ix, CHAN in enumerate(CHANNELS):
        
        print 'Processing Channel %s...' % CHAN
        
        WFf = os.path.join(datapath, WFList[ix])
        
        timex, rV = load_WF(WFf, chkNsamp=1.E5, chkSampInter=SampInter)
        rmV = rV * 1000.

        fmV = filter_Voltage(rmV, filt_kernel)

        # EXTRACTING INJECTED VOLTAGE LEVELS

        mv_levels = find_discrete_voltages_inwaveform(rmV, RTlevels, filtered=fmV,
                                                     debug=True) # TESTS
        emv_levels = np.ones_like(mv_levels)*np.abs(mv_levels)*0.1 # 10% uncertainty
        

        if doBayes:

            bayespick = os.path.join(
                respath, '%s_WFfit_%s_%s_bayes.pick' % (CHAN, Injector, Date))

            # ... DOING BAYESIAN ANALYSIS OF WAVEFORM

            mVrange = [min(mv_levels), max(mv_levels)]
            
            print 'Doing Bayesian Analysis...'

            bayresults = WaveBay.forwardModel(fmV, mVrange, pixTx0, Nlevels, burn=1000, run=2000, cores=8,
                                              doPlot=True, figkey='%s_%s_%s' % (CHAN, Injector, Date), 
                                             figspath=respath,debug=True)

            cPickleDumpDictionary(bayresults, bayespick)

            ai = []
            eai = []
            for i in range(1, Nlevels+1):
                ai.append(bayresults['a%i' % i])
                eai.append(bayresults['ea%i' % i])

            ixsort = np.argsort(ai)
            mv_levels = np.array(ai)[ixsort]
            emv_levels = np.array(eai)[ixsort]
            
        
        relerr_mv_levels = emv_levels/mv_levels
        
        # voltage grounding
        mv_levels -= mv_levels.min()

        data[CHAN]['RT_mV'] = mv_levels
        data[CHAN]['RT_emV'] = emv_levels
            
        Np2plot = pixTx0*Nlevels*2

        pldataset1 = dict(filtered=dict(x=timex[0:Np2plot],
                                        y=fmV[0:Np2plot],
                                        ls='-',
                                        color='b'))
        
        figs['WF_%s' % CHAN] = os.path.join(
            respath, '%s_%s_%s_WAVEFORM.png' % (CHAN, Injector, Date))

        NL_lib.myplot_NL(pldataset1, xlabel='t[s]', ylabel='mV', ylim=[],
                         title='ROE-TAB Injection. CHAN=%s' % CHAN,
                         figname=figs['WF_%s' % CHAN])

        #  V-DN calibration and NON-LINEARITY of ROE-TAB
        
        ixgood = np.where(RTlevels>0.)

        RT_pol_NL, mv_RT_data_NL = NL_lib.find_NL_pol(RTlevels[ixgood], mv_levels[ixgood],
                                       deg=degRT, sigma=relerr_mv_levels[ixgood],
                                        Full=True)
        
        data[CHAN]['RT_NL_pc'][ixgood] = mv_RT_data_NL.copy() * 100.
        
        RT_pol_DN2V = np.polyfit(RTlevels,mv_levels,deg=degRT,w=relerr_mv_levels)
        
        xRT = np.linspace(RTlevels[1], RTlevels[-1], 1000)
        
        mv_RT_bestfit_NL = np.polyval(RT_pol_NL,xRT) * 100.
        
        pldatasetRTNL = dict(
            data=dict(x=RTlevels[1:], y=mv_RT_data_NL*100., marker='o', color='k'),
            bestfit=dict(x=xRT, y=mv_RT_bestfit_NL, ls='--', color='r'))
        
        figs['NL_RT_%s' % CHAN] = os.path.join(
            respath, '%s_%s_%s_NL.png' % (CHAN, Injector, Date))
        
        NL_lib.myplot_NL(pldatasetRTNL, xlabel='', ylabel='NL[pc]', ylim=[-10., 10.],
                         title='ROE-TAB NonLin. CHAN=%s' % CHAN,
                         figname=figs['NL_RT_%s' % CHAN])
        
        for ip in range(degRT):
            data['RT_ADU_2_mV']['P%i' % ip][ix] = RT_pol_DN2V[ip]
        
        for ip in range(degRT):
            data['RT_NONLIN_pol']['P%i' % ip][ix] = RT_pol_NL[ip]
        

        ixmaxNL = np.argmax(np.abs(mv_RT_bestfit_NL))
        MaxNL = mv_RT_bestfit_NL[ixmaxNL]
        MaxNL_DN = xRT[ixmaxNL]

        data['RT_NONLIN']['MaxNL_pc'][ix] = MaxNL
        data['RT_NONLIN']['MaxNL_DN'][ix] = MaxNL_DN
            
    figs['keys'] = figs.keys()
    figs['jump'] = 26

    dddf = OrderedDict()
    for key in data.keys():
        dddf[key] = pd.DataFrame.from_dict(data[key])
    
    cdpRTNL = cdp.Tables_CDP()
    cdpRTNL.rootname = outfileroot
    cdpRTNL.path = respath
    
    cdpRTNL.ingest_inputs(
            data = dddf.copy(),
            meta = meta.copy(),
            header  = CDP_header.copy(),
            figs = figs.copy()
            )
    
    cdpRTNL.init_wb_and_fillAll(header_title='ROE-TAB NON LINEARITY Calibration Report')
    cdpRTNL.savetopickle()
    cdpRTNL.savehardcopy()


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

    banner = '\n\n#########################################################\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#       running ROE-TAB NON-LINEARITY ANALYSIS          #\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#########################################################\n'

    print banner
    

    run_ROETAB_LinCalib(inputsfile, incat, datapath=datapath, respath=respath,
                        doBayes=doBayes)
