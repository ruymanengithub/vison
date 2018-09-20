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
import pandas as pd

from vison.roe_fft import NL_lib
from vison.datamodel.ccd import CCD as CCDClass
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import vjson
#from vison.support.excel import ReportXL
from vison.datamodel import cdp
from vison.support import utils
from vison import __version__ as vison_version

from pylab import plot,show
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



def run_ROE_LinCalib(inputsfile, incatfile, datapath='', respath='', doExtractFits=True):
    """ """

    run_ttag = (datetime.datetime.now()).strftime('%d%b%y_%H%M%S')
    degROE = 10

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

    #CHANNELS = np.array([CHANNELS[0]])  # TESTS

    # Initialisations
    
    function, module = utils.get_function_module()
    CDP_header = OrderedDict()
    CDP_header['function'] = function
    CDP_header['module'] = module
    CDP_header['DATE'] = run_ttag
    CDP_header['vison version'] = vison_version

    outfile = os.path.join('ROE_NLCAL_%s_%s' % (ROE, Date))

    adu_levels_dict = OrderedDict()

    meta = OrderedDict(
            degROE=degROE)
    meta.update(inputs)
    meta['CHANNELS'] = CHANNELS.tolist().__repr__()
    
    figs = OrderedDict(keys=CHANNELS,
                       jump=26)
    
    data = OrderedDict()
    
    data['ROE_NONLIN'] = OrderedDict()
    data['ROE_NONLIN']['CHANNEL'] = CHANNELS.copy()
    data['ROE_NONLIN']['MaxNL_pc'] = []
    data['ROE_NONLIN']['MaxNL_mV'] = []
    
    data['ROE_NONLIN_POLS'] = OrderedDict()
    data['ROE_NONLIN_POLS']['CHANNEL'] = CHANNELS.copy()
    for i in range(degROE+1):
        data['ROE_NONLIN_POLS']['P%i' % i] = []
    
    
    for CHAN in CHANNELS:
        data[CHAN] = OrderedDict()
        data[CHAN]['RT_DN'] = RTlevels
        data[CHAN]['RT_mV'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
        data[CHAN]['ADC_ADUS'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
        data[CHAN]['fit_LS'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
        data[CHAN]['NLpc_LS'] = np.zeros_like(RTlevels,dtype='float32')+np.nan
        
    # Main Loop
    
    injCHANs = InjectorCal['CHANNELS']
    
    for ix, CHAN in enumerate(CHANNELS):
        
        print 'Processing Channel %s...' % CHAN
        
        adus_pickf = os.path.join(respath, '%s_%s_%s_NL_ADUs.pick' % (ROE, Date, CHAN))

        Q = CHAN[1]

        ixRT = injCHANs.index(CHAN)

        RT_pol_DN2V = InjectorCal['RT_ADU_2_V'][ixRT]

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
            
        #Vfitlevels = NL_lib.f_pol_byorigin(RTlevels, *RT_pol_V2DN)
        
        mVfitlevels = np.polyval(RT_pol_DN2V,RTlevels) * 1000.

        # NON-LINEARITY OF ROE

        ixgood = np.where((adu_levels < 2.**16-1.) & (mVfitlevels>10.))

        #adu_levels -= adu_levels[0]  # BIAS subtraction
        
        data[CHAN]['RT_mV'] = mVfitlevels
        data[CHAN]['ADC_ADUS'] = adu_levels

        R_pol_NL, R_data_NL = NL_lib.find_NL_pol(mVfitlevels[ixgood], 
                adu_levels[ixgood], deg=degROE, Full=True, debug=False)
        
        xR = np.linspace(mVfitlevels[ixgood][0], mVfitlevels[ixgood][-1], 1000)
        
        R_bestfit_NL = np.polyval(R_pol_NL,xR) * 100.
        
        #data[CHAN]['fit_LS'][ixgood] = ??
        data[CHAN]['NLpc_LS'][ixgood] = R_data_NL*100.
        
        for ip in range(degROE+1):
            data['ROE_NONLIN_POLS']['P%i' % ip].append(R_pol_NL[ip])

        ixmaxNL = np.argmax(np.abs(R_bestfit_NL))
        MaxNL = R_bestfit_NL[ixmaxNL]
        MaxNL_mV = xR[ixmaxNL]

        data['ROE_NONLIN']['MaxNL_pc'].append(MaxNL)
        data['ROE_NONLIN']['MaxNL_mV'].append(MaxNL_mV)


        pldatasetRNL = dict(
            data=dict(x=mVfitlevels[ixgood], y=R_data_NL*100., marker='o', color='k'),
            bestfit=dict(x=xR, y=R_bestfit_NL, ls='--', color='r')
        )

        figs[CHAN] = os.path.join(
            respath, '%s_%s_%s_ROE_NL.png' % (CHAN, ROE, Date))

        NL_lib.myplot_NL(pldatasetRNL, xlabel='mV', ylabel='NL[pc]', #ylim=[-10., 10.],
                         title='ROE NonLin. CHAN=%s' % CHAN,
                         figname=figs[CHAN])# TESTS

    
    dddf = OrderedDict()
    for key in data.keys():
        dddf[key] = pd.DataFrame.from_dict(data[key])
              
    cdpRNL = cdp.Tables_CDP()
    cdpRNL.rootname = outfile
    cdpRNL.path = respath
    stop()
    
    cdpRNL.ingest_inputs(
            data = dddf.copy(),
            meta = meta.copy(),
            header  = CDP_header.copy(),
            figs = figs.copy()
            )
    
    cdpRNL.init_wb_and_fillAll(header_title='ROE NON LINEARITY Calibration Report')
    cdpRNL.savetopickle()
    cdpRNL.savehardcopy()
    



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
