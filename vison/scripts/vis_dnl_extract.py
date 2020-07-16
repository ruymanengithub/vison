#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Exploring the histograms of counts in images.
Do we have "missing codes" due to larger-than-expected DNL?
DNL: Differential Non Linearity

This issue was brought to attention of VIS IDT by Ralf Kohley on
Jan 10th 2019.


:History:
Created on Tue Feb  5 09:35:25 2019

:author: raf

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from glob import glob
import os
import string as st
from collections import OrderedDict
import copy
from astropy.io import fits as fts
import itertools
from optparse import OptionParser
import sys
from vison.support import vjson
import multiprocessing as mp

from pylab import plot, show

from vison.support.files import cPickleRead, cPickleDumpDictionary
from vison.datamodel import ccd as ccdmod
from vison.pipe import lib as pilib
from vison.support import vistime


import matplotlib.pyplot as plt
plt.switch_backend('TkAgg')
# END IMPORT

Quads = ['E', 'F', 'G', 'H']


def load_all_explogs(parentpath):

    allexplogs = glob(os.path.join(parentpath, '*/EXP_LOG_*.txt'))

    alldatapaths = [os.path.split(item)[0] for item in allexplogs]

    explog = pilib.loadexplogs(allexplogs, elvis='7.5.X',
                               addpedigree=True, datapath=alldatapaths)

    explog['time'] = np.array(
        list(map(vistime.get_dtobj, explog['date']))).copy()

    explog.sort('time')

    return explog


def explog_selector(explog, testkeys):
    """ """

    N = len(explog)

    selbool = np.zeros(N, dtype='bool')

    for i in range(N):
        testname = explog['test'][i]
        selbool[i] = np.any([tkey in testname for tkey in testkeys])

    selix = np.where(selbool)

    explog = explog[selix]

    return explog


def extract_histograms(infile_name, i, N):
    """ """

    print(('Extracting histos from image %i/%i' % (i + 1, N)))

    ccdobj = ccdmod.CCD(infile_name)

    outfile_name = os.path.join(respath, '%s_H.pick' % infile_name)

    bins = np.arange(2**16)

    histos = OrderedDict()
    histos['bins'] = bins

    for Q in Quads:

        quaddata = ccdobj.get_quad(Q, canonical=True, extension=-1)

        Qhist, _ = np.histogram(quaddata, bins=bins)

        histos[Q] = Qhist

    cPickleDumpDictionary(histos, outfile_name)


def _wrap_extract_histograms(args):
    extract_histograms(*args)


def batch_extract_histos(explog, respath, multithread=1):
    """ """

    N = len(explog)

    arglist = []

    for i in range(N):

        datapath = explog['datapath'][i]
        infile_name = explog['File_name'][i]
        fits_name = os.path.join(datapath, '%s.fits' % infile_name)

        if os.path.exists(fits_name):
            arglist.append((fits_name, i, N))
        else:
            print('%s is MISSING' % fits_name)

    pool = mp.Pool(processes=multithread)
    pool.map(_wrap_extract_histograms, arglist)
    pool.close()


def interpol_histo_linear(histo):
    """ """

    #N = len(histo)
    thresh = 100  # code counts

    ihisto = np.zeros_like(histo) + np.nan
    #ahisto = np.zeros_like(histo)+np.nan

    def _f_interp(x, y):
        if np.any(y < thresh):
            return np.nan
        else:
            return x * (y[1] - y[0]) / 2. + y[0]

    # for i in range(1,N-1):
    #    y = np.array([histo[i-1],histo[i+1]])
    #    ahisto[i] = _f_interp(1.,y)

    ihisto[1:-1] = (histo[2:] - histo[0:-2]) / 2. + histo[0:-2]
    nonvalid = np.where(((histo[2:] < thresh) | (histo[0:-2] < thresh)))
    ihisto[1:-1][nonvalid] = np.nan

    try:
        ihisto[0] = _f_interp(-1., np.array([histo[1], histo[2]]))
    except BaseException:
        stop()
    try:
        ihisto[-1] = _f_interp(+3., np.array([histo[-3], histo[-2]]))
    except BaseException:
        stop()
    #ahisto[0] = _f_interp(-1.,np.array(histo[1],histo[2]))
    #ahisto[-1] = _f_interp(+3.,np.array(histo[-3],histo[-2]))

    return ihisto


def interpol_histo_cubic(histo):
    """NOT WORKING PROPERLY, DO NOT USE"""

    #N = len(histo)
    thresh = 100  # code counts

    N = len(histo)

    ihisto = np.zeros(N, dtype='float32') + np.nan

    x = [0, 1, 3, 4]
    for i in range(2, N - 2):
        y = [histo[i - 2], histo[i - 1], histo[i + 1], histo[i + 2]]
        if np.any(y < thresh):
            ihisto[i] = np.nan
            continue
        p = np.polyfit(x, y, 3)
        ihisto[i] = np.polyval(p, 2)

    return ihisto


def extract_dnls(histofile_name, dnl_name, i, N, interptype):

    print(('Extracting DNLs from image %i/%i' % (i + 1, N)))

    histodict = cPickleRead(histofile_name)

    bins = histodict['bins']

    dnldict = OrderedDict()
    dnldict['bins'] = bins.copy()

    if interptype == 'linear':
        interpoler = interpol_histo_linear
    elif interptype == 'cubic':
        interpoler = interpol_histo_cubic

    for Q in Quads:

        dnldict[Q] = OrderedDict()

        qh = histodict[Q].astype('float32').copy()

        dnldict[Q]['histo'] = qh.copy()

        interp_qh = interpoler(qh)

        qdnl = qh / interp_qh - 1.

        dnldict[Q]['dnl'] = qdnl.copy()

    cPickleDumpDictionary(dnldict, dnl_name)


def _wrap_extract_dnls(args):
    extract_dnls(*args)


def batch_extract_dnls(explog, respath, interptype, multithread=1):
    """ """

    N = len(explog)

    arglist = []

    for i in range(N):

        infile_name = explog['File_name'][i]
        histofile_name = os.path.join(respath, '%s_H.pick' % infile_name)
        dnl_name = os.path.join(respath, '%s_DNL.pick' % infile_name)

        if os.path.exists(histofile_name):
            arglist.append((histofile_name, dnl_name, i, N, interptype))
        else:

            print(('%s is MISSING' % histofile_name))

    pool = mp.Pool(processes=multithread)
    pool.map(_wrap_extract_dnls, arglist)
    pool.close()


def save_to_fits(array, fitsfile):
    """ """
    fts.writeto(fitsfile, array, overwrite=True)


def get_average_dnls(explog, respath, tag):
    """ """

    CCDs = [1, 2, 3]

    avgdnldict = OrderedDict()

    for CCD in CCDs:
        CCDk = 'CCD%i' % CCD
        avgdnldict[CCDk] = OrderedDict()
        for Q in Quads:
            avgdnldict[CCDk][Q] = OrderedDict()

    for CCD in CCDs:

        CCDk = 'CCD%i' % CCD

        selCCD = np.where(explog['CCD'] == CCDk)

        _explog = copy.deepcopy(explog[selCCD])

        N = len(_explog)

        dnlarraydict = OrderedDict()
        for Q in Quads:
            dnlarraydict[Q] = np.zeros((N, 2**16 - 1), dtype='float32') + np.nan

        for i in range(N):

            infile_name = _explog['File_name'][i]
            dnl_name = os.path.join(respath, '%s_DNL.pick' % infile_name)

            dnldict = cPickleRead(dnl_name)
            if i == 0:
                bins = dnldict['bins'][0:-1].copy()

            for Q in Quads:

                dnlarraydict[Q][i, :] = dnldict[Q]['dnl'].copy()

        for Q in Quads:

            fitsfile = os.path.join(respath, '%s_codes_%s_%s.fits' % (tag, CCDk, Q))

            save_to_fits(dnlarraydict[Q], fitsfile)

        for Q in Quads:

            avgdnldict[CCDk][Q]['bins'] = bins
            avgdnldict[CCDk][Q]['avgdnl'] = np.nanmean(dnlarraydict[Q], axis=0)

    return avgdnldict


def show_average_DNLs(data, filename='', ylim=[-2., 2.],
                      xlim=[0., 2**16 - 1]):

    CCDs = [1, 2, 3]

    fig, axsarr = plt.subplots(
        2, 6, sharex=True, sharey=True, figsize=(12, 8))

    axs = dict()

    # initialisation of self.axs

    for CCD in CCDs:
        CCDkey = 'CCD%i' % CCD
        axs[CCDkey] = dict()
        for Q in Quads:
            axs[CCDkey][Q] = None

    axs['CCD1']['E'] = axsarr[0, 0]

    plotlist = [item for item in itertools.product(CCDs, ['E', 'F'])] +\
               [item for item in itertools.product(CCDs, ['H', 'G'])]

    for k in range(1, len(plotlist) + 1):
        CCDkey = 'CCD%i' % plotlist[k - 1][0]
        Q = plotlist[k - 1][1]
        axs[CCDkey][Q] = axsarr.flatten()[k - 1]

    for CCD in CCDs:
        CCDkey = 'CCD%i' % CCD
        for Q in Quads:

            ax = axs[CCDkey][Q]
            CQdict = data[CCDkey][Q]

            ax.plot(CQdict['bins'], CQdict['avgdnl'], marker='.', linestyle='', color='k', ms=1,
                    alpha=0.3)

            if Q in ['E', 'H']:
                ax.text(0.05, 0.9, Q, horizontalalignment='left',
                        transform=axs[CCDkey][Q].transAxes)
            elif Q in ['F', 'G']:
                ax.text(0.9, 0.9, Q, horizontalalignment='right',
                        transform=axs[CCDkey][Q].transAxes)

            if Q == 'E':
                ax.set_title(CCDkey, x=1)

            if Q in ['H', 'G']:
                ax.set_xlabel('code')
            if Q in ['E', 'H'] and CCD == 1:
                ax.set_ylabel(r'$<DNL>$')

            ax.set_ylim(ylim)
            ax.set_xlim(xlim)

    for CCD in CCDs:
        for Q in ['E', 'F']:
            plt.setp(axs['CCD%i' % CCD]
                     [Q].get_xticklabels(), visible=False)
        if CCD != 1:
            for Q in Quads:
                plt.setp(axs['CCD%i' % CCD]
                         [Q].get_yticklabels(), visible=False)
        else:
            for Q in ['F', 'G']:
                plt.setp(axs['CCD%i' % CCD]
                         [Q].get_yticklabels(), visible=False)

    plt.subplots_adjust(hspace=0.0)
    plt.subplots_adjust(wspace=0.0)

    plt.margins(0.05)

    plt.subplots_adjust(top=0.85)

    if filename != '':
        plt.savefig(filename)
    else:
        plt.show()


def run(inpath, respath, figspath, doHistos, doDNLs, doAvgDNLs,
        tag, testkeys, interptype, multithread=1):

    if not os.path.exists(respath):
        os.system('mkdir %s' % respath)
    if not os.path.exists(figspath):
        os.system('mkdir %s' % figspath)

    avgdnlspick = os.path.join(respath, 'AVG_DNLs_%s.pick' % tag)
    avgdnlfig = os.path.join(figspath, 'AVG_DNLs_%s.png' % tag)

    explog = load_all_explogs(parentpath)
    explog = explog_selector(explog, testkeys)

    if doHistos:
        batch_extract_histos(explog, respath, multithread=multithread)

    if doDNLs:
        batch_extract_dnls(explog, respath, interptype, multithread=multithread)

    if doAvgDNLs:
        avgdnlsdict = get_average_dnls(explog, respath, tag)
        cPickleDumpDictionary(avgdnlsdict, avgdnlspick)

    avgdnlsdict = cPickleRead(avgdnlspick)

    show_average_DNLs(avgdnlsdict, filename=avgdnlfig, ylim=[-2, 2],
                      xlim=[0, 2**16])


if __name__ == '__main__':

    #parentpath = '../atCALDATA/data'
    #respath = 'results_atCALDATA/misscodesresults'
    #figspath = 'results_atCALDATA/misscodesfigs'
    #doExtractHistos = True
    #doEstimateDNLs = True
    #doAverageDNLs = True
    #tag = 'BORN'
    #testkeys = ['BIAS','CHINJ','FLAT','NL','PTC']

    parser = OptionParser()
    parser.add_option("-j", "--json", dest="json",
                      default='', help="json file with inputs")
    parser.add_option("-H", "--doHistos", dest="H", action="store_true",
                      default=False, help="Do Histograms of codes for each image.")
    parser.add_option("-D", "--doDNLs", dest="DNL", action="store_true",
                      default=False, help="Do DNLs for each image.")
    parser.add_option("-A", "--doAverageDNLs", dest="AVGDNL", action="store_true",
                      default=False, help="Do Average DNL for each quadrant.")
    parser.add_option("-t", "--tag", dest="tag", default="",
                      help="Tag to be added to file names for ease of identification. Optional.")
    parser.add_option(
        "-m",
        "--multithread",
        dest="multithread",
        default=1,
        help="Use multithreading? Number of threads / cores to use. Default=1 [single thread]. Number of threads must be " +
        "< number of available cores. Only some tasks are parallelized.")
    parser.add_option(
        "-i",
        "--inter",
        dest="interptype",
        default="linear",
        help="Histogram interpolation type: linear/cubic. ONLY LINEAR AVAILABLE BY NOW.")

    (options, args) = parser.parse_args()

    if options.json == '':
        parser.print_help()
        sys.exit()

    verbose_inputs = vjson.load_jsonfile(options.json, useyaml=True)

    parentpath = verbose_inputs['parentpath']
    respath = verbose_inputs['respath']
    figspath = verbose_inputs['figspath']
    testkeys = verbose_inputs['testkeys']

    doHistos = options.H
    doDNLs = options.DNL
    doAvgDNLs = options.AVGDNL
    tag = options.tag
    multithread = int(options.multithread)
    interptype = options.interptype

    assert interptype == 'linear'

    run(parentpath, respath, figspath, doHistos, doDNLs,
        doAvgDNLs, tag, testkeys, interptype, multithread)
