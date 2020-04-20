#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

ROE-TAB Calibration
Module to fit the Non-Linearity waveform, bayesian style.

Created on Mon Mar 19 13:28:52 2018

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from multiprocessing import Pool
import emcee
import corner
import os

import matplotlib
# matplotlib.use('Agg')
#matplotlib.rc('text', usetex=True)
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['image.interpolation'] = 'none'
import matplotlib.pyplot as plt

from pylab import plot, show
# END IMPORT


def get_priors(Vrange, pixTx, Nlevels):
    """ """
    # theta: phase, ground, wref, wpix, a_i {i=1,N}
    vmin, vmax = Vrange
    avmin = np.abs(vmin)

    #priors = [[-(Nlevels)*pixTx/2,0],[-2.*avmin,2.*avmin],[pixTx/4,3*pixTx/4],[0.8*pixTx,1.2*pixTx]]
    priors = [[-0.5 * pixTx, 0.5 * pixTx], [-1.2 * avmin, 1.2 * avmin],
              [pixTx / 4, 3 * pixTx / 4], [0.95 * pixTx, 1.05 * pixTx]]
    for i in range(Nlevels - 1, -1, -1):
        #expectV = (vmax-vmin)/Nlevels * i
        #priors += [[max(0,expectV-0.2),expectV+0.2]]
        priors += [[0., 1.2 * (vmax - vmin)]]

    return priors

#==============================================================================
# def waveform_generator_old(theta,N):
#     """ """
#     # theta: phase, ground, wref, wpix, a_i {i=1,N}
#     Nlevels = len(theta)-4
#
#     phi, gnd, wref, wpix = theta[0:4]
#     phi = int(phi)
#     wref = int(wref)
#     wpix = int(wpix)
#
#     ai = theta[4:]
#     Nlevels = len(ai)
#
#     Npercycle = wpix * Nlevels
#
#     cycle = np.zeros(Npercycle,dtype='float32') + gnd
#
#     for i in range(Nlevels):
#         cycle[i*(wpix)+wref:i*(wpix)+wref+(wpix-wref)] += ai[i]
#
#     indices = (np.arange(N)+phi) % Npercycle
#
#     waveform = cycle[indices]
#
#
#     return waveform
#
#==============================================================================


def waveform_generator(theta, N):
    """ """
    # theta: phase, ground, wref, wpix, a_i {i=1,N}
    Nlevels = len(theta) - 4

    phi, gnd, wref, wpix = theta[0:4]

    ai = np.array(theta[4:])
    Nlevels = len(ai)

    Npercycle = wpix * Nlevels

    fraction_pixsignal = wref / float(wpix)

    stepcount = np.arange(N, dtype='float32')

    cyclecount = (stepcount + phi) / (Npercycle)

    cyclephase = np.modf(cyclecount)[0] * Nlevels
    modf_cyclephase = np.modf(cyclephase)

    pixelid = modf_cyclephase[1].astype('int32')
    #ixneg = np.where(pixelid < 0)
    pixfrac = modf_cyclephase[0]

    negfrac = np.where((pixfrac < 0) | (pixelid < 0))

    pixelid[negfrac] = Nlevels - 1 + pixelid[negfrac]

    pixfrac[negfrac] = 1. + pixfrac[negfrac]

    hassignal = (np.abs(pixfrac) > fraction_pixsignal).astype('int32')

    waveform = hassignal * ai[(pixelid,)]

    waveform += gnd

    return waveform


def log_prior(theta, priors):
    """
    Generic model

    Priors, limit the values to a range but otherwise flat.
    """

    areinrange = [(priors[i][0] <= theta[i] <= priors[i][1]) for i in range(len(theta))]
    disc = np.prod(areinrange)

    if disc:
        return 0.
    else:
        return -np.inf


def log_posterior(theta, Varr, var, priors):
    """
    Posterior probability: combines the prior and likelihood.
    """

    lp = log_prior(theta, priors)

    if not np.isfinite(lp):
        return -np.inf

    return lp + log_likelihood(theta, Varr, var)


def log_likelihood(theta, Varr, var):
    """
    Logarithm of the likelihood function.
    """

    model = waveform_generator(theta, len(Varr))
    plot(model)
    show()
    stop()

    lnL = - 0.5 * np.sum((Varr - model)**2. / var)
    #print lnL
    return lnL


def forwardModel(Varr, Vrange, pixTx, Nlevels, burn=500, run=700, cores=8, figkey='', doPlot=False,
                 figspath='', debug=False):
    """
    Forward models the Non-Linearity waveform: extracts reliable estimates of level voltages.

    Notes:
    - emcee is run three times as it is important to have a good starting point for the final run.

    """
    var = 0.04**2.
    nwalkers = 1000

    # theta: phase, ground, wref, wpix, a_i {i=1,N}
    parnames = ['phi', 'gnd', 'wref', 'wpix']
    for i in range(1, Nlevels + 1):
        parnames += ['a%i' % i]

    priors = get_priors(Vrange, pixTx, Nlevels)

    results = dict()

    print("Let's go Bayes...")

    ndim = len(parnames)

    # Choose an initial set of positions for the walkers - fairly large area not to bias the results
    p0 = np.zeros((nwalkers, ndim))

    for ip in range(ndim):
        try:
            p0[:, ip] = np.random.uniform(priors[ip][0], priors[ip][1], nwalkers)
        except BaseException:
            stop()

    #p0[:,parnames.index('wpix')] = pixTx - p0[:,parnames.index('wref')]

    # for intpar in ['phi','wref','wpix']:
    #    ix = parnames.index(intpar)
    #    p0[:,ix] = p0[:,ix].astype('int32')

    if cores > 1:
        # A hack Dan gave me to not have ghost processes running as with threads keyword
        pool = Pool(cores)
        # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[xx,
        # yy, data, var, peakrange, spot.shape],
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                                        args=[Varr, var, priors],
                                        pool=pool)
    else:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                                        args=[Varr, var, priors])

    # Run a burn-in and set new starting position
    print("Burning-in...")
    pos, prob, state = sampler.run_mcmc(p0, burn)
    maxprob_index = np.argmax(prob)
    params_fit = pos[maxprob_index]
    print(("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction)))
    print(('Estimate:', params_fit))
    sampler.reset()

    print("Running MCMC...")
    pos, prob, state = sampler.run_mcmc(pos, run, rstate0=state)
    print(("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction)))

    # Get the index with the highest probability
    maxprob_index = np.argmax(prob)

    if cores > 1:
        pool.close()

    # Get the best parameters and their respective errors and print best fits
    params_fit = pos[maxprob_index]
    errors_fit = [sampler.flatchain[:, i].std() for i in range(ndim)]

    #_printResults(params_fit, errors_fit)

    # Best fit model
    # theta: phase, ground, wref, wpix, a_i {i=1,N}

    model = waveform_generator(params_fit, len(Varr))

    if doPlot:

        if debug:
            figname1 = ''
        else:
            figname1 = os.path.join(figspath, '%s_WaveformBayes.png' % figkey)

        maxix = min(50000, len(Varr))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(Varr[0:maxix], 'b-', label='data')
        ax.plot(model[0:maxix], 'k--', label='model')
        ax.set_ylabel('V')
        ax.set_title('Bayesian Fit of Waveform: %s' % figkey)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc='best')

        plt.tight_layout()

        plt.savefig(figname1)
        plt.close()

    # packing results

    for ip, parname in enumerate(parnames):
        results[parname] = params_fit[ip]
        results['e%s' % parname] = errors_fit[ip]

    #results['GoF'] = gof

    # plot
    samples = sampler.chain.reshape((-1, ndim))

    if doPlot:

        figname2 = os.path.join(figspath, '%s_NLwaveform_Triangle.png' % figkey)
        fig2 = corner.corner(samples, labels=parnames)
        fig2.suptitle('Bayesian Fit: Non-Lin Waveform (%s)' % figkey)
        fig2.savefig(figname2)
        plt.close()

    return results


def test0():
    """ """

    import ROE_LinCalib as RLC
    import os
    from astropy.io import ascii

    Nlevels = 17
    pixT = 14.245E-6  # s
    SampInter = 100.E-9
    datapath = 'CH1_Top_E_F_Pix_bounce/With_17_levels/'
    datafile = 'LinCalib_FQM_ROE_15Mar18.txt'

    pixTx = int(np.round(pixT / SampInter))
    toi = pixT / 3.
    # SampInter = 100.E-9 # s

    filt_kernel = int(np.round(toi / SampInter) / 4.)

    indata = ascii.read(datafile)

    WFList = indata['WAVEFORM'].data.copy()

    WFf = os.path.join(datapath, WFList[0])

    timex, rV = RLC.load_WF(WFf, chkNsamp=1.E5, chkSampInter=SampInter)

    fV = RLC.filter_Voltage(rV, filt_kernel)

    vmin, vmax = (fV.min(), fV.max())
    theta = [-600, -0.33, 68., 142.25]
    for i in range(Nlevels - 1, -1, -1):
        iexpectV = (vmax - vmin) / Nlevels * i
        theta += [iexpectV]

    modV = waveform_generator(theta, len(fV))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(fV)
    ax.plot(modV, 'r--')
    ax.set_xlim([0, 7000])

    plt.show()


if __name__ == '__main__':

    test0()
