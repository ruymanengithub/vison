#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in PTC analysis.

Created on Thu Sep 14 16:29:36 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from collections import OrderedDict
#from pylab import plot,show

from vison.support import flags as flmod
# END IMPORT


fitPTC_flags = OrderedDict()
fitPTC_flags['EXCEPTION'] = [2**0]
fitPTC_flags['POORFIT'] = [2**1]
fitPTC_flags['BADERRORS'] = [2**2]


def _prune_smear(means, var):
    """ """
    mask = np.ones_like(means, dtype='bool')
    bins = np.linspace(2**16 * 0.5, 2**16, 50)
    threshold = 0.1

    for i in range(len(bins) - 1):
        sel = np.where((means >= bins[i]) & (means < bins[i + 1]))
        if len(sel[0]) > 0:
            relrange = (var[sel].max() - var[sel].min()) / var[sel].mean()
            if relrange > threshold:
                mask[sel] = False
    return mask


def fitPTC(means, var, debug=False):
    """Fits Photon Transfer Curve to obtain gain."""
    from sklearn.linear_model import RANSACRegressor
    from sklearn.metrics import mean_squared_error
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.pipeline import make_pipeline

    if debug:
        from matplotlib import pyplot as plt

    sigmathresh = 5.  # sigma clipping
    # flulims = [1.E3,2**16-1.E3] # FLUENCE LIMITS
    flulims = [1.e3, 4.E4]  # FLUENCE LIMITS
    poldeg = 2
    flags = flmod.Flags(fitPTC_flags)

    order = np.argsort(means)
    means = np.array(means)[order]
    var = np.array(var)[order]

    smearmask = _prune_smear(means, var)
    presel = smearmask & (means >= flulims[0]) & (means <= flulims[1])

    means = means[np.where(presel)]
    var = var[np.where(presel)]

    #ixmaxvar = np.argmax(var)
    #maxvar = var[ixmaxvar]
    #ixsel = np.where((means < means[ixmaxvar]) & (var < maxvar*0.95))

    try:

        pipe = make_pipeline(PolynomialFeatures(poldeg, interaction_only=True),
                             RANSACRegressor())
        pipe.fit(np.expand_dims(means, 1), np.expand_dims(var, 1))
        robpredict = np.squeeze(pipe.predict(np.expand_dims(means, 1)))

        ixsel = np.where((np.abs(var - robpredict) / np.sqrt(robpredict) < sigmathresh))

        res = np.polyfit(means[ixsel], var[ixsel], poldeg, full=False, cov=True)

        p = res[0]
        V = res[1]

    except BaseException:
        p = np.zeros(poldeg + 1)
        p[0] = 0.01
        V = np.zeros((poldeg + 1, poldeg + 1), dtype='float32')

        flags.add('EXCEPTION')

    ep = np.sqrt(np.diag(V))

    if np.any((ep == 0.) | np.isinf(ep) | np.isnan(ep)):
        flags.add('BADERRORS')

    gain = 1. / p[1]

    if ep[1] / gain > 1.e-3:
        flags.add('POORFIT')

    quality = flags.value
    fitresults = dict(fit=p, efit=ep, gain=gain,
                      quadterm=p[0], rn=p[2],
                      quality=quality)

    if debug:
        fig = plt.figure(figsize=(8, 4))
        ax1 = fig.add_subplot(121)
        ax1.plot(means, var, 'r.')
        ax1.plot(means[ixsel], var[ixsel], 'b.')
        ax1.plot(means, np.polyval(p, means), 'k-')
        ax2 = fig.add_subplot(122)
        ax2.plot(means, var - np.polyval(p, means), 'k.')
        plt.show()
        stop()

    return fitresults


def foo_bloom_advanced_demoted(means, var, _fit, debug=False):
    """
    Finds blooming limit (where variance drops, if it does...).

    """

    bloom_ADU = np.nan

    var_bf = np.polyval(_fit['fit'], means)

    var_res = var - var_bf

    withinconfidence = np.where((means > 2.E3) & (means < 3.E4))
    var_mad = np.median(np.abs(var_res[withinconfidence]))  # median absolute deviation

    Nbins = 40
    bins = np.linspace(2**16 * 0.1, 2**16, Nbins)
    thresholdfactor = 15.

    #avgpopbin = len(means[means>=bins[0]])/float(Nbins)

    mask = np.ones(Nbins)
    binmeans = np.zeros(Nbins)
    binstds = np.zeros(Nbins)

    for i in range(Nbins - 1):
        sel = np.where((means >= bins[i]) & (means < bins[i + 1]))
        binmeans[i] = means[sel].mean()

        if len(sel[0]) > 9:
            #local_mad = np.median(np.abs(var[sel]-np.median(var[sel])))
            local_resstd = np.median(var[sel] - var_bf[sel])
            binstds[i] = local_resstd
            if local_resstd < -thresholdfactor * var_mad:
                mask[i] = 0

    if np.any(mask == 0.):
        selector = mask[1:] + mask[0:-1]

        try:
            bloom_ADU = binmeans[np.where(selector == 0)[0][0]]
        except BaseException:
            bloom_ADU = -means[-1]
    else:
        bloom_ADU = -means[-1]

    res = dict(bloom_ADU=bloom_ADU)

    if debug:

        from matplotlib import pyplot as plt
        fig = plt.figure(figsize=(12, 4))
        ax1 = fig.add_subplot(131)
        ax1.plot(binmeans, binstds, 'k.')
        ax1.axhline(y=-thresholdfactor * var_mad, ls='--', color='r')
        ax1.axvline(x=bloom_ADU, ls='-', color='g')
        ax2 = fig.add_subplot(132)
        ax2.plot(means, var, 'b.')
        ax2.axvline(x=bloom_ADU, ls='-', color='g')
        ax2.set_ylim([0, 1.2E4])
        ax3 = fig.add_subplot(133)
        ax3.plot(means, var_res, 'b.')
        ax3.axvline(x=bloom_ADU, ls='-', color='g')
        plt.show()
        stop()

    return res


def foo_bloom_advanced(means, var, _fit, debug=False):
    """
    Finds blooming limit (where variance drops, if it does...).

    """

    bloom_ADU = np.nan

    var_bf = np.polyval(_fit['fit'], means)

    var_res = var - var_bf

    withinconfidence = np.where((means > 2.E3) & (means < 3.E4))
    var_mad = np.median(np.abs(var_res[withinconfidence]))  # median absolute deviation,
    # a robust estimate of std
    thresholdfactor = 10.
    thresholdval = 3.E4

    mask = var_res < -thresholdfactor * var_mad

    bloom_ADU = -means[-1]
    for i in range(len(means) - 5):
        #print('%i %s %s' % (means[i],mask[i:i+3].__repr__(),np.all(mask[i:i+4])))
        if means[i] > thresholdval and np.all(mask[i:i + 5]):
            bloom_ADU = means[i]
            break

    res = dict(bloom_ADU=bloom_ADU)

    if debug:

        from matplotlib import pyplot as plt
        fig = plt.figure(figsize=(12, 4))
        ax1 = fig.add_subplot(121)
        ax1.plot(means, var, 'b.')
        ax1.axvline(x=bloom_ADU, ls='-', color='g')
        ax1.set_ylim([0, 1.3E4])
        ax2 = fig.add_subplot(122)
        ax2.plot(means, var_res, 'b.')
        ax2.axhline(y=-thresholdfactor * var_mad, ls='--', color='r')
        ax2.axhline(y=thresholdfactor * var_mad, ls='--', color='b')
        ax2.axvline(x=bloom_ADU, ls='-', color='g')
        plt.show()
        stop()

    return res
