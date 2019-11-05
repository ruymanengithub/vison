# -*- coding: utf-8 -*-
"""
Created on Tue May  2 15:33:17 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import warnings
from scipy import interpolate

from vison.pipe import lib as pilib
from vison.point import lib as polib

import matplotlib
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
import matplotlib.pyplot as plt

# END IMPORT


def fit_focus_single(x, y, yerror=None, degree=1, doplot=False):
    """ """
    assert len(x) == len(y)

    if yerror is not None:
        assert len(y) == len(yerror)
        weights = 1. / yerror
    else:
        weights = np.ones_like(y)

    warnings.simplefilter('ignore', np.RankWarning)
    coeffs, Vcoeffs = np.polyfit(x, y, degree, w=weights, full=False, cov=True)
    # try:
    #    coeffs, Vcoeffs = np.polyfit(x, y, degree, w=weights, full=False, cov=True)
    # except ValueError:
    #    res = dict(coeffs=np.zeros(degree+1)+np.nan,
    #               ecoeffs=np.zeros(degree+1)+np.nan,
    #               focus=np.nan)
    #    return res

    ecoeffs = np.diag(Vcoeffs)**0.5

    pol = np.poly1d(coeffs)

    xfocus = -coeffs[1] / (2. * coeffs[0])
    # exfocus = np.abs(xfocus)*np.sqrt((ecoeffs[1]/coeffs[1])**2.+
    #     1/2.*(ecoeffs[0]/coeffs[0])**2.) # WRONG?

    if doplot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if yerror is None:
            ax.plot(x, y, 'bo')
        else:
            ax.errorbar(x, y, yerr=yerror, fmt='o', color='b')
        ax.plot(x, pol(x), 'r--')
        ax.axvline(x=xfocus, ls='--', color='k')
        ax.set_xlabel('Mirr_pos')
        ax.set_ylabel('FWHM')

        plt.show()

    res = dict(coeffs=coeffs, ecoeffs=ecoeffs, focus=xfocus)

    return res


def build_fwhm_map_CQ(delta_fwhm, x, y, xlims, ylims):

    N = 10

    img = np.zeros((N, N), dtype='float32') + np.nanmean(delta_fwhm)

    return img

    #x0, x1 = xlims
    #y0, y1 = ylims

    # flin = interpolate.interp2d(x, y, delta_fwhm, kind='linear',
    #                            bounds_error=False,
    #                            fill_value=np.mean(delta_fwhm))

    #xQ = np.linspace(x0, x1, N)
    #yQ = np.linspace(y0, y1, N)

    #xxQ, yyQ = np.meshgrid(xQ, yQ)

    #img = flin(xQ, yQ)

    # return img
