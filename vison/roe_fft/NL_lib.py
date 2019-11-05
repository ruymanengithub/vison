#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Some utilities for analysis of Non-Linearity of ROE and ROE-TAB.

Created on Wed Mar 21 16:32:19 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from scipy import optimize as opt
import copy

from pylab import plot, show
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

# END IMPORT


def myplot_NL(datasets, xlabel='', ylabel='', ylim=[], title='', figname=''):
    """ """

    fig = plt.figure()
    ax = fig.add_subplot(111)

    for dkey in list(datasets.keys()):
        dataset = datasets[dkey]
        x = dataset['x']
        y = dataset['y']
        try:
            ls = dataset['ls']
        except BaseException:
            ls = ''
        try:
            color = dataset['color']
        except BaseException:
            color = None
        try:
            marker = dataset['marker']
        except BaseException:
            marker = None

        ax.plot(x, y, color=color, ls=ls, marker=marker, label=dkey)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if len(ylim) > 0:
        ax.set_ylim(ylim)

    handles, labels = ax.get_legend_handles_labels()

    ax.legend(handles, labels, loc='best')

    ax.set_title(title)

    plt.tight_layout()

    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()


def find_NL_pol(x, y, deg=2, sigma=None, Full=False, debug=False):
    """ """

    pfitLin = np.polyfit(x, y, deg=1, w=sigma)
    yLin = np.polyval(pfitLin, x)
    pfitLin_throughzero = copy.deepcopy(pfitLin)
    pfitLin_throughzero[-1] = 0.
    yLin_throughzero = np.polyval(pfitLin_throughzero, x)

    NLin = (y - yLin) / yLin_throughzero

    pfitNLin = np.polyfit(x, NLin, deg=deg, w=sigma)

    if debug:
        plotset = dict(Nlin=dict(x=x, y=y - yLin,
                                 ls='-', color='b'))
        myplot_NL(plotset, xlabel='mV', ylabel='y-yLin', title='', figname='')

    if Full:
        return pfitNLin, (x, y, yLin, yLin_throughzero, NLin)
    else:
        return pfitNLin
