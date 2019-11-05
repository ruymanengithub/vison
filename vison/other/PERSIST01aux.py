#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxilisary Functions and resources to PERSIST01.


Created on Thu Feb  7 18:28:47 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd
from matplotlib import cm
import copy

from vison.plot import figclasses, baseplotclasses
from vison.plot import trends
from vison.datamodel import cdp
from vison.support.files import cPickleRead


# END IMPORT

check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='PERSIST01_offset_vs_time.png',
                          caption='PERSIST01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='PERSIST01-checks: offsets',
                                    ylim=trends.offset_lims))
check_deltaoff_dict = dict(
    stats=[
        'deltaoff_pre',
        'deltaoff_ove'],
    trendaxis='time',
    figname='PERSIST01_deltaoff_vs_time.png',
    caption='PERSIST01: $\delta$offset vs. time. Offset value in each frame minus the average value.',
    meta=dict(
            doLegend=True,
            doNiceXDate=True,
            suptitle='PERSIST01-checks: delta-offsets',
            ylim=[
                -10.,
                10.]))


check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='PERSIST01_std_vs_time.png',
                      caption='PERSIST01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='PERSIST01-checks: std',
                                ylim=trends.RON_lims)
                      )

check_avgflu_dict = dict(stats=['chk_avgflu_img'],
                         trendaxis='time',
                         figname='PERSIST01_avgbgd_vs_time.png',
                         caption='PERSIST01: Background level vs. time.',
                         meta=dict(doLegend=True,
                                   doNiceXDate=True,
                                   suptitle='PERSIST01-checks: background',
                                   )
                         )

P01satmask_dict = dict(
    figname='PERSIST01_SATMASK_2Dimgshow.png',
    caption='PERSIST01: Saturation Masks for the CCDs.',
    meta=dict(doLegend=False,
              doColorbar=True,
              suptitle='PERSIST01: Saturation Masks',
              corekwargs=dict(cmap=cm.gray, aspect='auto', norm=None,
                              origin='lower left'))
)


P01whiskers_dict = dict(
    figname='PERSIST01_whiskersplot.png',
    caption='Persistence "whisker-box" plot. Statistics of pixel values in saturated regions before and after the saturation. ' +
    'x-axis labeled in seconds before/after the hard-latent. All exposures shown are 565 seconds long. The box extends ' +
    'from the lower (25\%) to upper quartile (75\%) values of the data, with an orange line at the median. The whiskers ' +
    'extend from the mean value to +/- 1 standard deviation.',
    meta=dict(
        ylim=[-10., 10.],
        doRotateXLabels=True),
    data=None
)


class P01WhiskersPlot(baseplotclasses.BeamPlotYvX):
    """ """

    def _ax_core_funct(self, ax, CQdict):
        """ """

        stats_trans = dict(label='deltasec',
                           mean='mean',
                           med='p50',
                           q1='p25',
                           q3='p75')

        Nref = len(CQdict['REF']['ObsID'])
        Nlat = len(CQdict['LAT']['ObsID'])

        def feed_stats(CQdict, tag, i):
            istats = dict()
            for key in list(stats_trans.keys()):
                istats[key] = CQdict[tag][stats_trans[key]][i]
            istats['label'] = '%i' % CQdict[tag]['deltasec'][i]
            istats['whislo'] = CQdict[tag]['mean'][i] -\
                CQdict[tag]['std'][i]
            istats['whishi'] = CQdict[tag]['mean'][i] +\
                CQdict[tag]['std'][i]
            istats['fliers'] = []

            return istats

        statsREF = []

        for i in range(Nref):
            # REFERENCE
            statsREF.append(feed_stats(CQdict, 'REF', i))
        # ax.bxp(statsREF)

        statsLAT = []
        for j in range(Nlat):
            statsLAT.append(feed_stats(CQdict, 'LAT', j))

        stats = statsREF + statsLAT

        ax.bxp(stats)

        ax.axvline(x=Nref + 0.5, linestyle='--', color='r')

        handle, label = None, None
        return handle, label


class Fig_BeamP01Whiskers(figclasses.Fig):
    plotclass = P01WhiskersPlot


def get_P01figs():
    P01figs = dict()
    P01figs['P01checks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
    P01figs['P01checks_deltaoff'] = [trends.Fig_Basic_Checkstat, check_deltaoff_dict]
    #P01figs['B01checks_stds'] = [plB01check,dict(stat='std')]
    P01figs['P01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    P01figs['P01checks_avgflu'] = [trends.Fig_Basic_Checkstat, check_avgflu_dict]
    P01figs['P01satmasks'] = [figclasses.Fig_BeamImgShow, P01satmask_dict]
    P01figs['P01whiskers'] = [Fig_BeamP01Whiskers, P01whiskers_dict]
    P01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return P01figs


def get_CDP_lib():

    satarea_cdp = cdp.Tables_CDP()
    satarea_cdp.rootname = 'PERSIST01_SATAREAS'

    CDP_lib = OrderedDict()
    CDP_lib['SATAREA_TB'] = satarea_cdp

    CDP_lib['PERSIST_STATS'] = cdp.CDP()
    CDP_lib['PERSIST_STATS'].rootname = 'PERSIST01_STATS'

    return CDP_lib


def _get_stats(pick_list, satmaskobj):
    """ """

    Quads = ['E', 'F', 'G', 'H']
    stats = OrderedDict()
    Nobs = len(pick_list)

    statkeys = ['mean', 'p25', 'p50', 'p75', 'std', 'N']

    qmasks = OrderedDict()
    for Q in Quads:
        stats[Q] = OrderedDict()

        qmasks[Q] = satmaskobj.get_quad(Q, canonical=False, extension=-1).copy()

        for statkey in statkeys:
            stats[Q][statkey] = np.zeros(Nobs, dtype='float32') + np.nan

    for i in range(Nobs):

        ccdpick = pick_list[i]
        ccdobj = cPickleRead(ccdpick)

        bgdobj = copy.deepcopy(ccdobj)

        bgdobj.get_mask(satmaskobj.extensions[-1].data)

        for Q in Quads:

            qimg = ccdobj.get_quad(Q, canonical=False, extension=-1)

            bgd_stats = bgdobj.get_stats(Q, sector='img', statkeys=['median', 'mean'],
                                         ignore_pover=True, extension=-1)
            bgd = bgd_stats[1]

            qimg -= bgd

            data = qimg.data[np.where(qmasks[Q] == 1)].copy()

            stats[Q]['mean'][i] = np.mean(data)
            stats[Q]['p25'][i] = np.percentile(data, 25)
            stats[Q]['p50'][i] = np.percentile(data, 50)
            stats[Q]['p75'][i] = np.percentile(data, 75)
            stats[Q]['std'][i] = np.std(data)
            stats[Q]['N'][i] = np.size(data)

    return stats
