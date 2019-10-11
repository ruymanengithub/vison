#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to NL01.

Created on Tue Jan 30 18:43:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT

def get_CDP_lib():
    nl_tb_cdp = cdp.Tables_CDP()
    nl_tb_cdp.rootname = 'NL01_RESULTS_TABLE'
    
    curves_cdp = cdp.CDP()
    curves_cdp.rootname = 'NL01_CURVES'
    
    CDP_lib = dict(NL_TB=nl_tb_cdp,
                   NL_CURVES=curves_cdp)
    return CDP_lib


check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='NL01_offset_vs_time.png',
                          caption='NL01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='NL01-checks: offsets',
                                    ylim=trends.offset_lims))

check_deltaoff_dict = dict(stats=['deltaoff_pre', 'deltaoff_ove'],
                          trendaxis='time',
                          figname='NL01_deltaoff_vs_time.png',
                          caption='NL01: $\delta$offset vs. time. Offset value in each frame minus the average value.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='NL01-checks: delta-offsets',
                                    ylim=[-10.,10.]))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='NL01_std_vs_time.png',
                      caption='NL01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='NL01-checks: std',
                                ylim=trends.RON_lims))


check_img_flu_dict = dict(stats=['flu_med_img'],
                          trendaxis='exptime',
                          figname='NL01_flu_vs_exptime.png',
                          caption='NL01: Fluence vs. exposure time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=False,
                                    suptitle='NL01-checks: Fluence',
                                    ylabel='[ADU]',
                                    xlabel='exptime [s]'))

check_img_std_dict = dict(stats=['flu_std_img'],
                          trendaxis='exptime',
                          figname='NL01_imgstd_vs_exptime.png',
                          caption='NL01: Image area STD vs. exposure time.',
                          meta=dict(doLegend=False,
                                    doNiceXDate=False,
                                    suptitle='NL01-checks: Image STD',
                                    ylabel='[ADU]',
                                    xlabel='exptime [s]'))

NL_curves_dict = dict(
    figname='NL01_curves.png',
    caption='NL01: Relative non-linearity vs. fluence and best fit curves.',
    meta=dict(doLegend=True,
              ylabel='Z [percentage]',
              xlabel=r'$Y_{NL} [kADU]$',
              suptitle='NL01: Non-Linearity Curves.',
              ylim=[-10.,10.],
                   corekwargs=dict(data=dict(marker='.',linestyle='',color='b'),
                                   fit=dict(marker='',linestyle='--',color='r')))
    )


def get_NL01figs():
    NL01figs = dict()
    NL01figs['NL01checks_offsets'] = [
        trends.Fig_Basic_Checkstat, check_offsets_dict.copy()]
    NL01figs['NL01checks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, check_deltaoff_dict.copy()]
    NL01figs['NL01checks_stds'] = [
        trends.Fig_Basic_Checkstat, check_std_dict.copy()]
    NL01figs['NL01checks_flu'] = [
        trends.Fig_Basic_Checkstat, check_img_flu_dict.copy()]
    NL01figs['NL01checks_imgstd'] = [
        trends.Fig_Basic_Checkstat, check_img_std_dict.copy()]
    NL01figs['NL01_fit_curves'] = [
            figclasses.Fig_Beam2DPlot, NL_curves_dict.copy()]
    NL01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return NL01figs