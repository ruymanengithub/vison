#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to TP01.

Created on Wed Jan 31 15:11:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.plot import figclasses
from vison.plot import trends

# END IMPORT


check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='TP01_offset_vs_time.png',
                          caption='TP01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='TP01-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='TP01_std_vs_time.png',
                      caption='TP01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='TP01-checks: std',
                                ylim=trends.RON_lims))

check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='TP01_injlevel_vs_time.png',
                           caption='TP01: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='TP01-checks: Injection Level'))

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='TP01_injstd_vs_time.png',
                         caption='TP01: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='TP01-checks: Injection STD')
                         )


def get_TP01figs():
    TP01figs = dict()
    TP01figs['TP01checks_offsets'] = [
        trends.Fig_Basic_Checkstat, check_offsets_dict]
    TP01figs['TP01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    TP01figs['TP01checks_injlevel'] = [
        trends.Fig_Basic_Checkstat, check_injlevel_dict]
    TP01figs['TP01checks_injstd'] = [trends.Fig_Basic_Checkstat, check_injstd_dict]
    TP01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return TP01figs
