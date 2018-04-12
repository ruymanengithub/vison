#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to TP02.

Created on Wed Jan 31 14:59:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.plot import baseclasses as plbaseclasses
from vison.plot import trends

# END IMPORT


check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='TP02_offset_vs_time.png',
                          caption='TP02: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='TP02-checks: offsets'))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='TP02_std_vs_time.png',
                      caption='TP02: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='TP02-checks: std'))

check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='TP02_injlevel_vs_time.png',
                           caption='TP02: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='TP02-checks: Injection Level'))

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='TP02_injstd_vs_time.png',
                         caption='TP02: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='TP02-checks: Injection STD')
                         )


TP02figs = dict()
TP02figs['TP02checks_offsets'] = [
    trends.pl_basic_checkstat, check_offsets_dict]
TP02figs['TP02checks_stds'] = [trends.pl_basic_checkstat, check_std_dict]
TP02figs['TP02checks_injlevel'] = [
    trends.pl_basic_checkstat, check_injlevel_dict]
TP02figs['TP02checks_injstd'] = [trends.pl_basic_checkstat, check_injstd_dict]
TP02figs['BlueScreen'] = [plbaseclasses.BlueScreen, dict()]
