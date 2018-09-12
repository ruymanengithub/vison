#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to CHINJ02.

Created on Tue Jan 30 14:28:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

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
                          figname='CHINJ02_offset_vs_time.png',
                          caption='CHINJ02: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='CHINJ02-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='CHINJ02_std_vs_time.png',
                      caption='CHINJ02: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='CHINJ02-checks: std',
                                ylim=trends.RON_lims)
                      )
check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='CHINJ02_injlevel_vs_time.png',
                           caption='CHINJ02: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='CHINJ02-checks: Injection Level')
                           )

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='CHINJ02_injstd_vs_time.png',
                         caption='CHINJ02: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='CHINJ02-checks: Injection STD')
                         )


CH02figs = dict()
CH02figs['CH02checks_offsets'] = [
    trends.Fig_Basic_Checkstat, check_offsets_dict]
CH02figs['CH02checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
CH02figs['CH02checks_injlevel'] = [
    trends.Fig_Basic_Checkstat, check_injlevel_dict]
CH02figs['CH02checks_injstd'] = [trends.Fig_Basic_Checkstat, check_injstd_dict]
CH02figs['BlueScreen'] = [figclasses.BlueScreen, dict()]

CDP_lib = dict()