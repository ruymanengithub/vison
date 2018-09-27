#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to DARK01.

Created on Sat Jan 27 16:11:00 2018

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
                          figname='DARK01_offset_vs_time.png',
                          caption='DARK01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='DARK01-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='DARK01_std_vs_time.png',
                      caption='DARK01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='DARK01-checks: std',
                                ylim=trends.RON_lims)
                      )

check_flu_dict = dict(stats=['chk_flu_img'],
                      trendaxis='time',
                      figname='DARK01_FLU_vs_time.png',
                      caption='DARK01: Fluence vs. time.',
                      meta=dict(doLegend=False,
                                doNiceXDate=True,
                                suptitle='DARK01-checks: Image Area Fluence.')
                      )


D01figs = dict()
D01figs['D01checks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
D01figs['D01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
D01figs['D01checks_flu'] = [trends.Fig_Basic_Checkstat, check_flu_dict]
D01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
