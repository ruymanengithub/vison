#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to NL01.

Created on Tue Jan 30 18:43:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st

from vison.plot import baseclasses as plbaseclasses
from vison.plot import trends

# END IMPORT

check_offsets_dict= dict(stats=['offset_pre','offset_ove'],
                          trendaxis='time',
                          figname='NL01_offset_vs_time.png',
                          caption='NL01: offset vs. time.',
                          meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='NL01-checks: offsets'))

check_std_dict = dict(stats=['std_pre','std_ove'],
                      trendaxis='time',
                          figname='NL01_std_vs_time.png',
                          caption='NL01: std vs. time.',
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='NL01-checks: std' ))


check_img_flu_dict = dict(stats=['flu_med_img'],
                        trendaxis='exptime',
                          figname='NL01_flu_vs_exptime.png',
                          caption='NL01: Fluence vs. exposure time.',
                          meta=dict(doLegend=True,
                                doNiceXDate=False,
                                suptitle='NL01-checks: Fluence',
                                ylabel='[ADU]',
                                xlabel='exptime [s]'))

check_img_var_dict = dict(stats=['flu_var_img'],
                        trendaxis='exptime',
                          figname='NL01_var_vs_exptime.png',
                          caption='NL01: Variance vs. exposure time.',
                          meta=dict(doLegend=False,
                                doNiceXDate=False,
                                suptitle='NL01-checks: Variance',
                                ylabel='[ADU^2]',
                                xlabel='exptime [s]'))

NL01figs = dict()
NL01figs['NL01checks_offsets'] = [trends.pl_basic_checkstat,check_offsets_dict.copy()]
NL01figs['NL01checks_stds'] = [trends.pl_basic_checkstat,check_std_dict.copy()]
NL01figs['NL01checks_flu'] = [trends.pl_basic_checkstat,check_img_flu_dict.copy()]
NL01figs['NL01checks_var'] = [trends.pl_basic_checkstat,check_img_var_dict.copy()]
NL01figs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]

