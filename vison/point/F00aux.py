#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to FOCUS00.

Created on Tue Nov 14 13:54:34 2017

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

def gt_check_offsets_dict(wave):
    return dict(stats=['offset_pre','offset_ove'],
                          trendaxis='time',
                          figname='FOCUS00_%i_offset_vs_time.png' % wave,
                          caption='FOCUS00\_%i: offset vs. time.' % wave,
                          meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='FOCUS01\_%i-checks: offsets' % wave))

def gt_check_std_dict(wave):
    return dict(stats=['std_pre','std_ove'],
                      trendaxis='time',
                          figname='FOCUS00_%i_std_vs_time.png' % wave,
                          caption='FOCUS00\_%i: std vs. time.' % wave,
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='FOCUS00\_%i-checks: std' % wave))
 
def gt_F00figs(wave):
    F00figs = dict()
    F00figs['F00checks_offsets'] = [trends.pl_basic_checkstat,gt_check_offsets_dict(wave)]
    F00figs['F00checks_stds'] = [trends.pl_basic_checkstat,gt_check_std_dict(wave)]
    F00figs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]
    return F00figs

