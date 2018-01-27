#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BIAS01.

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


check_offsets_dict = dict(stats=['offset_pre','offset_img','offset_ove'],
                          figname='BIAS01_offset_vs_time.png',
                          caption='BIAS01: offset vs. time.',
                          meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='BIAS01-checks: offsetes'))

check_std_dict = dict(stats=['std_pre','std_img','std_ove'],
                          figname='BIAS01_std_vs_time.png',
                          caption='BIAS01: std vs. time.',
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='BIAS01-checks: std')
                          )
 
B01figs = dict()
B01figs['B01checks_offsets'] = [trends.pl_basic_checkstat,check_offsets_dict]
#B01figs['B01checks_stds'] = [plB01check,dict(stat='std')]
B01figs['B01checks_stds'] = [trends.pl_basic_checkstat,check_std_dict]
B01figs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]