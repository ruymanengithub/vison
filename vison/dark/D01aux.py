#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to DARK01.

Created on Sat Jan 27 16:11:00 2018

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


check_offsets_dict = dict(stats=['offset_pre','offset_ove'],
                          figname='DARK01_offset_vs_time.png',
                          caption='DARK01: offset vs. time.',
                          meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='DARK01-checks: offsets'))

check_std_dict = dict(stats=['std_pre','std_ove'],
                          figname='DARK01_std_vs_time.png',
                          caption='DARK01: std vs. time.',
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='DARK01-checks: std')
                          )
 
check_flu_dict = dict(stats=['chk_flu_img'],
                          figname='DARK01_FLU_vs_time.png',
                          caption='DARK01: Fluence vs. time.',
                          meta=dict(doLegend=False,
                                doNiceXDate=True,
                                suptitle='DARK01-checks: Image Area Fluence.')
                          )
    

    
D01figs = dict()
D01figs['D01checks_offsets'] = [trends.pl_basic_checkstat,check_offsets_dict]
D01figs['D01checks_stds'] = [trends.pl_basic_checkstat,check_std_dict]
D01figs['D01checks_flu'] = [trends.pl_basic_checkstat,check_flu_dict]
D01figs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]