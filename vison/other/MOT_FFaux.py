#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to MOT_FF

Created on Tue Jul 31 17:50:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import copy

from vison.flat import BF01aux
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT



def gt_MOT_FF_figs(test):
    BF01_figs = BF01aux.gt_BF01figs(test)
    
    MOT_FF_figs = OrderedDict()
    
    for key in BF01_figs.keys():
        if 'BF01' in key:
            nkey = st.replace(key,'BF01','MOT_FF')
        else:
            nkey = key
        MOT_FF_figs[nkey] = copy.deepcopy(BF01_figs[key])
    
    return MOT_FF_figs
