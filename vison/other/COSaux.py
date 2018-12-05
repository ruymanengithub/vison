#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to COSMETICS00

Created on Wed Dec 06 18:04:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import copy

from vison.datamodel import cdp
from vison.flat import BF01aux
from vison.plot import figclasses
from vison.plot import trends
#from MOT_FF import extract_overscan_profiles
# END IMPORT





def get_COS_figs():
    COS_figs = OrderedDict()
    
    COS_figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return COS_figs


def get_CDP_lib():
       
    CDP_lib = OrderedDict()
    
    
    return CDP_lib
