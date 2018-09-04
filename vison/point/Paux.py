#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to PointTask.

Created on Tue Nov 14 13:54:34 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT




lock_cdp = cdp.Tables_CDP()
lock_cdp.rootname = 'LOCK_TB_%s'

CDP_lib = dict(LOCK_TB=lock_cdp)