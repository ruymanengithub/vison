#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 15:38:37 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop

from vison.datamodel import core
from vison import __version__
# END IMPORT

def DataDict_builder(explog, inputs):
    """ """

    # Build DataDict

    dd = core.FpaDataDict()
    # Load Metadata
    dd.meta = dict(inputs=inputs, vison=__version__)
    # Load Exposure Log
    dd.loadExpLog(explog)

    return dd