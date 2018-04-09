#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:25:43 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import json
#from decimal import Decimal
# END IMPORT


def load_jsonfile(jsonfile):
    """ """
    with open(jsonfile) as json_data:
        pydict = json.load(json_data,encoding='ascii',parse_float=float,parse_int=int)
    return pydict

