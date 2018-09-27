#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

json files handling utilities.

Created on Tue Mar 27 14:25:43 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import json
#import ast
import yaml
#from decimal import Decimal
# END IMPORT


def load_jsonfile(jsonfile, useyaml=False):
    """ """

    with open(jsonfile) as json_data:

        if useyaml:
            json_str = json_data.read()
            pydict = yaml.safe_load(json_str)
        else:
            pydict = json.load(json_data, encoding='ascii',
                               parse_float=float, parse_int=int)
    return pydict


def dumps_to_json(pydict):
    """ """
    return json.dumps(pydict, indent=4)


def save_jsonfile(pydict, jsonfile):
    """ """
    json = dumps_to_json(pydict)
    with open(jsonfile, 'wa') as f:
        print >> f, json
