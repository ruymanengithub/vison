#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 14:25:43 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import json
# END IMPORT


def load_metafile(metafile):
    """ """
    with open(metafile) as json_data:
        meta = json.load(json_data,encoding='ascii')['metadata']
    return meta
    