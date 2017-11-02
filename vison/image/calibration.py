#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:54:28 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

def load_CDPs(FDict,dataclass):

    # Load File for each CCD
    
    data = dict()
    for CCDindex in range(1,4):
        CCDkey = 'CCD%i' % CCDindex
        if CCDkey in FDict.keys():
            data[CCDkey] = dataclass(fitsfile=FDict[CCDkey])
    return data