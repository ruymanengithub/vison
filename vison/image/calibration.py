#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:54:28 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""


def load_FITS_CDPs(FDict, dataclass, **kwargs):
    """Dummy function to load CDPs for all 3 CCDs.
    Input is of type dict(CCD1='',CCD2='',CCD3='')"""
    # Load File for each CCD

    data = dict()
    for CCDindex in range(1, 4):
        CCDkey = 'CCD%i' % CCDindex
        if CCDkey in FDict.keys():
            data[CCDkey] = dataclass(FDict[CCDkey], **kwargs)
    return data
