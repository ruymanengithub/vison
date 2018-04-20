#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:55:48 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""
#IMPORT STUFF
import os
#END IMPORT

def get_bitsize_file(item):
    """ """
    return os.stat(item).st_size

def get_bitsize(filelist):
    """ """
    bitsize = 0
    for item in filelist:
        bitsize += get_bitsize_file(item)
    return bitsize
