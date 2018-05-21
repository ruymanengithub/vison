#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Performance parameters of the ROE+CCDs.
Compilation of CCD offsets, offset gradients, RONs... used for 
checks.

Created on Wed Nov  1 09:57:44 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from collections import OrderedDict
# ENDIMPORT

CCDs = ['CCD1','CCD2','CCD3']
Quads = ['E','F','G','H']


offsets_margins = [-200, 200]  # ADU


def get_offsets(BLOCKID=None):
     return dict(CCD1=2000, CCD2=2000, CCD3=2000)  # ADU

def get_offsets_lims(offsets,offsets_margins):
    """ """
    offsets_lims = OrderedDict()
    for CCD in CCDs:
        offsets_lims[CCD] = OrderedDict()
        for Q in Quads:
            offsets_lims[CCD][Q] = [offsets[CCD][Q] +
                            offsets_margins[0], offsets[CCD][Q]+offsets_margins[1]]
    return offsets_lims


# offsets_gradients = dict(CCDX=[[prescan-prescan+-margin],[image-prescan+-margin],[overscan-prescan+-margin]])
offsets_gradients = dict(CCD1=dict(pre=[0, 0], img=[5, 10], ove=[5, 10]), 
                         CCD2=dict(pre=[0, 0], img=[5, 10], ove=[5, 10]),
                         CCD3=dict(pre=[0, 0], img=[5, 10], ove=[5, 10]))  # ADU
RONs = dict(CCD1=1.4, CCD2=1.4, CCD3=1.4)  # ADU, STD
RONs_margins = [-0.2, 0.2]  # ADU, STD
RONs_lims = dict()
for CCD in CCDs:
    RONs_lims[CCD] = [RONs[CCD] +
                         RONs_margins[0], RONs[CCD]+RONs_margins[1]]


gains = dict(CCD1=3.3, CCD2=3.3, CCD3=3.3)



def_perf_rdout = dict(
#                offsets=offsets, 
                 offsets_margins=offsets_margins,
#                  offsets_lims=offsets_lims,
                 offsets_gradients=offsets_gradients,
                 RONs=RONs, 
                 RONs_margins=RONs_margins,
                 RONs_lims=RONs_lims)

def get_perf_rdout(BLOCKID):
    """ """
    
    perf_rdout = OrderedDict()
    perf_rdout.update(def_perf_rdout)
    
    offsets = get_offsets(BLOCKID)
    
    perf_rdout['offsets'] = offsets
    perf_rdout['offsets_lims'] = get_offsets_lims(offsets,offsets_margins)
    
    return perf_rdout
