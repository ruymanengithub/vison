#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Performance parameters of the ROE+CCDs.
Compilation of CCD offsets, offset gradients, RONs... used for 
checks.

Created on Wed Nov  1 09:57:44 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from collections import OrderedDict
import os
from pdb import set_trace as stop

from vison.support import vjson
from vison import blocks
# ENDIMPORT

CCDs = ['CCD1', 'CCD2', 'CCD3']
Quads = ['E', 'F', 'G', 'H']

offsets_UNK = OrderedDict(CCD1=OrderedDict(E=2800.,
                                           F=2800.,
                                           G=2800.,
                                           H=2800.))
offsets_UNK['CCD2'] = offsets_UNK['CCD1'].copy()
offsets_UNK['CCD3'] = offsets_UNK['CCD1'].copy()

offsets_margins = [-200, 200]  # ADU

blockspath = blocks.__path__[0]


def get_offsets(BLOCKID=None):

    if BLOCKID is None:
        return offsets_UNK

    offsets_json = '%s_offsets.json' % BLOCKID
    fpath = os.path.join(blockspath, offsets_json)
    try:
        offsets = vjson.load_jsonfile(fpath)
    except:
        stop()
    return offsets  # ADU


def get_offsets_lims(offsets, offsets_margins):
    """ """
    offsets_lims = OrderedDict()
    for CCD in CCDs:
        offsets_lims[CCD] = OrderedDict()
        for Q in Quads:
            offsets_lims[CCD][Q] = [offsets[CCD][Q] +
                                    offsets_margins[0], offsets[CCD][Q]+offsets_margins[1]]
    return offsets_lims


# offsets_gradients = dict(CCDX=[[prescan-prescan+-margin],[image-prescan+-margin],[overscan-prescan+-margin]])
offsets_gradients = dict(CCD1=dict(pre=[0, 0], img=[-10, 10], ove=[-10, 10]),
                         CCD2=dict(pre=[0, 0], img=[-10, 10], ove=[-10, 10]),
                         CCD3=dict(pre=[0, 0], img=[-10, 10], ove=[-10, 10]))  # ADU

#RONs = dict(CCD1=1.4, CCD2=1.4, CCD3=1.4)  # ADU, STD
#RONs_margins = [-0.2, 0.2]  # ADU, STD

RONs_lims = dict()
for CCD in CCDs:
    #RONs_lims[CCD] = [RONs[CCD] +
    #                  RONs_margins[0], RONs[CCD]+RONs_margins[1]]
    RONs_lims[CCD] = [0.86,1.44]

gains = dict(CCD1=3.3, CCD2=3.3, CCD3=3.3)


def_perf_rdout = dict(
    #                offsets=offsets,
    offsets_margins=offsets_margins,
    #                  offsets_lims=offsets_lims,
    offsets_gradients=offsets_gradients,
    #RONs=RONs,
    #RONs_margins=RONs_margins,
    RONs_lims=RONs_lims)


def get_perf_rdout(BLOCKID):
    """ """

    perf_rdout = OrderedDict()
    perf_rdout.update(def_perf_rdout)

    offsets = get_offsets(BLOCKID)
    perf_rdout['offsets'] = offsets
    perf_rdout['offsets_lims'] = get_offsets_lims(offsets, offsets_margins)

    return perf_rdout
