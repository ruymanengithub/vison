#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to MOT_FF

Created on Tue Jul 31 17:50:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import copy

from vison.flat import BF01aux
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT


def extract_overscan_profiles(ccdobj, thresholds, direction='serial'):
    """ """

    ixjump = 10
    Qs = ['E', 'F', 'G', 'H']

    if direction == 'serial':
        detedge = ccdobj.NAXIS1 / 2 - ccdobj.overscan
    elif direction == 'parallel':
        detedge = ccdobj.NAXIS2 / 2 - ccdobj.voverscan

    x = np.arange(25) + detedge - ixjump + 1

    profiles = dict()
    profiles['ixjump'] = ixjump

    for Q in Qs:

        imgdata = ccdobj.get_quad(Q, canonical=True)

        if direction == 'serial':
            strip = imgdata[-ccdobj.overscan - ixjump:-ccdobj.overscan + 15, :].copy()
        elif direction == 'parallel':
            strip = imgdata[ccdobj.prescan:-ccdobj.overscan,
                            ccdobj.NrowsCCD - ixjump:ccdobj.NrowsCCD + 15].transpose().copy()

        injection = np.mean(strip[0:ixjump, :], axis=0)
        bias = np.mean(strip[ixjump + 3:, :], axis=0)

        ixgood = np.where((injection <= thresholds[1]) & (injection >= thresholds[0]))
        #print '%i rows averaged' % (len(ixgood[0]),)

        strip = strip[:, ixgood[0]].copy()

        nstrip = (strip - bias[ixgood]) / (injection[ixgood] - bias[ixgood])

        profile = np.mean(nstrip, axis=1)

        if isinstance(profile, np.ma.masked_array):
            ixgood2 = np.where(profile.mask == False)
            profiles[Q] = dict(y=profile[ixgood2].data.copy(),
                               x=x[ixgood2].copy())
        else:
            profiles[Q] = dict(y=profile.copy(),
                               x=x.copy())

    return profiles


prof_HER_ser_dict = dict(
    figname='MOT_FF_profs_HER_ser_allOBSIDs.png',
    caption='MOT\_FF: HER profiles, serial direction.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim=[-0.002, 0.005],
              suptitle='MOT\_FF: Serial HER.')
)

prof_HER_ver_dict = dict(
    figname='MOT_FF_profs_HER_ver_allOBSIDs.png',
    caption='MOT\_FF: HER profiles, vertical/parallel direction.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='Col. [pix]',
              ylim=[-0.002, 0.005],
              suptitle='MOT\_FF: Vertical/Parallel HER.')
)


def gt_MOT_FF_figs(test):

    BF01_figs = BF01aux.gt_BF01figs(test)

    MOT_FF_figs = OrderedDict()

    for key in list(BF01_figs.keys()):
        #        if 'BF01' in key:
        #            nkey = key.replace('BF01','MOT_FF')
        #        else:
        #            nkey = key
        MOT_FF_figs[key] = copy.deepcopy(BF01_figs[key])

    MOT_FF_figs['MOTFF_HER_ser'] = [
        figclasses.Fig_Beam2DPlot, prof_HER_ser_dict]

    MOT_FF_figs['MOTFF_HER_ver'] = [
        figclasses.Fig_Beam2DPlot, prof_HER_ver_dict]

    return MOT_FF_figs
