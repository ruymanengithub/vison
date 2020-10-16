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
    """Extract profiles across over scan transition from image area.
    

    Parameters
    ----------
    
    :param thresholds: thresholds in image fluence to extract the profiles. The profiles are
                        relative to the image area fluence, and that's why these thresholds
                        are important.
    :type thresholds: list of 2 floats, min/max thresholds.
    :param ccdobj: CCD object. The Data.
    :type ccdobj: instance of vison.datamodel.ccd.CCD.

    :param direction: one of {'serial', 'parallel'}
    :type direction: str

    """

    ixjump = 10 
    Npixprof = 25

    Qs = ccdobj.Quads # Quadrants

    if direction == 'serial':
        detedge = ccdobj.NAXIS1 // 2 - ccdobj.overscan
    elif direction == 'parallel':
        detedge = ccdobj.NAXIS2 // 2 - ccdobj.voverscan

    x = np.arange(Npixprof) + detedge - ixjump + 1

    profiles = dict()
    profiles['ixjump'] = ixjump

    for Q in Qs:

        imgdata = ccdobj.get_quad(Q, canonical=True)

        if direction == 'serial':
            strip = imgdata[-ccdobj.overscan - ixjump:-ccdobj.overscan + Npixprof-ixjump, :].copy()
        elif direction == 'parallel':
            strip = imgdata[ccdobj.prescan:-ccdobj.overscan,
                ccdobj.NrowsCCD - ixjump:ccdobj.NrowsCCD - ixjump + Npixprof].transpose().copy()

        injection = np.mean(strip[0:ixjump, :], axis=0)
        bias = np.mean(strip[ixjump + 5:, :], axis=0)

        ixgood = np.where((injection <= thresholds[1]) & (injection >= thresholds[0])) # fluence selection
        #print '%i rows averaged' % (len(ixgood[0]),)


        nstrip = (strip[:, ixgood[0]] - bias[ixgood]) / (injection[ixgood] - bias[ixgood]) # normalization

        profile = np.mean(nstrip, axis=1) # stacking

        stop()

        if isinstance(profile, np.ma.masked_array):
            ixgood2 = np.where(profile.mask == False)
            profiles[Q] = dict(y=profile[ixgood2].data.copy(),
                               x=x[ixgood2].copy())
        else:
            profiles[Q] = dict(y=profile.copy(),
                               x=x.copy())

    return profiles


def extract_prescan_profiles(ccdobj, thresholds):
    """Extract profiles across prescan transition to image area (serial direction).
    

    Parameters
    ----------
    
    :param thresholds: thresholds in image fluence to extract the profiles. The profiles are
                        relative to the image area fluence, and that's why these thresholds
                        are important.
    :type thresholds: list of 2 floats, min/max thresholds.
    :param ccdobj: CCD object. The Data.
    :type ccdobj: instance of vison.datamodel.ccd.CCD.

    """

    ixjump = 10 
    Npixprof = 25

    Qs = ccdobj.Quads # Quadrants

    detedge = ccdobj.prescan

    x = np.arange(Npixprof) + detedge - ixjump + 1

    profiles = dict()
    profiles['ixjump'] = ixjump

    for Q in Qs:

        imgdata = ccdobj.get_quad(Q, canonical=True)
        
        strip = imgdata[ccdobj.prescan - ixjump:ccdobj.prescan + Npixprof-ixjump, :].copy()

        injection = np.mean(strip[ixjump+3:, :], axis=0)
        bias = np.mean(strip[0:ixjump-1, :], axis=0)
        

        ixgood = np.where((injection >= thresholds[0]) & (injection <= thresholds[1])) # fluence selection
        #print '%i rows averaged' % (len(ixgood[0]),)

        nstrip = (strip[:, ixgood[0]] - bias[ixgood]) / (injection[ixgood] - bias[ixgood]) # normalization

        profile = np.mean(nstrip, axis=1)-1. # stacking

        if isinstance(profile, np.ma.masked_array):
            ixgood2 = np.where(profile.mask == False)
            profiles[Q] = dict(y=profile[ixgood2].data.copy(),
                               x=x[ixgood2].copy())
        else:
            profiles[Q] = dict(y=profile.copy(),
                               x=x.copy())

    return profiles

def extract_transcan_profiles(ccdobj, thresholds, direction='serial', scan='over'):
    """Extract profiles across pre / over scan transitions to / from image area.
    

    Parameters
    ----------
    
    :param thresholds: thresholds in image fluence to extract the profiles. The profiles are
                        relative to the image area fluence, and that's why these thresholds
                        are important.
    :type thresholds: list of 2 floats, min/max thresholds.
    :param ccdobj: CCD object. The Data.
    :type ccdobj: instance of vison.datamodel.ccd.CCD.

    :param direction: one of {'serial', 'parallel'}. 
        If scan=='pre', direction can only be 'pre'.
    :type direction: str
    :param scan: one of {'over', 'pre'}
    :type scan: str

    """

    if scan == 'over':
        return extract_overscan_profiles(ccdobj, thresholds, direction)
    elif scan == 'pre':
        assert direction == 'serial'
        return extract_prescan_profiles(ccdobj, thresholds)

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
