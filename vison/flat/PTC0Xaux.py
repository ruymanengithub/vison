#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to PTC0X.

Created on Tue Jan 30 17:51:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import pandas as pd
import copy
from matplotlib import cm
from scipy import stats

from vison.plot import figclasses
from vison.plot import trends
from vison.datamodel import cdp
from vison.datamodel.ccd import CCD
from vison.datamodel import ccd_aux
# END IMPORT


class CCDclone(CCD):

  def __init__(self, ccd2clone):
    """ """

    attrs2clone = []

    for attr in dir(ccd2clone): 
        if (not callable(getattr(ccd2clone, attr)) and (attr[0:2] != '__')):
            attrs2clone.append(attr)

    for attr in attrs2clone:
        setattr(self, attr, getattr(ccd2clone, attr))

  def get_tiles_stats(self, Quad, tile_coos, statkey=None, estimator=None, 
        extension=-1, binfactor=1, clipsigma=-1):
    """Returns statistics on a list of tiles.

    :param Quad: Quadrant where to take the cutouts from.
    :type Quad: str

    :param tile_coos: A dictionary with tiles coordinates, as output by 
        get_tile_coos.
    :type tile_coos: dict()

    :param statkey: stat to retrieve (one of mean, median, std)
    :type statkey: str

    :param estimator: function to retrieve stat (alternative to statkey)
    :type estimator: obj

    :param extension: Extension index, last is -1
    :type extension: int

    :param binfactor: apply binning of binfactor x binfactor before extracting stat
    :type binfactor: int

    :return: A 1D numpy array with the stat values for the tiles. 

    """

    if isinstance(self.extensions[extension].data, np.ma.masked_array):
        stat_dict = dict(
            mean=np.ma.mean, median=np.ma.median, std=np.ma.std)
    else:
        stat_dict = dict(mean=np.mean, median=np.median, std=np.std)

    if estimator is None and statkey is not None:
        estimator = stat_dict[statkey]

    assert (binfactor>1) # if not, what are you doing, you dummy?
    assert isinstance(binfactor,int)


    tiles = self.get_tiles(Quad, tile_coos, extension=extension)

    if binfactor > 1:

      _tiles = copy.deepcopy(tiles)

      newshape = tuple([item//binfactor for item in tiles[0].shape])
      window = tuple([item*binfactor for item in newshape]) 

      tiles = []
      for i in range(len(_tiles)):
        tiles.append(ccd_aux.rebin(_tiles[i][0:window[0],0:window[1]],newshape,stat='mean'))

    
    if clipsigma <0:
      vals = np.array(list(map(estimator, tiles)))
    else:
      vals = []
      for i in range(len(tiles)):
        ivals = stats.sigmaclip(tiles[i].flatten(), clipsigma, clipsigma).clipped
        vals.append(ivals)
      vals = np.array(vals)
    
    return vals



class HER_CDP(cdp.Tables_CDP):

    def ingest_inputs(self, mx_dct, CCDs, Quads, meta=None, header=None, figs=None):

        keys = list(mx_dct.keys())

        _data = dict()

        for key in keys:

            cbe = mx_dct[key].copy()
            cbe_dict = OrderedDict()

            for jCCD, CCDk in enumerate(CCDs):
                cbe_dict[CCDk] = OrderedDict()
                for kQ, Q in enumerate(Quads):
                    cbe_dict[CCDk][Q] = cbe[jCCD, kQ]

            df = pd.DataFrame.from_dict(cbe_dict)
            _data[key] = df

        super(HER_CDP, self).ingest_inputs(
            _data, meta=meta, header=header, figs=figs)


def get_CDP_lib(test):
    gain_tb_cdp = cdp.Tables_CDP()
    gain_tb_cdp.rootname = '%s_GAIN_TB' % test

    HER_cdp = HER_CDP()
    HER_cdp.rootname = '%s_HER' % test

    HER_profiles_cdp = cdp.CDP()
    HER_profiles_cdp.rootname = 'HER_profiles_%s' % test

    ptccurves_cdp = cdp.CDP()
    ptccurves_cdp.rootname = 'PTC_curves_%s' % test

    bloomMaps_cdp = cdp.FitsTables_CDP()
    bloomMaps_cdp.rootname = 'BloomMaps_%s' % test

    CDP_lib = dict(GAIN_TB=gain_tb_cdp,
                   CURVES_PTC=ptccurves_cdp,
                   HER=HER_cdp,
                   HER_PROFILES=HER_profiles_cdp,
                   BLOOM_MAPS=bloomMaps_cdp)
    return CDP_lib


def gt_PTC_curves_dict(test):

    nicetest = test.replace('_', '\_')

    return dict(
        figname='%s_PTC_curves.png' % test,
        caption='%s: PTC curves and best fits.' % nicetest,
        meta=dict(doLegend=True,
                  ylabel='VAR',
                  xlabel='MED',
                  xlim=[0., 2**16],
                  ylim=[0., 2.**16 / 3.],
                  suptitle='%s: PTC Curves.' % nicetest,
                  corekwargs=dict(data=dict(marker='.', linestyle='', color='b'),
                                  fit=dict(marker='', linestyle='-', color='r'),
                                  bloom=dict(marker='', linestyle='--', color='k')))
    )


def gt_check_offsets_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='%s_offset_vs_time.png' % (test,),
                caption='%s: offset vs. time.' % (ntest,),
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest,
                          ylim=trends.offset_lims))


def gt_check_deltaoff_dict(test):
    ntest = test.replace('_', '\_')
    return dict(
        stats=[
            'deltaoff_pre', 'deltaoff_ove'], trendaxis='time', figname='%s_deltaoff_vs_time.png' %
        (test,), caption='%s: $\delta$offset vs. time. Offset value in each frame minus the average value.' %
        (ntest,), meta=dict(
            doLegend=True, doNiceXDate=True, suptitle='%s-checks: delta-offsets' %
            ntest, ylim=[
                -10., 10.]))


def gt_check_std_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='%s_std_vs_time.png' % test,
                caption='%s: std vs. time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: std' % ntest,
                          ylim=trends.RON_lims))


def gt_check_img_flu_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['flu_med_img'],
                trendaxis='exptime',
                figname='%s_flu_vs_exptime.png' % test,
                caption='%s: Fluence vs. exposure time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=False,
                          suptitle='%s-checks: Fluence' % ntest,
                          xlabel='exptime[s]',
                          ylabel='[ADU]'))


def gt_check_img_std_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['flu_std_img'],
                trendaxis='exptime',
                figname='%s_imgstd_vs_exptime.png' % test,
                caption='%s: STD vs. exposure time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=False,
                          suptitle='%s-checks: Image Area STD' % ntest,
                          xlabel='exptime[s]',
                          ylabel='[ADU]'))


def gt_HER_dict(test):
    ntest = test.replace('_', '\_')
    return dict(
        figname='%s_profs_HER_ser.png' % test,
        caption='%s: HER profiles, serial direction.' % ntest,
        meta=dict(doLegend=True,
                  ylabel=r'$\delta ADU/Cliff Value$',
                  xlabel='Row [pix]',
                  ylim=[-0.002, 0.005],
                  suptitle='%s: Serial HER.' % ntest,
                  corekwargs=dict(marker='.', linestyle='-')))

def gt_BM_dict(test):
  ntest = test.replace('_', '\_')
  return dict(
        figname='%s_BLOOMINGMAP_2D.png' % test,
        caption='%s: Blooming Map of the CCDs. Regions that did not reach blooming are in white.' % ntest,
        meta=dict(doLegend=False,
                  doColorbar=True,
                  suptitle='%s: Blooming Map [ADU].' % ntest,
                  corekwargs=dict(cmap=cm.rainbow, aspect='auto',  norm=None,
                                  origin='lower left')))

def gt_PTC0Xfigs(test):
    PTC0Xfigs = dict()
    PTC0Xfigs['PTC0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, gt_check_offsets_dict(test)]
    PTC0Xfigs['PTC0Xchecks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, gt_check_deltaoff_dict(test)]
    PTC0Xfigs['PTC0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, gt_check_std_dict(test)]
    PTC0Xfigs['PTC0Xchecks_flu'] = [
        trends.Fig_Basic_Checkstat, gt_check_img_flu_dict(test)]
    PTC0Xfigs['PTC0Xchecks_imgstd'] = [
        trends.Fig_Basic_Checkstat, gt_check_img_std_dict(test)]
    PTC0Xfigs['PTC0X_PTC_curves'] = [
        figclasses.Fig_Beam2DPlot, gt_PTC_curves_dict(test)]
    PTC0Xfigs['PTC0X_HER'] = [figclasses.Fig_Beam2DPlot, gt_HER_dict(test)]
    PTC0Xfigs['BLOOM_MAPS'] = [figclasses.Fig_BeamImgShow, gt_BM_dict(test)]
    PTC0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return PTC0Xfigs
