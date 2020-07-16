#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BF01.

Created on Tue Jul 31 17:50:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import pandas as pd

from vison.flat import PTC0Xaux
from vison.plot import figclasses
from vison.plot import trends

from vison.datamodel import cdp

from vison.point import spot as spotmod

# END IMPORT


def get_CDP_lib():

    covtable_cdp = cdp.Tables_CDP()
    covtable_cdp.rootname = 'BF01_COVTABLE'

    bftable_cdp = cdp.Tables_CDP()
    bftable_cdp.rootname = 'BF01_G15TABLE'

    bfFITtable_cdp = cdp.Tables_CDP()
    bfFITtable_cdp.rootname = 'BF01FIT_G15TABLE'

    profscov_cdp = cdp.CDP()
    profscov_cdp.rootname = 'profs_COV1D_BF01'

    profsker_cdp = cdp.CDP()
    profsker_cdp.rootname = 'profs_KER1D_BF01'

    CDP_lib = dict(COVTABLE=covtable_cdp,
                   PROFSCOV1D=profscov_cdp,
                   BFTABLE=bftable_cdp,
                   BFfitTABLE=bfFITtable_cdp,
                   PROFSKER1D=profsker_cdp)
    return CDP_lib


prof_COV_ver_dict = dict(
    figname='BF01_COV_profs_ver.png',
    caption='BF01: COV 1D profiles, vertical/parallel direction.',
    meta=dict(doLegend=True,
              ylabel='COV/VAR, ADIM.',
              xlabel='Y',
              ylim=[-0.005,0.07],
              corekwargs=dict(),
              suptitle='BF01: COV 1D Profile, Vertical/Parallel.')
)

prof_COV_ser_dict = dict(
    figname='BF01_COV_profs_ser.png',
    caption='BF01: COV 1D profiles, serial direction.',
    meta=dict(doLegend=True,
              ylabel='COV/VAR, ADIM.',
              xlabel='X',
              ylim=[-0.005,0.07],
              corekwargs=dict(),
              suptitle='BF01: COV 1D Profile, Serial.')
)

prof_KER_ver_dict = dict(
    figname='BF01_KER_profs_ver.png',
    caption='BF01: KERNEL 1D profiles, vertical/parallel direction.',
    meta=dict(doLegend=True,
              ylabel='log([ADU])',
              xlabel='Y [pix]',
              #ylim=[-0.03, 0.1],
              corekwargs=dict(),
              suptitle='BF01: Kernels 1D Profile, Vertical/Parallel.')
)

prof_KER_ser_dict = dict(
    figname='BF01_KER_profs_ser.png',
    caption='BF01: KERNEL 1D profiles, serial direction.',
    meta=dict(doLegend=True,
              ylabel='log([ADU])',
              xlabel='X [pix]',
              #ylim=[-0.03, 0.1],
              corekwargs=dict(),
              suptitle='BF01: Kernels 1D Profile, Serial.')
)

FWHMx_v_flu_dict = dict(
    figname='BF01_FWHMx_v_flu.png',
    caption='BF01: FWHM(x) vs. Fluence.',
    meta=dict(doLegend=True,
              ylabel='FWHM(x), [um]',
              xlabel='ADU',
              ylim = [5.,15.],
              corekwargs=dict(
                  data=dict(marker='o', linestyle=''),
                  fit=dict(marker='', linestyle='--')),
              suptitle='BF01: FWHMx in microns vs. Fluence')
)

FWHMy_v_flu_dict = dict(
    figname='BF01_FWHMy_v_flu.png',
    caption='BF01: FWHM(y) vs. Fluence.',
    meta=dict(doLegend=True,
              ylabel='FWHM(y), [um]',
              xlabel='ADU',
              ylim = [5.,15.],
              corekwargs=dict(
                  data=dict(marker='o', linestyle=''),
                  fit=dict(marker='', linestyle='--')),
              suptitle='BF01: FWHMy in microns vs. Fluence')
)


def gt_PTC_curves_dict(test, BFEcorr='no'):

  nicetest = st.replace(test, '_', '\_')

  if BFEcorr=='no':
    BFEtag = 'wBFE'
    BFEcap = 'BFE left in'
  elif BFEcorr == 'fludep':
    BFEtag = 'noBFE'
    BFEcap = 'BFE removed'
  elif BFEcorr == 'flufix':
    BFEtag = 'noBFEfixflu'
    BFEcap = 'BFE removed, mid-fluence Axij'

  
  return dict(
      figname='%s_PTC_curves_%s.png' % (test,BFEtag),
      caption='%s: PTC curves, %s. Theoretical line has a fixed gain of 3.5.' % (nicetest,BFEcap),
      meta=dict(doLegend=True,
              ylabel='VAR',
              xlabel='MED',
              xlim=[0., 2**16],
              ylim=[0., 2.**16 / 3.],
              suptitle='%s: PTC Curves.' % nicetest,
              corekwargs=dict(data=dict(marker='.', linestyle='', color='b'),
                      theo=dict(marker='', linestyle='-', color='r'))
        ))


def gt_BF01figs(test):

    BF01figs = dict()
    BF01figs['BF01checks_offsets'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_offsets_dict(test)]
    BF01figs['BF01checks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_deltaoff_dict(test)]
    BF01figs['BF01checks_stds'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_std_dict(test)]
    BF01figs['BF01checks_flu'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_img_flu_dict(test)]
    BF01figs['BF01checks_imgstd'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_img_std_dict(test)]
    BF01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    # renaming to BF01... HACK

    keys_to_rename = ['caption', 'suptitle', 'figname']
    for figkey in list(BF01figs.keys()):
        _dict = BF01figs[figkey][1]
        try:
            _mdict = BF01figs[figkey][1]['meta']
            hasmeta = True
        except KeyError:
            hasmeta = False
        for key in keys_to_rename:
            if key in _dict:
                _dict[key] = st.replace(_dict[key], test, 'BF01')
            if hasmeta:
                if key in _mdict:
                    _mdict[key] = st.replace(_mdict[key], test, 'BF01')

    BF01figs['BF01_COV_ver'] = [
        figclasses.Fig_Beam2DPlot, prof_COV_ver_dict]

    BF01figs['BF01_COV_hor'] = [
        figclasses.Fig_Beam2DPlot, prof_COV_ser_dict]

    BF01figs['BF01_KER_ver'] = [
        figclasses.Fig_Beam2DPlot, prof_KER_ver_dict]

    BF01figs['BF01_KER_hor'] = [
        figclasses.Fig_Beam2DPlot, prof_KER_ser_dict]

    BF01figs['BF01_fwhmx_v_flu'] = [
        figclasses.Fig_Beam2DPlot, FWHMx_v_flu_dict]

    BF01figs['BF01_fwhmy_v_flu'] = [
        figclasses.Fig_Beam2DPlot, FWHMy_v_flu_dict]

    BF01figs['BF01_PTC_BFE'] = [
        figclasses.Fig_Beam2DPlot, gt_PTC_curves_dict(test, BFEcorr='no')]

    BF01figs['BF01_PTC_NOBFE'] = [
        figclasses.Fig_Beam2DPlot, gt_PTC_curves_dict(test, BFEcorr='fludep')]

    BF01figs['BF01_PTC_NOBFEALT'] = [
        figclasses.Fig_Beam2DPlot, gt_PTC_curves_dict(test, BFEcorr='flufix')]



    return BF01figs

