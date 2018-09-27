#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to FOCUS00.

Created on Tue Nov 14 13:54:34 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends
from vison.point import Paux

# END IMPORT


def gt_check_offsets_dict(wave):
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='FOCUS00_%i_offset_vs_time.png' % wave,
                caption='FOCUS00\_%i: offset vs. time.' % wave,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='FOCUS01\_%i-checks: offsets' % wave,
                          ylim=trends.offset_lims))


def gt_check_std_dict(wave):
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='FOCUS00_%i_std_vs_time.png' % wave,
                caption='FOCUS00\_%i: std vs. time.' % wave,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='FOCUS00\_%i-checks: std' % wave,
                          ylim=trends.RON_lims))


def gt_check_bgd_dict(wave):
    return dict(stats=['bgd_img'],
                trendaxis='time',
                figname='FOCUS00_%i_bgd_vs_time.png' % wave,
                caption='FOCUS00\_%i: Background vs. time.' % wave,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='FOCUS00\_%i-checks: BGD' % wave))


def gt_check_flu_dict(wave):
    return dict(stats=['chk_fluence'],
                trendaxis='mirr_pos',
                figname='FOCUS00_%i_flu_vs_mirr.png' % wave,
                caption='FOCUS00\_%i: Fluence vs. Mirror Position.' % wave,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='FOCUS00\_%i-checks: Fluence' % wave,
                          xlabel='Mirr [mm]',
                          ylabel='Flu.[ADU]'))


def gt_check_fwhmx_dict(wave):
    return dict(stats=['chk_fwhmx'],
                trendaxis='mirr_pos',
                figname='FOCUS00_%i_fwhmx_vs_mirr.png' % wave,
                caption='FOCUS00\_%i: FWHM(x) vs. Mirror Position.' % wave,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='FOCUS00\_%i-checks: FWHM(x)' % wave,
                          xlabel='Mirr [mm]',
                          ylabel='FWHMx [pix]'))


def gt_check_fwhmy_dict(wave):
    return dict(stats=['chk_fwhmy'],
                trendaxis='mirr_pos',
                figname='FOCUS00_%i_fwhmy_vs_mirr.png' % wave,
                caption='FOCUS00\_%i: FWHM(y) vs. Mirror Position.' % wave,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='FOCUS00\_%i-checks: FWHM(y)' % wave,
                          xlabel='Mirr [mm]',
                          ylabel='FWHMy [pix]'))

def gt_meta_deltafwhm_dict(wave):
    return dict(figname='FOCUS00_%i_deltafwhm_map.png' % wave,
                caption='FOCUS00\_%i: Delta-FWHM for best focus position' % wave,
                meta=dict(suptitle='FOCUS00\_%i-checks: FWHM(y)' % wave))

def gt_F00figs(wave):
    F00figs = dict()
    
    # CHECKS
    
    F00figs['F00checks_offsets'] = [
        trends.Fig_Basic_Checkstat, gt_check_offsets_dict(wave)]
    F00figs['F00checks_stds'] = [
        trends.Fig_Basic_Checkstat, gt_check_std_dict(wave)]
    F00figs['F00checks_bgd'] = [
        trends.Fig_Basic_Checkstat, gt_check_bgd_dict(wave)]
    F00figs['F00checks_fluence'] = [
        trends.Fig_Basic_Checkstat, gt_check_flu_dict(wave)]
    F00figs['F00checks_fwhmx'] = [
        trends.Fig_Basic_Checkstat, gt_check_fwhmx_dict(wave)]
    F00figs['F00checks_fwhmy'] = [
        trends.Fig_Basic_Checkstat, gt_check_fwhmy_dict(wave)]
    
    # META
    
    F00figs['F00meta_deltafwhm'] = [
        figclasses.Fig_BeamImgShow, gt_meta_deltafwhm_dict(wave)]
    
    F00figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return F00figs



focus_cdp = cdp.Tables_CDP()
focus_cdp.rootname = 'FOCUS00'

CDP_lib = dict(FOCUS=focus_cdp)
CDP_lib.update(Paux.CDP_lib)
