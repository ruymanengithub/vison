#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to DARK01.

Created on Sat Jan 27 16:11:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
from matplotlib import cm

from vison.plot import figclasses
from vison.plot import trends
from vison.datamodel import cdp

# END IMPORT


def get_CDP_lib():
    
    dk_tb_cdp = cdp.Tables_CDP()
    dk_tb_cdp.rootname = 'DARK01_TB'
    
    MB_profiles_cdp = cdp.CDP()
    MB_profiles_cdp.rootname = 'MB_profiles_DARK01'
    
    CDP_lib = dict(DARK_TB=dk_tb_cdp,
                   MB_profiles=MB_profiles_cdp)
    return CDP_lib

check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='DARK01_offset_vs_time.png',
                          caption='DARK01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='DARK01-checks: offsets',
                                    ylim=trends.offset_lims))

check_deltaoff_dict = dict(stats=['deltaoff_pre', 'deltaoff_ove'],
                          trendaxis='time',
                          figname='DARK01_deltaoff_vs_time.png',
                          caption='DARK01: offset-<offset> vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='DARK01-checks: delta-offsets',
                                    ylim=[-10.,10.]))


check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='DARK01_std_vs_time.png',
                      caption='DARK01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='DARK01-checks: std',
                                ylim=trends.RON_lims)
                      )

check_flu_dict = dict(stats=['chk_flu_img'],
                      trendaxis='time',
                      figname='DARK01_FLU_vs_time.png',
                      caption='DARK01: Fluence vs. time.',
                      meta=dict(doLegend=False,
                                doNiceXDate=True,
                                suptitle='DARK01-checks: Image Area Fluence.',
                                ylim=[-10.,20.])
                      )

meta_prof1Dhor_dict = dict(
    figname='DARK01_profs1D_hor_MASTERDARK.png',
    caption='DARK01: Average profiles across columns of Master Dark.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Column [pix]',
              ylim = [-20., 20.],
              suptitle='DARK01/Master: Profiles across columns.')
)

meta_prof1Dver_dict = dict(
    figname='DARK01_profs1D_ver_MASTERDARK.png',
    caption='BIAS01: Average profiles across rows of Master Dark.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim = [-20., 100.],
              suptitle='DARK01/Master: Profiles across rows.')
)

meta_MDK2D_dict = dict(
    figname='DARK01_MASTERDARK_2Dimgshow.png',
    caption='DARK01: Master Dark for the CCDs.',
    meta=dict(doLegend=False,
              doColorbar=True,
              suptitle='DARK01/Master:Quadrant Images',
    corekwargs=dict(cmap=cm.rainbow,aspect='auto',norm=None,origin='lower left',
                    )))


def get_D01figs():
    D01figs = dict()
    D01figs['D01checks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
    D01figs['D01checks_deltaoff'] = [trends.Fig_Basic_Checkstat, check_deltaoff_dict]
    D01figs['D01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    D01figs['D01checks_flu'] = [trends.Fig_Basic_Checkstat, check_flu_dict]
    D01figs['D01meta_prof1D_hor'] = [
            figclasses.Fig_Beam2DPlot, meta_prof1Dhor_dict]
    D01figs['D01meta_prof1D_ver'] = [
            figclasses.Fig_Beam2DPlot, meta_prof1Dver_dict]
    D01figs['D01meta_MasterDark_2D'] = [
            figclasses.Fig_BeamImgShow, meta_MDK2D_dict]
    D01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return D01figs


