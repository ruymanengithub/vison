#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to TP01.

Created on Wed Jan 31 15:11:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
from matplotlib import cm
import pandas as pd

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT

CCDs = ['CCD1','CCD2','CCD3']

check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='TP01_offset_vs_time.png',
                          caption='TP01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='TP01-checks: offsets',
                                    ylim=trends.offset_lims))

check_deltaoff_dict = dict(stats=['deltaoff_pre', 'deltaoff_ove'],
                          trendaxis='time',
                          figname='TP01_deltaoff_vs_time.png',
                          caption='TP01: $\delta$offset vs. time. Offset value in each frame minus the average value.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='TP01-checks: delta-offset',
                                    ylim=[-10.,10.]))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='TP01_std_vs_time.png',
                      caption='TP01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='TP01-checks: std',
                                ylim=trends.RON_lims))

check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='TP01_injlevel_vs_time.png',
                           caption='TP01: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='TP01-checks: Injection Level'))

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='TP01_injstd_vs_time.png',
                         caption='TP01: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='TP01-checks: Injection STD')
                         )

def gt_meta_heatmaps(mkey):
    return dict(
            figname='TP01_PcTau_Heatmap_%s.png' % mkey,
            caption='TP01: log10(tau[us])-log10(Pc) Heatmap for pumping mode %s.' % mkey,
            meta=dict(doLegend=False,
              doColorbar=False,
              suptitle='TP01: %s' % mkey,
              xlabel='log10(tau [us])',
              ylabel='log10(Pc)',
              corekwargs=dict(cmap=cm.inferno_r,
                               aspect='auto',
                               norm=None,
                               origin='lower left',
                               extent=[]))
                )
    
def get_TP01figs():
    TP01figs = dict()
    TP01figs['TP01checks_offsets'] = [
        trends.Fig_Basic_Checkstat, check_offsets_dict]
    TP01figs['TP01checks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, check_deltaoff_dict]
    TP01figs['TP01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    TP01figs['TP01checks_injlevel'] = [
        trends.Fig_Basic_Checkstat, check_injlevel_dict]
    TP01figs['TP01checks_injstd'] = [trends.Fig_Basic_Checkstat, check_injstd_dict]
    
    
    for mkey in ['m123','m234','m341','m412']:
        
        TP01figs['TP01meta_%s' % mkey] = [figclasses.Fig_BeamImgShow,
                gt_meta_heatmaps(mkey)]
    
    
    TP01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return TP01figs



def get_CDP_lib():
    """ """
    
    CDP_lib = OrderedDict()
    
    for CCD in CCDs:
    
        mastercat = cdp.CDP()
        mastercat.rootname = 'TP01_MasterCat_%s' % CCD
        CDP_lib['MASTERCAT_%s' % CCD] = mastercat
               
        mergedcat = cdp.CDP()
        mergedcat.rootname = 'TP01_MergedCat_%s' % CCD
        CDP_lib['MERGEDCAT_%s' % CCD] = mergedcat
        
               
        chinjcharact_cdp = cdp.CDP()
        chinjcharact_cdp.rootname = 'TP01_CHINJCHAR'
        CDP_lib['CHINJCHARACT'] = chinjcharact_cdp
        
        chinj_cdp = cdp.Tables_CDP()
        chinj_cdp.rootname = 'TP01_CHINJ'
        
        CDP_lib['CHINJ'] = chinj_cdp
    
    return CDP_lib
    