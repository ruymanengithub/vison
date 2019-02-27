#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to TP02.

Created on Wed Jan 31 14:59:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT


CCDs = ['CCD1','CCD2','CCD3']

check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='TP02_offset_vs_time.png',
                          caption='TP02: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='TP02-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='TP02_std_vs_time.png',
                      caption='TP02: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='TP02-checks: std',
                                ylim=trends.RON_lims))

check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='TP02_injlevel_vs_time.png',
                           caption='TP02: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='TP02-checks: Injection Level'))

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='TP02_injstd_vs_time.png',
                         caption='TP02: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='TP02-checks: Injection STD')
                         )

def gt_meta_Pctau_dict(mkey):
    return dict(
            figname='TP02_PcTau_%s.png' % mkey,
            caption='TP02: log10(tph[us])-log10(Pc) plot for pumping mode %s. ' % mkey+\
                'tph is the sum of the dwell time and the serial toi, 4.75 us.',
            meta=dict(doLegend=True,
              doNiceXDate=False,
              suptitle='TP02: %s' % mkey,
              xlabel='log10(tph [us])',
              ylabel='log10(Pc)',
              corekwargs=dict(North=dict(marker='.',linestyle='',color='b'),
                              South=dict(marker='.',linestyle='',color='r'))
                ))

def get_TP02figs():
    TP02figs = dict()
    TP02figs['TP02checks_offsets'] = [
        trends.Fig_Basic_Checkstat, check_offsets_dict]
    TP02figs['TP02checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    TP02figs['TP02checks_injlevel'] = [
        trends.Fig_Basic_Checkstat, check_injlevel_dict]
    TP02figs['TP02checks_injstd'] = [trends.Fig_Basic_Checkstat, check_injstd_dict]
    
    
    for mkey in ['m23','m31']:
        
        TP02figs['TP02meta_%s' % mkey] = [figclasses.Fig_Beam2DPlot,
                gt_meta_Pctau_dict(mkey)]
    
    TP02figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return TP02figs




def get_CDP_lib():
    """ """
    
    mastercat = cdp.CDP()
    mastercat.rootname = 'TP02_MasterCat'
    
    CDP_lib = OrderedDict()
           
    for CCD in CCDs:
    
        mastercat = cdp.CDP()
        mastercat.rootname = 'TP01_MasterCat_%s' % CCD
        CDP_lib['MASTERCAT_%s' % CCD] = mastercat
               
        mergedcat = cdp.CDP()
        mergedcat.rootname = 'TP01_MergedCat_%s' % CCD
        CDP_lib['MERGEDCAT_%s' % CCD] = mergedcat
          
    return CDP_lib