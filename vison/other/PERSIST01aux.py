#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxilisary Functions and resources to PERSIST01.


Created on Thu Feb  7 18:28:47 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd
from matplotlib import cm
import copy

from vison.plot import figclasses
from vison.plot import trends
from vison.datamodel import cdp


# END IMPORT

check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='PERSIST01_offset_vs_time.png',
                          caption='PERSIST01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='PERSIST01-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='PERSIST01_std_vs_time.png',
                      caption='PERSIST01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='PERSIST01-checks: std',
                                ylim=trends.RON_lims)
                      )

check_avgflu_dict = dict(stats=['chk_avgflu_img'],
                      trendaxis='time',
                      figname='PERSIST01_avgbgd_vs_time.png',
                      caption='PERSIST01: Background level vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='PERSIST01-checks: background',
                                )
                      )

P01satmask_dict = dict(
        figname='PERSIST01_SATMASK_2Dimgshow.png',
            caption='PERSIST01: Saturation Masks for the CCDs.',
            meta=dict(doLegend=False,
              doColorbar=True,
              suptitle='PERSIST01: Saturation Masks',
              corekwargs=dict(cmap=cm.gray,aspect='auto',norm=None,
                                origin='lower left'))
        )


def get_P01figs():
    P01figs = dict()
    P01figs['P01checks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
    #P01figs['B01checks_stds'] = [plB01check,dict(stat='std')]
    P01figs['P01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    P01figs['P01checks_avgflu'] = [trends.Fig_Basic_Checkstat, check_avgflu_dict]
    P01figs['P01satmasks'] = [figclasses.Fig_BeamImgShow, P01satmask_dict]
    P01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return P01figs


def get_CDP_lib():
    
    satarea_cdp = cdp.Tables_CDP()
    satarea_cdp.rootname = 'PERSIST01_SATAREAS'
    
    CDP_lib = OrderedDict()    
    CDP_lib['SATAREA_TB'] = satarea_cdp
    
    return CDP_lib
    