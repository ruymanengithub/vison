#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to CHINJ02.

Created on Tue Jan 30 14:28:00 2018

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


check_offsets_dict = dict(stats=['offset_pre', 'offset_ove'],
                          trendaxis='time',
                          figname='CHINJ02_offset_vs_time.png',
                          caption='CHINJ02: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='CHINJ02-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='CHINJ02_std_vs_time.png',
                      caption='CHINJ02: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='CHINJ02-checks: std',
                                ylim=trends.RON_lims)
                      )
check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='CHINJ02_injlevel_vs_time.png',
                           caption='CHINJ02: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='CHINJ02-checks: Injection Level')
                           )

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='CHINJ02_injstd_vs_time.png',
                         caption='CHINJ02: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='CHINJ02-checks: Injection STD')
                         )

extract_alcol_dict = dict(
    figname='CHINJ02_along_columns_profiles.png',
    caption='CHINJ02: Average along-columns profiles.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='ROW [REL.]',
              ylim=[0,2**16],
              suptitle='CHINJ02: ALONG-COLUMNs.')
    )

extract_alrow_dict = dict(
    figname='CHINJ02_along_rows_profiles.png',
    caption='CHINJ02: Average along-rows profiles.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='COL. [REL.]',
              ylim=[0,2**16],
              suptitle='CHINJ02: ALONG-ROWs.')
    )


CH02figs = dict()
CH02figs['CH02checks_offsets'] = [
    trends.Fig_Basic_Checkstat, check_offsets_dict]
CH02figs['CH02checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
CH02figs['CH02checks_injlevel'] = [
    trends.Fig_Basic_Checkstat, check_injlevel_dict]
CH02figs['CH02checks_injstd'] = [trends.Fig_Basic_Checkstat, check_injstd_dict]
CH02figs['CH02_alrow'] = [
        figclasses.Fig_Beam2DPlot, extract_alrow_dict]
CH02figs['CH02_alcol'] = [
        figclasses.Fig_Beam2DPlot, extract_alcol_dict]
CH02figs['BlueScreen'] = [figclasses.BlueScreen, dict()]

extract_cdp = cdp.Tables_CDP()
extract_cdp.rootname = 'CHINJ02_EXTRACTION_TABLE'

meta_cdp = cdp.Tables_CDP()
meta_cdp.rootname = 'CHINJ01_META_TABLE'


CDP_lib = dict(EXTRACT=extract_cdp,
               META=meta_cdp)
