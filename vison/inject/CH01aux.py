#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to CHINJ01.

Created on Tue Jan 30 14:09:00 2018

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
                          figname='CHINJ01_offset_vs_time.png',
                          caption='CHINJ01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='CHINJ01-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_ove'],
                      trendaxis='time',
                      figname='CHINJ01_std_vs_time.png',
                      caption='CHINJ01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='CHINJ01-checks: std',
                                ylim=trends.RON_lims)
                      )
check_injlevel_dict = dict(stats=['chk_mea_inject', 'chk_med_inject'],
                           trendaxis='time',
                           figname='CHINJ01_injlevel_vs_time.png',
                           caption='CHINJ01: Injection Level vs. time.',
                           meta=dict(doLegend=True,
                                     doNiceXDate=True,
                                     suptitle='CHINJ01-checks: Injection Level')
                           )

check_injstd_dict = dict(stats=['chk_std_inject'],
                         trendaxis='time',
                         figname='CHINJ01_injstd_vs_time.png',
                         caption='CHINJ01: Injection STD vs. time.',
                         meta=dict(doLegend=False,
                                   doNiceXDate=True,
                                   suptitle='CHINJ01-checks: Injection STD')
                         )

extract_alcol_dict = dict(
    figname='CHINJ01_along_columns_profiles.png',
    caption='CHINJ01: Average along-columns profiles.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='ROW [REL.]',
              ylim=[0,2**16],
              suptitle='CHINJ01: ALONG-COLUMNs.')
    )

extract_alrow_dict = dict(
    figname='CHINJ01_along_rows_profiles.png',
    caption='CHINJ01: Average along-rows profiles.',
    meta=dict(doLegend=True,
              ylabel='ADU',
              xlabel='COL. [REL.]',
              ylim=[0,2**16],
              suptitle='CHINJ01: ALONG-ROWs.')
    )

chinj01_meta_dict = dict(
    figname='CHINJ01_Injection_profiles.png',
    caption='CHINJ01: Injection profiles.',
    meta=dict(doLegend=True,
              ylabel='UNKNOWN',
              xlabel='IG1 [V]',
              suptitle='CHINJ01: INJECTION PROFILE.')
    )

CH01figs = dict()
CH01figs['CH01checks_offsets'] = [
    trends.Fig_Basic_Checkstat, check_offsets_dict]
CH01figs['CH01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
CH01figs['CH01checks_injlevel'] = [
    trends.Fig_Basic_Checkstat, check_injlevel_dict]
CH01figs['CH01checks_injstd'] = [trends.Fig_Basic_Checkstat, check_injstd_dict]
CH01figs['CH01_alrow'] = [
        figclasses.Fig_Beam2DPlot, extract_alrow_dict]
CH01figs['CH01_alcol'] = [
        figclasses.Fig_Beam2DPlot, extract_alcol_dict]
CH01figs['CH01_meta'] = [
        figclasses.Fig_Beam2DPlot, chinj01_meta_dict]
CH01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]


extract_cdp = cdp.Tables_CDP()
extract_cdp.rootname = 'CHINJ01_EXTRACTION_TABLE'

meta_cdp = cdp.Tables_CDP()
meta_cdp.rootname = 'CHINJ01_META_TABLE'

CDP_lib = dict(EXTRACT=extract_cdp,
               META=meta_cdp)
