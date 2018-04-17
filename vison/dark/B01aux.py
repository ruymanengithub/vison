#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BIAS01.

Created on Tue Nov 14 13:54:34 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.plot import figclasses
from vison.plot import trends
from vison.datamodel import cdp

# END IMPORT


check_offsets_dict = dict(stats=['offset_pre', 'offset_img', 'offset_ove'],
                          trendaxis='time',
                          figname='BIAS01_offset_vs_time.png',
                          caption='BIAS01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='BIAS01-checks: offsets'))

check_std_dict = dict(stats=['std_pre', 'std_img', 'std_ove'],
                      trendaxis='time',
                      figname='BIAS01_std_vs_time.png',
                      caption='BIAS01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='BIAS01-checks: std')
                      )

basic_prof1Dhor_dict = dict(
        figname='BIAS01_profs1D_hor_allOBSIDs.png',
        caption='BIAS01: Average profiles across columns.',
        meta=dict(doLegend=False,
                  ylabel='ADU',
                  xlabel='Column [pix]',
                  suptitle='BIAS01: Profiles across columns.')
        )

basic_prof1Dver_dict = dict(
        figname='BIAS01_profs1D_ver_allOBSIDs.png',
        caption='BIAS01: Average profiles across rows.',
        meta=dict(doLegend=False,
                  ylabel='ADU',
                  xlabel='Row [pix]',
                  suptitle='BIAS01: Profiles across rows.')
        )

basic_histosRON_dict = dict(
        figname='BIAS01_RON_distro_allOBSIDs.png',
        caption='BIAS01: RON distribution',
        meta = dict(doLegend=False,
                ylabel='N',
                xlabel='RON [ADU]',
                suptitle='BIAS01: RON Distribution'),
        )


B01figs = dict()
B01figs['B01checks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
#B01figs['B01checks_stds'] = [plB01check,dict(stat='std')]
B01figs['B01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
B01figs['B01basic_prof1D_hor'] = [figclasses.Fig_Beam2DPlot, basic_prof1Dhor_dict]
B01figs['B01basic_prof1D_ver'] = [figclasses.Fig_Beam2DPlot, basic_prof1Dver_dict]
B01figs['B01basic_histosRON'] = [figclasses.Fig_Beam1DHist, basic_histosRON_dict]
B01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]


RON_CDP = cdp.Tables_CDP()
RON_CDP.rootname = 'RON_BIAS01'

CDP_lib = dict(RON=RON_CDP)
