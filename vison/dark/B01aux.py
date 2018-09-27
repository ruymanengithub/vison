#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BIAS01.

Created on Tue Nov 14 13:54:34 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd

from vison.plot import figclasses
from vison.plot import trends
from vison.datamodel import cdp
from matplotlib import cm

# END IMPORT


check_offsets_dict = dict(stats=['offset_pre', 'offset_img', 'offset_ove'],
                          trendaxis='time',
                          figname='BIAS01_offset_vs_time.png',
                          caption='BIAS01: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='BIAS01-checks: offsets',
                                    ylim=trends.offset_lims))

check_std_dict = dict(stats=['std_pre', 'std_img', 'std_ove'],
                      trendaxis='time',
                      figname='BIAS01_std_vs_time.png',
                      caption='BIAS01: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='BIAS01-checks: std',
                                ylim=trends.RON_lims)
                      )

basic_prof1Dhor_dict = dict(
    figname='BIAS01_profs1D_hor_allOBSIDs.png',
    caption='BIAS01: Average profiles across columns.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Column [pix]',
              ylim = [0., 2.**16],
              suptitle='BIAS01: Profiles across columns.')
)

basic_prof1Dver_dict = dict(
    figname='BIAS01_profs1D_ver_allOBSIDs.png',
    caption='BIAS01: Average profiles across rows.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim = [0., 2.**16],
              suptitle='BIAS01: Profiles across rows.')
)



basic_histosRON_dict = dict(
    figname='BIAS01_RON_distro_allOBSIDs.png',
    caption='BIAS01: RON distribution',
    meta=dict(doLegend=True,
              ylabel='N',
              xlabel='RON [ADU]',
              xlim=[0.,4.],
              suptitle='BIAS01: RON Distribution')
)


meta_prof1Dhor_dict = dict(
    figname='BIAS01_profs1D_hor_MASTERBIAS.png',
    caption='BIAS01: Average profiles across columns of Master Bias.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Column [pix]',
              ylim = [0., 2.**16],
              suptitle='BIAS01/Master: Profiles across columns.')
)

meta_prof1Dver_dict = dict(
    figname='BIAS01_profs1D_ver_MASTERBIAS.png',
    caption='BIAS01: Average profiles across rows of Master Bias.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim = [0., 2.**16],
              suptitle='BIAS01/Master: Profiles across rows.')
)

meta_MB2D_dict = dict(
    figname='BIAS01_MASTERBIAS_2Dimgshow.png',
    caption='BIAS01: Master Bias for the CCDs.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim = [0., 2.**16],
              suptitle='BIAS01/Master: Profiles across rows.'),
    corekwargs=dict(cmap=cm.rainbow,aspect='auto',norm=None,origin='lower left',
                    ))

B01figs = dict()
B01figs['B01checks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
#B01figs['B01checks_stds'] = [plB01check,dict(stat='std')]
B01figs['B01checks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
B01figs['B01basic_prof1D_hor'] = [
    figclasses.Fig_Beam2DPlot, basic_prof1Dhor_dict]
B01figs['B01basic_prof1D_ver'] = [
    figclasses.Fig_Beam2DPlot, basic_prof1Dver_dict]
B01figs['B01basic_histosRON'] = [
    figclasses.Fig_Beam1DHist, basic_histosRON_dict]
B01figs['B01meta_prof1D_hor'] = [
    figclasses.Fig_Beam2DPlot, meta_prof1Dhor_dict]
B01figs['B01meta_prof1D_ver'] = [
    figclasses.Fig_Beam2DPlot, meta_prof1Dver_dict]
B01figs['B01meta_MasterBias_2D'] = [
    figclasses.Fig_BeamImgShow, meta_MB2D_dict]
B01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]


class RON_CDP(cdp.Tables_CDP):

    def ingest_inputs(self, mx_dct, CCDs, Quads, meta=None, header=None, figs=None):
        
        keys = mx_dct.keys()
        
        _data = dict()
        
        for key in keys:
            
            mx = mx_dct[key].copy()

            msk = np.zeros_like(mx, dtype='int32')
            msk[np.where((np.isclose(mx, 0.)) | (np.isnan(mx)))] = 1
            cbe = np.ma.median(np.ma.masked_array(
                mx, mask=msk), axis=1).data.copy()  # best estimate of RON
    
            cbe_dict = OrderedDict()
            for jCCD, CCDk in enumerate(CCDs):
                cbe_dict[CCDk] = OrderedDict()
                for kQ, Q in enumerate(Quads):
                    cbe_dict[CCDk][Q] = cbe[jCCD, kQ]
            df = pd.DataFrame.from_dict(cbe_dict)  # PENDING
    
            _data[key] = df
        
        
        super(RON_CDP, self).ingest_inputs(
            _data, meta=meta, header=header, figs=figs)


ron_cdp = RON_CDP()
ron_cdp.rootname = 'RON_BIAS01'
off_cdp = RON_CDP()
off_cdp.rootname = 'OFFSET_BIAS01'

CDP_lib = dict(RON=ron_cdp,
               OFF=off_cdp)
