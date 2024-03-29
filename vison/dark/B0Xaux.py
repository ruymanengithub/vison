#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BIAS0X.

Created on Tue Nov 14 13:54:34 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd
from matplotlib import cm

from vison.plot import figclasses
from vison.plot import trends
from vison.datamodel import cdp


# END IMPORT


check_offsets_dict = dict(stats=['offset_pre', 'offset_img', 'offset_ove'],
                          trendaxis='time',
                          figname='BIAS0X_offset_vs_time.png',
                          caption='BIAS0X: offset vs. time.',
                          meta=dict(doLegend=True,
                                    doNiceXDate=True,
                                    suptitle='BIAS0X-checks: offsets',
                                    ylim=trends.offset_lims))

check_deltaoff_dict = dict(
    stats=[
        'deltaoff_pre',
        'deltaoff_img',
        'deltaoff_ove'],
    trendaxis='time',
    figname='BIAS0X_deltaoff_vs_time.png',
    caption='BIAS0X: $\delta$offset vs. time. Offset value in each frame minus the average value.',
    meta=dict(
            doLegend=True,
            doNiceXDate=True,
            suptitle='BIAS0X-checks: delta-offsets',
            ylim=[
                -10.,
                10.]))


check_std_dict = dict(stats=['std_pre', 'std_img', 'std_ove'],
                      trendaxis='time',
                      figname='BIAS0X_std_vs_time.png',
                      caption='BIAS0X: std vs. time.',
                      meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='BIAS0X-checks: std',
                                ylim=trends.RON_lims)
                      )

basic_prof1Dhor_dict = dict(
    figname='BIAS0X_profs1D_hor_allOBSIDs.png',
    caption='BIAS0X: Average profiles across columns, avg. offset subtracted. '+\
            'Output node on the left.' +\
            'From col. "1" to col. 51+2048+20="2119" and from vstart to vend rows.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Column [pix]',
              ylim=[-20., 20.],
              suptitle='BIAS0X: Profiles across columns.')
)

basic_prof1Dver_dict = dict(
    figname='BIAS0X_profs1D_ver_allOBSIDs.png',
    caption='BIAS0X: Average profiles across rows, avg. offset subtracted.' +\
            'Output node on the left.'+\
            ' From vstart to vend rows and from col. "1" to col. 51+2048+20="2119".',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim=[-20., 20.],
              suptitle='BIAS0X: Profiles across rows.')
)

basic_prof1Dstdver_dict = dict(
    figname='BIAS0X_profs1Dstd_ver_allOBSIDs.png',
    caption='BIAS0X: Average STDDEV profiles across rows.' +
            ' From vstart to vend rows, and from col. "52" to col. "2100".',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim=[0.5, 2.5],
              suptitle='BIAS0X: STDDEV Profiles across rows.')
)

ronrange = [0.5, 2]

basic_histosRON_dict = dict(
    figname='BIAS0X_RON_distro_allOBSIDs.png',
    caption='BIAS0X: RON distribution',
    meta=dict(doLegend=True,
              ylabel='N',
              xlabel='RON [ADU]',
              xlim=ronrange,
              suptitle='BIAS0X: RON Distribution',
              corekwargs=dict())
)


meta_prof1Dhor_dict = dict(
    figname='BIAS0X_profs1D_hor_MASTERBIAS.png',
    caption='BIAS0X: Average profiles across columns of Master Bias, avg. offset subtracted. ' +
            'Output node on the left.'+\
            ' From col. "1" to col. 51+2048+20="2119", and from vstart to vend rows.',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Column [pix]',
              ylim=[-20., 20.],
              suptitle='BIAS0X/Master: Profiles across columns.')
)

meta_prof1Dver_dict = dict(
    figname='BIAS0X_profs1D_ver_MASTERBIAS.png',
    caption='BIAS0X: Average profiles across rows of Master Bias, avg. offset subtracted. ' +
            'Output node on the left.'+\
            ' From vstart to vend rows and from col. "1" to col. 51+2048+20="2119".',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim=[-20., 20.],
              suptitle='BIAS0X/Master: Profiles across rows.')
)

meta_prof1Dstdver_dict = dict(
    figname='BIAS0X_profs1Dstd_ver_MASTERBIAS.png',
    caption='BIAS0X: Average profiles of STDDEV across rows of Master Bias.' +
            ' From vstart to vend rows, and from col. "52" to col. "2100".',
    meta=dict(doLegend=False,
              ylabel='ADU',
              xlabel='Row [pix]',
              ylim=[0.0, 2.5],
              suptitle='BIAS0X/Master: STDDEV Profiles across rows.')
)

meta_MB2D_dict = dict(
    figname='BIAS0X_MASTERBIAS_2Dimgshow.png',
    caption='BIAS0X: Master Bias for the CCDs.',
    meta=dict(doLegend=False,
              doColorbar=True,
              suptitle='BIAS0X/Master:Quadrant Images',
              corekwargs=dict(cmap=cm.gray, aspect='auto', norm=None,
                              origin='lower left',
                              ))
)


meta_std_vs_mean_dict = dict(
        figname='BIAS0X_std_vs_offset_prescan.png',
        caption='Plot of standard deviation vs. mean in each column of the images, in the pre-scan region. '+\
        'Each colour corresponds to a different frame / acquisition. In each quadrant, the "means" '+\
        'have the offset value of the first image subtracted for ease of representation.',
        meta=dict(doLegend=False,
              ylabel='STD [ADU]',
              xlabel='MEAN [ADU]',
              xlim=[-2., 10.],
              suptitle='BIAS0X: STD vs. OFFSET, prescan')
        )

def get_B0Xfigs():
    B0Xfigs = dict()
    B0Xfigs['B0Xchecks_offsets'] = [trends.Fig_Basic_Checkstat, check_offsets_dict]
    B0Xfigs['B0Xchecks_deltaoff'] = [trends.Fig_Basic_Checkstat, check_deltaoff_dict]
    B0Xfigs['B0Xchecks_stds'] = [trends.Fig_Basic_Checkstat, check_std_dict]
    B0Xfigs['B0Xbasic_prof1D_hor'] = [
        figclasses.Fig_Beam2DPlot, basic_prof1Dhor_dict]
    B0Xfigs['B0Xbasic_prof1D_ver'] = [
        figclasses.Fig_Beam2DPlot, basic_prof1Dver_dict]
    B0Xfigs['B0Xbasic_prof1Dstd_ver'] = [
        figclasses.Fig_Beam2DPlot, basic_prof1Dstdver_dict]
    B0Xfigs['B0Xbasic_histosRON'] = [
        figclasses.Fig_Beam1DHist, basic_histosRON_dict]
    B0Xfigs['B0Xmeta_prof1D_hor'] = [
        figclasses.Fig_Beam2DPlot, meta_prof1Dhor_dict]
    B0Xfigs['B0Xmeta_prof1D_ver'] = [
        figclasses.Fig_Beam2DPlot, meta_prof1Dver_dict]
    B0Xfigs['B0Xmeta_prof1Dstd_ver'] = [
        figclasses.Fig_Beam2DPlot, meta_prof1Dstdver_dict]
    B0Xfigs['B0Xmeta_MasterBias_2D'] = [
        figclasses.Fig_BeamImgShow, meta_MB2D_dict]
    B0Xfigs['B0Xmeta_std_vs_mean'] = [
        figclasses.Fig_Beam2DPlot, meta_std_vs_mean_dict]
    B0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return B0Xfigs


class RON_CDP(cdp.Tables_CDP):

    def ingest_inputs(self, mx_dct, CCDs, Quads, meta=None, header=None, figs=None):

        keys = list(mx_dct.keys())

        _data = dict()

        for key in keys:

            mx = mx_dct[key].copy()

            msk = np.zeros_like(mx, dtype='int32')
            msk[np.where((np.isclose(mx, 0.)) | (np.isnan(mx)))] = 1
            cbe = np.ma.median(np.ma.masked_array(
                mx, mask=msk), axis=0).data.copy()  # best estimate of RON

            cbe_dict = OrderedDict()
            for jCCD, CCDk in enumerate(CCDs):
                cbe_dict[CCDk] = OrderedDict()
                for kQ, Q in enumerate(Quads):
                    cbe_dict[CCDk][Q] = cbe[jCCD, kQ]
            df = pd.DataFrame.from_dict(cbe_dict)  # PENDING

            _data[key] = df

        super(RON_CDP, self).ingest_inputs(
            _data, meta=meta, header=header, figs=figs)


def get_CDP_lib():
    """ """
    ron_cdp = RON_CDP()
    ron_cdp.rootname = 'RON_BIAS0X'
    off_cdp = RON_CDP()
    off_cdp.rootname = 'OFFSET_BIAS0X'

    MB_profiles_cdp = cdp.CDP()
    MB_profiles_cdp.rootname = 'MB_profiles_BIAS0X'

    CDP_lib = dict(RON=ron_cdp,
                   OFF=off_cdp,
                   MB_profiles=MB_profiles_cdp)
    return CDP_lib
