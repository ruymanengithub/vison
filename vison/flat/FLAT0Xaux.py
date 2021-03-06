#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to FLAT0X.

Created on Tue Jan 30 17:39:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
from matplotlib import cm

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends


# END IMPORT

def get_CDP_lib():
    prnu_tb_cdp = cdp.Tables_CDP()
    prnu_tb_cdp.rootname = '%s_%snm_PRNU_TB'

    CDP_lib = dict(PRNU_TB=prnu_tb_cdp)
    return CDP_lib


def gt_check_offsets_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='%s_offset_vs_time.png' % (test,),
                caption='%s: offset vs. time.' % (ntest,),
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest),
                ylim=trends.offset_lims)


def gt_check_deltaoff_dict(test):
    ntest = test.replace('_', '\_')
    return dict(
        stats=[
            'deltaoff_pre', 'deltaoff_ove'], trendaxis='time', figname='%s_deltaoff_vs_time.png' %
        (test,), caption='%s: $\delta$offset vs. time. Offset value in each frame minus the average value.' %
        (ntest,), meta=dict(
            doLegend=True, doNiceXDate=True, suptitle='%s-checks: delta-offsets' %
            ntest), ylim=[
                -10., 10.])


def gt_check_std_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='%s_std_vs_time.png' % test,
                caption='%s: std vs. time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: std' % ntest,
                          ylim=trends.RON_lims))


def gt_check_img_flu_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['flu_med_img'],
                figname='%s_flu_vs_time.png' % test,
                caption='%s: Fluence vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: Fluence' % ntest,
                          ylabel='[ADU]'))


def gt_check_img_std_dict(test):
    ntest = test.replace('_', '\_')
    return dict(stats=['flu_std_img'],
                figname='%s_imgstd_vs_time.png' % test,
                caption='%s: Image area STD vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: Image-Area STD' % ntest,
                          ylabel='[ADU]'))


def gt_indiv_prof1Dhor_dict(test):
    ntest = test.replace('_', '\_')
    return dict(
        #figname='%s_profs1D_hor_allOBSIDs.png' % test,
        caption='%s: Average profiles across columns, PLACEHOLDER.' % ntest,
        meta=dict(doLegend=False,
                  ylabel='ADU',
                  xlabel='Column [pix]',
                  suptitle='%s: Profiles across columns, PLACEHOLDER.' % ntest)
    )


def gt_indiv_prof1Dver_dict(test):
    ntest = test.replace('_', '\_')
    return dict(
        #figname='%s_profs1D_ver_allOBSIDs.png' % test,
        caption='%s: Average profiles across rows, PLACEHOLDER.' % ntest,
        meta=dict(doLegend=False,
                  ylabel='ADU',
                  xlabel='Row [pix]',
                  suptitle='%s: Profiles across rows, PLACEHOLDER.' % ntest)
    )


def gt_meta_MFF2D_dict(test):
    ntest = test.replace('_', '\_')
    return dict(
        figname='%s_MASTERFLATFIELD_2Dimgshow_PLACEHOLDER.png' % test,
        caption='%s: Master FlatField for the CCDs [PLACEHOLDER]. Smoothed with gaussian kernel and displayed using histogram equalization to highlight structure.' % ntest,
        meta=dict(doLegend=False,
                  doColorbar=False,
                  suptitle='%s/Master:Quadrant Images [PLACEHOLDER].' % ntest,
                  corekwargs=dict(cmap=cm.gray, aspect='auto',  # norm=None,
                                  origin='lower left'))
    )


def gt_FL0Xfigs(test):
    FL0Xfigs = dict()
    FL0Xfigs['FL0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, gt_check_offsets_dict(test)]
    FL0Xfigs['FL0Xchecks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, gt_check_deltaoff_dict(test)]
    FL0Xfigs['FL0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, gt_check_std_dict(test)]
    FL0Xfigs['FL0Xchecks_flu'] = [
        trends.Fig_Basic_Checkstat, gt_check_img_flu_dict(test)]
    FL0Xfigs['FL0Xchecks_imgstd'] = [
        trends.Fig_Basic_Checkstat, gt_check_img_std_dict(test)]

    FL0Xfigs['FL0Xindiv_prof1D_hor_generic'] = [
        figclasses.Fig_Beam2DPlot, gt_indiv_prof1Dhor_dict(test)]
    FL0Xfigs['FL0Xindiv_prof1D_ver_generic'] = [
        figclasses.Fig_Beam2DPlot, gt_indiv_prof1Dver_dict(test)]
    FL0Xfigs['FL0Xmeta_MFF_2D'] = [
        figclasses.Fig_BeamImgShow, gt_meta_MFF2D_dict(test)]

    FL0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return FL0Xfigs
