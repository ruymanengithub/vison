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

from vison.plot import figclasses
from vison.plot import trends

# END IMPORT


def gt_check_offsets_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='%s_offset_vs_time.png' % (test,),
                caption='%s: offset vs. time.' % (ntest,),
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest),
                          ylim=trends.offset_lims)


def gt_check_std_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='%s_std_vs_time.png' % test,
                caption='%s: std vs. time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: std' % ntest,
                          ylim=trends.RON_lims))


def gt_check_img_flu_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['flu_med_img'],
                figname='%s_flu_vs_time.png' % test,
                caption='%s: Fluence vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: Fluence' % ntest,
                          ylabel='[ADU]'))


def gt_check_img_std_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['flu_std_img'],
                figname='%s_imgstd_vs_time.png' % test,
                caption='%s: Image area STD vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: Image-Area STD' % ntest,
                          ylabel='[ADU]'))
def gt_indiv_prof1Dhor_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(
        #figname='%s_profs1D_hor_allOBSIDs.png' % test,
        caption='%s: Average profiles across columns, PLACEHOLDER.' % ntest,
        meta=dict(doLegend=False,
                  ylabel='ADU',
                  xlabel='Column [pix]',
                  suptitle='%s: Profiles across columns, PLACEHOLDER.' % ntest)
        )

def gt_indiv_prof1Dver_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(
        #figname='%s_profs1D_ver_allOBSIDs.png' % test,
        caption='%s: Average profiles across rows, PLACEHOLDER.' % ntest,
        meta=dict(doLegend=False,
                  ylabel='ADU',
                  xlabel='Row [pix]',
                  suptitle='%s: Profiles across rows, PLACEHOLDER.' % ntest)
        )



def gt_FL0Xfigs(test):
    FL0Xfigs = dict()
    FL0Xfigs['FL0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, gt_check_offsets_dict(test)]
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
    
    FL0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return FL0Xfigs
