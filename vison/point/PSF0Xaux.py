#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to PSF0X.

Created on Tue Jan 30 16:31:00 2018

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
from vison.point import Paux
# END IMPORT


def gt_check_offsets_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='%s_offset_vs_time.png' % (test,),
                caption='%s: offset vs. time.' % (ntest,),
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest,
                          ylim=trends.offset_lims))


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


def gt_check_bgd_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['bgd_img'],
                trendaxis='time',
                figname='%s_bgd_vs_time.png' % test,
                caption='%s: Background vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: BGD' % ntest))


def gt_check_flu_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['chk_fluence'],
                trendaxis='exptime',
                figname='%s_flu_vs_exptime.png' % test,
                caption='%s: Fluence vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: Fluence' % ntest,
                          xlabel='seconds',
                          ylabel='Flu.[ADU]'))


def gt_check_fwhmx_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['chk_fwhmx'],
                trendaxis='exptime',
                figname='%s_fwhmx_vs_exptime.png' % test,
                caption='%s: FWHM(x) vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: FWHM(x)' % ntest,
                          xlabel='seconds',
                          ylabel='FWHMx [pix]'))


def gt_check_fwhmy_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['chk_fwhmy'],
                trendaxis='exptime',
                figname='%s_fwhmy_vs_exptime.png' % test,
                caption='%s: FWHM(y) vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: FWHM(y)' % ntest,
                          xlabel='seconds',
                          ylabel='FWHMy [pix]'))


def gt_PSF0Xfigs(test):
    PSF0Xfigs = dict()
    PSF0Xfigs['PSF0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, gt_check_offsets_dict(test)]
    PSF0Xfigs['PSF0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, gt_check_std_dict(test)]
    PSF0Xfigs['PSF0Xchecks_bgd'] = [
        trends.Fig_Basic_Checkstat, gt_check_bgd_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fluence'] = [
        trends.Fig_Basic_Checkstat, gt_check_flu_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhmx'] = [
        trends.Fig_Basic_Checkstat, gt_check_fwhmx_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhmy'] = [
        trends.Fig_Basic_Checkstat, gt_check_fwhmy_dict(test)]
    PSF0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return PSF0Xfigs


def get_PSF01_PANCHRO_figs():
    PSF01_PANCHRO_figs = dict()
    return PSF01_PANCHRO_figs

def gt_CDP_lib(test):
    """ """
    CDP_lib = OrderedDict()
    CDP_lib.update(Paux.CDP_lib)