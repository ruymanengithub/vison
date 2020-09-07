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
import copy
from matplotlib import cm

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends
from vison.point import Paux
# END IMPORT


def get_check_offsets_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='%s_offset_vs_time.png' % (test,),
                caption='%s: offset vs. time.' % (ntest,),
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest,
                          ylim=trends.offset_lims))


def get_check_deltaoff_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(
        stats=[
            'deltaoff_pre', 'deltaoff_ove'], trendaxis='time', figname='%s_deltaoff_vs_time.png' %
        (test,), caption='%s: $\delta$offset vs. time. Offset value in each frame minus the average value.' %
        (ntest,), meta=dict(
            doLegend=True, doNiceXDate=True, suptitle='%s-checks: delta-offsets' %
            ntest, ylim=[
                -10., 10.]))


def get_check_std_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='%s_std_vs_time.png' % test,
                caption='%s: std vs. time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: std' % ntest,
                          ylim=trends.RON_lims))


def get_check_bgd_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(stats=['bgd_img'],
                trendaxis='time',
                figname='%s_bgd_vs_time.png' % test,
                caption='%s: Background vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: BGD' % ntest))


def get_check_flu_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(stats=['chk_fluence'],
                trendaxis='exptime',
                figname='%s_flu_vs_exptime.png' % test,
                caption='%s: Fluence vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: Fluence' % ntest,
                          xlabel='seconds',
                          ylabel='Flu.[ADU]'))


def get_check_fwhmx_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(stats=['chk_fwhmx'],
                trendaxis='exptime',
                figname='%s_fwhmx_vs_exptime.png' % test,
                caption='%s: FWHM(x) vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: FWHM(x)' % ntest,
                          xlabel='seconds',
                          ylabel='FWHMx [pix]'))


def get_check_fwhmy_dict(test):
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    return dict(stats=['chk_fwhmy'],
                trendaxis='exptime',
                figname='%s_fwhmy_vs_exptime.png' % test,
                caption='%s: FWHM(y) vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: FWHM(y)' % ntest,
                          xlabel='seconds',
                          ylabel='FWHMy [pix]'))


def get_crosstalk_dict(test, figtype):
    tcaption = '%s: Cross-Talk [%s]. Green means positive cross-talk, red means negative cross-talk' +\
        ' (does not mean compliance/non-compliance). Pale colours mean less accurate results.'
    ntest = copy.deepcopy(test)
    ntest.replace('_', '\_')
    crosstalk_dict = dict(
        figname='%s_crosstalk_%s.png' % (test, figtype),
        caption=tcaption % (ntest, figtype),
        meta=dict(),
        data=None
    )
    return crosstalk_dict

def get_spotsposter_dict(test, BFE=True):
    """ """
    if BFE:
        figtype ='noBFE'
        tcaption = '%s: Spots Poster, BFE not corrected. Log scale.'
    else:
        figtype = 'withBFE'
        tcaption = '%s: Spots Poster, BFE corrected using G+15. Log scale.'

    ntest = copy.deepcopy(test)
    stop()
    ntest.replace('_', '\_')
    sp_dict = dict(
        figname='%s_spotsposter_%s.png' % (test, figtype),
        caption=tcaption % (ntest,),
        meta=dict(doColorbar=True,
            corekwargs=dict(cmap=cm.gray, 
                            aspect='auto',  
                            # norm=None,
                            origin='lower left')),
        data=None
        )
    return sp_dict

def get_PSF0Xfigs(test):
    PSF0Xfigs = dict()
    PSF0Xfigs['PSF0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, get_check_offsets_dict(test)]
    PSF0Xfigs['PSF0Xchecks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, get_check_deltaoff_dict(test)]
    PSF0Xfigs['PSF0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, get_check_std_dict(test)]
    PSF0Xfigs['PSF0Xchecks_bgd'] = [
        trends.Fig_Basic_Checkstat, get_check_bgd_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fluence'] = [
        trends.Fig_Basic_Checkstat, get_check_flu_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhmx'] = [
        trends.Fig_Basic_Checkstat, get_check_fwhmx_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhmy'] = [
        trends.Fig_Basic_Checkstat, get_check_fwhmy_dict(test)]
    PSF0Xfigs['PSF0X_crosstalk_ADU'] = [
        figclasses.Fig_Husk, get_crosstalk_dict(test, 'ADU')]
    PSF0Xfigs['PSF0X_crosstalk_RATIO'] = [
        figclasses.Fig_Husk, get_crosstalk_dict(test, 'RATIO')]
    PSF0Xfigs['SpotsPoster'] = [
        figclasses.Fig_ImgShow, get_spotsposter_dict(test, BFE=False)]
    PSF0Xfigs['SpotsPosterNOBFE'] = [
        figclasses.Fig_ImgShow, get_spotsposter_dict(test, BFE=True)]
    PSF0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return PSF0Xfigs


def get_PSF01_PANCHRO_figs():
    PSF01_PANCHRO_figs = dict()
    return PSF01_PANCHRO_figs


def get_CDP_lib(test):
    """ """
    CDP_lib = OrderedDict()
    CDP_lib['RAW_CTALK'] = cdp.CDP()
    CDP_lib['RAW_CTALK'].rootname = 'Raw_crosstalk'

    CDP_lib['CTALK'] = cdp.CDP()
    CDP_lib['CTALK'].rootname = 'crosstalk'

    CDP_lib.update(Paux.get_CDP_lib())

    return CDP_lib
