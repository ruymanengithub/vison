#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to PSF0X.

Created on Tue Jan 30 16:31:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st

from vison.plot import baseclasses as plbaseclasses
from vison.plot import trends

# END IMPORT

def gt_check_offsets_dict(test):
    ntest = st.replace(test,'_','\_')
    return dict(stats=['offset_pre','offset_ove'],
                          trendaxis='time',
                          figname='%s_offset_vs_time.png' % (test,),
                          caption='%s: offset vs. time.' % (ntest,),
                          meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest))

def gt_check_std_dict(test):
    ntest = st.replace(test,'_','\_')
    return dict(stats=['std_pre','std_ove'],
                      trendaxis='time',
                          figname='%s_std_vs_time.png' % test,
                          caption='%s: std vs. time.' % ntest,
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='%s-checks: std' % ntest))

def gt_check_bgd_dict(test):
    ntest = st.replace(test,'_','\_')
    return dict(stats=['bgd_img'],
                      trendaxis='time',
                          figname= '%s_bgd_vs_time.png' % test,
                          caption='%s: Background vs. time.' % ntest,
                          meta=dict(doLegend=False,
                                doNiceXDate=True,
                                suptitle='%s-checks: BGD' % ntest))

def gt_check_flu_dict(test):
    ntest = st.replace(test,'_','\_')
    return dict(stats=['chk_fluence'],
                      trendaxis='exptime',
                          figname='%s_flu_vs_mirr.png' % test,
                          caption='%s: Fluence vs. Mirror Position.' % ntest,
                          meta=dict(doLegend=False,
                                doNiceXDate=False,
                                suptitle='%s-checks: Fluence' % ntest,
                                xlabel='Mirr [mm]',
                                ylabel='Flu.[ADU]'))

def gt_check_fwhm_dict(test):
    ntest = st.replace(test,'_','\_')
    return dict(stats=['chk_fwhmx','chk_fwhmy'],
                      trendaxis='time',
                          figname='%s_fwhmxy_vs_mirr.png' % test,
                          caption='%s: FWHM(x,y) vs. Mirror Position.' % ntest,
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='%s-checks: FWHM(x,y)' % ntest,
                                ylabel='FWHM [pix]'))


def gt_PSF0Xfigs(test):
    PSF0Xfigs = dict()
    PSF0Xfigs['PSF0Xchecks_offsets'] = [trends.pl_basic_checkstat,gt_check_offsets_dict(test)]
    PSF0Xfigs['PSF0Xchecks_stds'] = [trends.pl_basic_checkstat,gt_check_std_dict(test)]
    PSF0Xfigs['PSF0Xchecks_bgd'] = [trends.pl_basic_checkstat,gt_check_bgd_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fluence'] = [trends.pl_basic_checkstat,gt_check_flu_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhm'] = [trends.pl_basic_checkstat,gt_check_fwhm_dict(test)]
    PSF0Xfigs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]
    return PSF0Xfigs

