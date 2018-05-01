#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to PTC0X.

Created on Tue Jan 30 17:51:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

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
                          suptitle='%s-checks: offsets' % ntest))


def gt_check_std_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='%s_std_vs_time.png' % test,
                caption='%s: std vs. time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: std' % ntest))


def gt_check_img_flu_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['flu_med_img'],
                trendaxis='exptime',
                figname='%s_flu_vs_exptime.png' % test,
                caption='%s: Fluence vs. exposure time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=False,
                          suptitle='%s-checks: Fluence' % ntest,
                          xlabel='exptime[s]',
                          ylabel='[ADU]'))


def gt_check_img_std_dict(test):
    ntest = st.replace(test, '_', '\_')
    return dict(stats=['flu_std_img'],
                trendaxis='exptime',
                figname='%s_imgstd_vs_exptime.png' % test,
                caption='%s: STD vs. exposure time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=False,
                          suptitle='%s-checks: Image Area STD' % ntest,
                          xlabel='exptime[s]',
                          ylabel='[ADU]'))


def gt_PTC0Xfigs(test):
    PTC0Xfigs = dict()
    PTC0Xfigs['PTC0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, gt_check_offsets_dict(test)]
    PTC0Xfigs['PTC0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, gt_check_std_dict(test)]
    PTC0Xfigs['PTC0Xchecks_flu'] = [
        trends.Fig_Basic_Checkstat, gt_check_img_flu_dict(test)]
    PTC0Xfigs['PTC0Xchecks_imgstd'] = [
        trends.Fig_Basic_Checkstat, gt_check_img_std_dict(test)]
    PTC0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return PTC0Xfigs
