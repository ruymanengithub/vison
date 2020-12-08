#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to NL01.

Created on Tue Jan 30 18:43:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT


def get_CDP_lib(testname):
    nl_tb_cdp = cdp.Tables_CDP()
    nl_tb_cdp.rootname = '{}_RESULTS_TABLE'.format(testname)

    curves_cdp = cdp.CDP()
    curves_cdp.rootname = '{}_CURVES'.format(testname)

    CDP_lib = dict(NL_TB=nl_tb_cdp,
                   NL_CURVES=curves_cdp)
    return CDP_lib


def get_check_offsets_dict(testname):
  return dict(stats=['offset_pre', 'offset_ove'],
      trendaxis='time',
      figname='{}_offset_vs_time.png'.format(testname),
      caption='{}: offset vs. time.'.format(testname),
      meta=dict(doLegend=True,
        doNiceXDate=True,
        suptitle='{}-checks: offsets'.format(testname),
        ylim=trends.offset_lims))

def get_check_deltaoff_dict(testname):
  return dict(
    stats=[
        'deltaoff_pre',
        'deltaoff_ove'],
    trendaxis='time',
    figname='{}_deltaoff_vs_time.png'.format(testname),
    caption='{}: $\delta$offset vs. time. '.format(testname)+\
      'Offset value in each frame minus the average value.',
    meta=dict(
      doLegend=True,
      doNiceXDate=True,
      suptitle='{}-checks: delta-offsets'.format(testname),
      ylim=[
          -10.,
          10.]))

def get_check_std_dict(testname):
  return dict(stats=['std_pre', 'std_ove'],
    trendaxis='time',
    figname='{}_std_vs_time.png'.format(testname),
    caption='{}: std vs. time.'.format(testname),
    meta=dict(doLegend=True,
      doNiceXDate=True,
      suptitle='{}-checks: std'.format(testname),
      ylim=trends.RON_lims))


def get_check_img_flu_dict(testname):
  return dict(stats=['flu_med_img'],
    trendaxis='exptime',
    figname='{}_flu_vs_exptime.png'.format(testname),
    caption='{}: Fluence vs. exposure time.'.format(testname),
    meta=dict(doLegend=True,
      doNiceXDate=False,
      suptitle='{}-checks: Fluence'.format(testname),
      ylabel='[ADU]',
      xlabel='exptime [s]'))

def get_check_img_flu_dict(testname):
  return dict(stats=['flu_std_img'],
    trendaxis='exptime',
    figname='{}_imgstd_vs_exptime.png'.format(testname),
    caption='{}: Image area STD vs. exposure time.'.format(testname),
    meta=dict(doLegend=False,
      doNiceXDate=False,
      suptitle='{}-checks: Image STD'.format(testname),
      ylabel='[ADU]',
      xlabel='exptime [s]'))

def get_NL_curves_dict(testname):
  return dict(
      figname='{}_curves.png'.format(testname),
      caption='{}: Relative non-linearity vs. '.format(testname)+\
        'fluence and best fit curves.',
      meta=dict(doLegend=True,
        ylabel='Z [percentage]',
        xlabel=r'$Y_{NL} [kADU]$',
        suptitle='{}: Non-Linearity Curves.'.format(testname),
        ylim=[-10., 10.],
        corekwargs=dict(data=dict(marker='.', linestyle='', color='b'),
                        fit=dict(marker='', linestyle='--', color='r')))
      )


def get_NL_singcurves_dict(testname):
  return dict(
    figname='{}_curves_single.png'.format(testname),
    caption='{}: Relative non-linearity vs. fluence and best fit curves.'.format(testname),
    meta=dict(doLegend=True,
      ylabel='Z [percentage]',
      xlabel=r'$Y_{NL} [kADU]$',
      title='{}: Non-Linearity Curves.'.format(testname),
      ylim=[-10., 10.],
      corekwargs = dict(linestyle='')
      #corekwargs=dict(data=dict(marker='.', linestyle='', color='b'),
      #                fit=dict(marker='', linestyle='--', color='r')))
      )
    )


def get_NL0Xfigs(testname):
    NL0Xfigs = dict()
    NL0Xfigs['NL0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, get_check_offsets_dict(testname)]
    NL0Xfigs['NL0Xchecks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, get_check_deltaoff_dict(testname)]
    NL0Xfigs['NL0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, get_check_std_dict(testname)]
    NL0Xfigs['NL0Xchecks_flu'] = [
        trends.Fig_Basic_Checkstat, get_check_img_flu_dict(testname)]
    NL0Xfigs['NL0Xchecks_imgstd'] = [
        trends.Fig_Basic_Checkstat, get_check_img_std_dict(testname)]
    NL0Xfigs['NL0X_fit_curves'] = [
        figclasses.Fig_Beam2DPlot, get_NL_curves_dict(testname)]
    NL0Xfigs['NL0X_fit_curves_single'] = [
        figclasses.Fig_XYPlot, get_NL_singcurves(testname)]
    NL0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return NL0Xfigs
