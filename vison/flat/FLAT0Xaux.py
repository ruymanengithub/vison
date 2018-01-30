#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to FLAT0X.

Created on Tue Jan 30 17:39:00 2018

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


def gt_check_img_fluvar_dict(test):
    ntest = st.replace(test,'_','\_')
    return dict(stats=['flu_med_img','flu_var_img'],
                          figname='%s_fluvar_vs_time.png' % test,
                          caption='%s: Fluence \& Variance vs. time.' % ntest,
                          meta=dict(doLegend=True,
                                doNiceXDate=True,
                                suptitle='%s-checks: Fluence \& Variance' % ntest,
                                ylabel='[ADU]'))


def gt_FL0Xfigs(test):
    FL0Xfigs = dict()
    FL0Xfigs['FL0Xchecks_offsets'] = [trends.pl_basic_checkstat,gt_check_offsets_dict(test)]
    FL0Xfigs['FL0Xchecks_stds'] = [trends.pl_basic_checkstat,gt_check_std_dict(test)]
    FL0Xfigs['FL0Xchecks_fluvar'] = [trends.pl_basic_checkstat,gt_check_img_fluvar_dict(test)]
    FL0Xfigs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]
    return FL0Xfigs

