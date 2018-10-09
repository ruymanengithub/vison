#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BF01.

Created on Tue Jul 31 17:50:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import pandas as pd

from vison.flat import PTC0Xaux
from vison.plot import figclasses
from vison.plot import trends

from vison.datamodel import cdp

# END IMPORT



covtable_cdp = cdp.Tables_CDP()
covtable_cdp.rootname = 'BF01_COVTABLE'

bftable_cdp = cdp.Tables_CDP()
bftable_cdp.rootname = 'BF01_G15TABLE'

CDP_lib = dict(COVTABLE=covtable_cdp,
               BFTABLE=bftable_cdp)


prof_COV_ver_dict = dict(
    figname='BF01_COV_profs_ver.png',
    caption='BF01: COV 1D profiles, vertical/parallel direction.',
    meta=dict(doLegend=True,
              ylabel='COV/VAR, ADIM.',
              xlabel='Y',
              ylim=[-0.005,0.03],
              suptitle='BF01: COV 1D Profile, Vertical/Parallel.')
)

prof_COV_ser_dict = dict(
    figname='BF01_COV_profs_ser.png',
    caption='BF01: COV 1D profiles, serial direction.',
    meta=dict(doLegend=True,
              ylabel='COV/VAR, ADIM.',
              xlabel='X',
              ylim=[-0.005,0.03],
              suptitle='BF01: COV 1D Profile, Serial.')
)



def gt_BF01figs(test):
    
    BF01figs = dict()
    BF01figs['BF01checks_offsets'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_offsets_dict(test)]
    BF01figs['BF01checks_stds'] = [
        trends.Fig_Basic_Checkstat, PTC0Xaux.gt_check_std_dict(test)]
    BF01figs['BF01checks_flu'] = [
        trends.Fig_Basic_Checkstat,  PTC0Xaux.gt_check_img_flu_dict(test)]
    BF01figs['BF01checks_imgstd'] = [
        trends.Fig_Basic_Checkstat,  PTC0Xaux.gt_check_img_std_dict(test)]
    BF01figs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    
    # renaming to BF01... HACK
    
    keys_to_rename = ['caption', 'suptitle', 'figname']
    for figkey in BF01figs.keys():
        _dict = BF01figs[figkey][1]
        try:
            _mdict = BF01figs[figkey][1]['meta']
            hasmeta=True
        except KeyError:
            hasmeta=False
        for key in keys_to_rename:
            if key in _dict:
                _dict[key] = st.replace(_dict[key],test,'BF01')
            if hasmeta:
                if key in _mdict:
                    _mdict[key] = st.replace(_mdict[key],test,'BF01')
    
    
    BF01figs['BF01_COV_ver'] = [
    figclasses.Fig_Beam2DPlot, prof_COV_ver_dict]
    
    BF01figs['BF01_COV_hor'] = [
    figclasses.Fig_Beam2DPlot, prof_COV_ser_dict]
    
    return BF01figs
