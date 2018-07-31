#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BF01.

Created on Tue Jul 31 17:50:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st

from vison.flat import PTC0Xaux
from vison.plot import figclasses
from vison.plot import trends

# END IMPORT





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
    return BF01figs
