#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Common Use Task for Trap-Pumping Analysis.

Created on Tue Jan  2 17:44:04 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import copy
import os

from vison.pipe.task import Task
from vison.datamodel import core, ccd
from vison.pipe import lib as pilib
#from vison.pipe.task import Task
from vison.inject.InjTask import InjTask
# END IMPORT


class PumpTask(InjTask):

    def __init__(self, *args, **kwargs):
        super(PumpTask, self).__init__(*args, **kwargs)

    def check_data(self, **kwargs):
        """ """
        test = self.inputs['test']
        if test == 'TP01':
            kwargs = dict(pattern=(2066, 0, 1),
                          figkeys=['TP01checks_offsets', 'TP01checks_stds',
                                   'TP01checks_injlevel', 'TP01checks_injstd'])
        elif test == 'TP02':
            kwargs = dict(pattern=(2066, 0, 1),
                          figkeys=['TP02checks_offsets', 'TP02checks_stds',
                                   'TP02checks_injlevel', 'TP02checks_injstd'])

        InjTask.check_data(self, **kwargs)

    def check_metrics_ST(self, **kwargs):
        """ 

        """
        super(PumpTask, self).check_metrics_ST(**kwargs)
