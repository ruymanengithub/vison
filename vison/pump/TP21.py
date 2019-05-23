#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP02

Trap-Pumping calibration (serial)

Created on Tue Aug 29 17:38:00 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os

from vison.pipe.task import HKKeys
#from vison.pipe.task import Task
from vison.pump import TP02
import TP02aux
# END IMPORT

isthere = os.path.exists


TP21_commvalues = TP02.TP02_commvalues.copy()
TP21_commvalues['test'] = 'TP21'
TP21_commvalues['IDL'] = 10.5





class TP21(TP02.TP02):
    """ """


    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(TP21, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug, 
                    cleanafter=cleanafter)
        self.name = 'TP21'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = TP02aux.get_TP02figs()
        self.CDP_lib = TP02aux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs', 
                   ccdpickles='ccdpickles',
                   products='products')
        
        self.commvalues = TP21_commvalues.copy()

