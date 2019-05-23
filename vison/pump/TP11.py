#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP11

Trap-Pumping calibration (vertical)

Created on Thu May 23 10:42:00 2019

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os

from vison.pipe.task import HKKeys
from vison.pump import TP01
import TP01aux
# END IMPORT

isthere = os.path.exists


TP11_commvalues = TP01.TP01_commvalues.copy()
TP11_commvalues['test'] = 'TP11'
TP11_commvalues['IDL'] = 10.5


class TP11(TP01.TP01):
    """ """

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('injection', self.charact_injection),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('debugtask',self.debugtask),
                         ('meta', self.meta_analysis)]
        super(TP11, self).__init__(inputs=inputs, log=log, 
            drill=drill, debug=debug, cleanafter=cleanafter)
        self.name = 'TP11'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = TP01aux.get_TP01figs()
        self.CDP_lib = TP01aux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       products='products')
        self.commvalues = TP11_commvalues.copy()
        
        
