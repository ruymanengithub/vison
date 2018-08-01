# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: MOT_FF

Brighter-Fatter Analysis
   Using data from test PTC01 (via BF01)
Hard Edge Response in serial / paralle
Bit Correlations (ADC health)

Created on Tue Jul 31 18:04:00 2018

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import warnings
import copy
import string as st
from collections import OrderedDict

from vison.pipe.task import Task
from vison.flat.BF01 import BF01
from vison.other import MOT_FFaux
from vison.datamodel import inputs
from vison.support.files import cPickleRead, cPickleDumpDictionary
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


class MOT_FF_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
        ('Npix', ([int], 'Number of Pixels (linear) to consider for Covariance Matrix')),
        ('surrogate', ([str], 'Test to use as surrogate'))
    ])))


class MOT_FF(BF01):
    """ """

    inputsclass = MOT_FF_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract_COV', self.extract_COV),
                         ('extract_BF', self.extract_BF),
                         ('extract_ADC', self.extract_ADC),
                         ('extract_HER', self.extract_HER),
                         ('meta', self.meta_analysis)]
        
        super(MOT_FF, self).__init__(inputs, log, drill, debug)
        self.name = 'MOT_FF'
        #self.type = 'Simple'
        
        #self.HKKeys = HKKeys
        self.figdict = MOT_FFaux.gt_MOT_FF_figs(self.inputs['surrogate'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles')

    def check_data(self):
        
        kwargs = dict(figkeys=['MOT_FFchecks_offsets', 'MOT_FFchecks_stds',
                                   'MOT_FFchecks_flu', 'MOT_FFchecks_imgstd'])
        
        Task.check_data(self, **kwargs)
    
    def extract_ADC(self):
        
        stop()
        
    
    def extract_HER(self):
        return