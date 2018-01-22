# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PSF01_PANCHRO

PSF vs. Fluence, Wavelength, nominal temperature
   Meta-test across Wavelengths

Tasks:

    - Gather results at several wavelengths
    - Do synthetic analysis of results.


Created on Thu Nov 30 16:38:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import string as st
import warnings
import copy
from collections import OrderedDict

from vison.datamodel import core
from vison.datamodel import ccd
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.ogse import ogse
#from vison.datamodel import HKtools
from vison.datamodel import scriptic as sc
from vison.flat import FlatFielding as FFing
#from vison.point import lib as polib
#from vison.support.report import Report
from vison.support import files
from vison.image import calibration
from vison.pipe.task import Task
import PSF0Xaux
from vison.image import performance
from vison.support.files import cPickleRead
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists


class PSF01_PANCHRO_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('inpath',([str],'Data path.')),
            ('wavelengths',([list],'Wavelengths.')),
            ])))


class PSF01_PANCHRO(Task):
    
    inputsclass = PSF01_PANCHRO_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(PSF01_PANCHRO,self).__init__(inputs,log,drill,debug)
        self.name = 'PSF01_PANCHRO'
        self.type = 'Meta'
        self.subtasks = [('check',self.check_data),
                         ('meta',self.meta)]
        #self.HKKeys = HKKeys
        self.figdict = PSF0Xaux.PSF0Xfigs
        self.perflimits.update(performance.perf_rdout)

    def ingest_data_MetaTest(self):
        """ """
        
        wavelengths = self.inputs['wavelengths']
        inpath = self.inputs['inpath']
        
        self.dd = dict()        
        for iw,wave in enumerate(wavelengths):
            wpath = 'PSF01_%i' % wave
            iDDpickfile = os.path.join(inpath,wpath,'PSF01_%i_DataDict.pick' % wave)
            self.dd['nm%i' % wave] = cPickleRead(iDDpickfile)
        
    
    def check_data(self):
        """ """
    
    def meta(self):
        """ """
    
