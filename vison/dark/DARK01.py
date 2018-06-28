#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: DARK01

"Dark Current" analysis script

Created on Tue Aug 29 17:21:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
from copy import deepcopy
from collections import OrderedDict
import copy

from vison.support import context
#from vison.pipe import lib as pilib
#from vison.point import lib as polib
from vison.datamodel import scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#from vison.pipe.task import Task
from DarkTask import DarkTask
import D01aux
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


DARK01_commvalues = dict(program='CALCAMP', test='DARK01',
                         flushes=7, shuttr=0,
                         siflsh=1, siflsh_p=500,
                         comments='DARK')


class DARK01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('N', ([int], 'Number of Frame Acquisitions.')),
        ('exptime', ([float], 'Exposure time.')),
    ])))


Flu_lims = OrderedDict(CCD1=OrderedDict(E=[-2., 5.]))
for Q in ['F', 'G', 'H']:
    Flu_lims['CCD1'][Q] = copy.deepcopy(Flu_lims['CCD1']['E'])
for i in [2, 3]:
    Flu_lims['CCD%i' % i] = copy.deepcopy(Flu_lims['CCD1'])


class DARK01(DarkTask):
    """ """

    inputsclass = DARK01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(DARK01, self).__init__(inputs, log, drill, debug)
        self.name = 'DARK01'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = D01aux.D01figs.copy()
        self.inputs['subpaths'] = dict(figs='figs')
        self.init_todo_flags()

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = dict(N=4, exptime=565)

    def set_perfdefaults(self, **kwargs):
        super(DARK01, self).set_perfdefaults(**kwargs)
        self.perfdefaults['Flu_lims'] = Flu_lims

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds DARK01 script structure dictionary.

        :param diffvalues: dict, opt, differential values.

        """

        N = self.inputs['N']
        exptime = self.inputs['exptime']
        DARK01_sdict = dict(col1=dict(frames=N, exptime=exptime))

        Ncols = len(DARK01_sdict.keys())
        DARK01_sdict['Ncols'] = Ncols

        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(DARK01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        DARK01_sdict = sc.update_structdict(
            DARK01_sdict, commvalues, diffvalues)

        return DARK01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(DARK01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                   wavedkeys=wavedkeys)

    def prep_data(self):
        """

        DARK01: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            [BIAS SUBTRACTION]
            cosmetics masking

        """
        super(DARK01, self).prepare_images(
            doExtract=True, doMask=True, doOffset=True, doBias=True)

    def basic_analysis(self):
        """ 

        DARK01: Basic analysis of data.

        **METACODE**

        ::

            f. e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        produce mask of hot pixels
                        count hot pixels / columns
                        produce a 2D poly model of masked-image, save coefficients
                        produce average profile along rows
                        produce average profile along cols
                        measure and save RON after subtracting large scale structure
                    save 2D model and profiles in a pick file for each OBSID-CCD

            plot average profiles f. each CCD and Q (color coded by time)

        """

        raise NotImplementedError

    def meta_analysis(self):
        """ 

        **METACODE**

        ::

            f. each CCD:
                f. e. Q:
                    stack all ObsIDs to produce Master Dark
                    produce mask of hot pixels / columns
                    count hot pixels / columns
                    measure average profile along rows
                    measure average profile along cols

            plot average profiles of Master Bias f. each CCD,Q
            show Master Dark (images), include in report
            report stats of defects, include in report
            save name of MasterDark to DataDict, report
            save name of Defects in Darkness Mask to DD, report


        """

        raise NotImplementedError
