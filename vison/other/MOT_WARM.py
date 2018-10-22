#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: MOT_WARM

Readiness verification: Warm Test before Cooling Down. 


Created on Mon Oct 22 17:11:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import unittest
from matplotlib.colors import Normalize


from vison.support import context
from vison.datamodel import scriptic as sc
from vison.datamodel import ccd
#import B01aux
from DarkTask import DarkTask
from vison.datamodel import inputs, cdp
from vison.support import utils
from vison.support.files import cPickleRead

# END IMPORT

isthere = os.path.exists


HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


MOT_WARM_commvalues = dict(program='CALCAMP', test='MOT_WARM', 
                         flushes=7, siflsh=1, siflsh_p=500,
                         inisweep=1,
                         vstart=0, vend=2086,
                         toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                         chinj=0,
                         s_tpump=0,
                         v_tpump=0,
                         exptime=0., 
                         shuttr=0, 
                         e_shuttr=0,
                         mirr_on=0,
                         wave=4,
                         motr_on=0,
                         source='flat',
                         comments='BIAS')


class MOT_WARM_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    #manifesto.update(OrderedDict(sorted([
    #    ('N', ([int], 'Number of Frame Acquisitions.')),
    #])))


class MOT_WARM(DarkTask):
    """ """

    inputsclass = MOT_WARM_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data), 
                         ('basic', self.basic_analysis)]

        super(MOT_WARM, self).__init__(inputs, log, drill, debug)
        self.name = 'MOT_WARM'
        self.type = 'Simple'
        
        
        self.HKKeys = HKKeys        
        #self.figdict = B01aux.B01figs.copy()
        #self.CDP_lib = B01aux.CDP_lib.copy()
        self.inputs['subpaths'] = dict(figs='figs',
                                       profiles='profiles', 
                                       products='products')
        

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds MOT_WARM script structure dictionary.
        
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """

        UP TO HERE
        
        MW_sdict = dict(col001=dict(frames=N, exptime=0))

        Ncols = len(BIAS01_sdict.keys())
        BIAS01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(BIAS01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        BIAS01_sdict = sc.update_structdict(
            BIAS01_sdict, commvalues, diffvalues)

        return BIAS01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(MOT_WARM, self).filterexposures(structure, explog, OBSID_lims, 
                                    colorblind=False,
                                    wavedkeys=wavedkeys)


    def basic_analysis(self):
        """ 

        MOT_WARM: Basic analysis of data.


        """



class Test(unittest.TestCase):
    """
    Unit tests for the BIAS01 class.
    """

    def setUp(self):

        inputs = dict()
        self.mw = MOT_WARM(inputs, log=None, drill=True, debug=False)

    def test_check_data(self):
        """

        :return: None
        """
        self.mw.check_data()


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(Test)
    unittest.TextTestRunner(verbosity=3).run(suite)
