#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ00 - used to test charge injection functionality using ELVIS.
  NOT intended for performance evaluation.


Created on Tue Jan 13 12:08:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from copy import deepcopy
from collections import OrderedDict

from vison.support import context
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#import datetime
from InjTask import InjTask
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

CHINJ00_commvalues = dict(program='CALCAMP', test='CHINJ00',
                          IG2_T=5.5, IG2_B=5.5,
                          flushes=7, siflsh=1,siflsh_p=500,
                          inisweep=1,
                          vstart=0, vend=2086,
                          chinj=1, chinj_on=30, chinj_of=100,
                          id_wid=60,
                          exptime=0., shuttr=0, e_shuttr=0,
                          mirr_on=0,
                          wave=4,
                          motr_on=0,
                          source='flat',
                          comments='')


class CHINJ00_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('IDH', ([float], 'Injection Drain High Voltage.')),
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
        ('chinj_on', ([int], 'Number of lines injected per cycle.')),
        ('chinj_of', ([int], 'Number of lines NON injected per cycle.'))
    ])))


class CHINJ00(InjTask):
    """ """

    inputsclass = CHINJ00_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data)]
        super(CHINJ00, self).__init__(inputs, log, drill, debug)
        self.name = 'CHINJ00'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = dict()
        self.inputs['subpaths'] = dict(figs='figs')
 
    def set_inpdefaults(self, **kwargs):
        """ """
        toi_chinj = 500

        self.inpdefaults = dict(
            IDH=18.,
            toi_chinj=toi_chinj,
            chinj_on=30,
            chinj_of=300)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """
        Builds CHINJ00 script structure dictionary.

        :param diffvalues: dict, opt, differential values.

        """

        id_wid = 60.  # us

        IDH = self.inputs['IDH']
        toi_chinj = self.inputs['toi_chinj']
        chinj_on = self.inputs['chinj_on']
        chinj_off = self.inputs['chinj_of']

        id_delays = np.array([2.5, 1.5]) * toi_chinj

        IG1s = [4., 6., 4., 6.]
        IG2s = [6., 4., 6., 4.]
        IDLs = [13, 13, 16., 16.]

        CCDs = [1, 2, 3]
        halves = ['T', 'B']

        assert len(id_delays) == 2

        CHINJ00_sdict = dict()

        colcounter = 1

        for j, id_delay in enumerate(id_delays):

            for i, IG1 in enumerate(IG1s):

                IDL = IDLs[i]
                IG2 = IG2s[i]

                colkey = 'col%03i' % colcounter
                #print colkey

                CHINJ00_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                             id_dly=id_delays[j],
                                             id_wid=id_wid,
                                             toi_ch=toi_chinj,
                                             chinj_on=chinj_on, chinj_off=chinj_off,
                                             comments='del%i_IG1%.1f' % (id_delay, IG1))

                for CCD in CCDs:
                    for half in halves:
                        CHINJ00_sdict[colkey]['IG1_%i_%s' % (CCD, half)] = IG1
                for half in halves:
                    CHINJ00_sdict[colkey]['IG2_%s' % (half,)] = IG2

                colcounter += 1

        Ncols = len(CHINJ00_sdict.keys())
        CHINJ00_sdict['Ncols'] = Ncols

        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(CHINJ00_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        CHINJ00_sdict = sc.update_structdict(
            CHINJ00_sdict, commvalues, diffvalues)

        return CHINJ00_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(CHINJ00, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                    wavedkeys=wavedkeys)
