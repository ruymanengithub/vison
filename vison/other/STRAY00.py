#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: STRAY00 - used to investigate STRAY-LIGHT sources in OGSE.
  NOT intended for performance evaluation.
  COMMISSIONING.


Created on Thu Feb 08 14:07:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from copy import deepcopy
from collections import OrderedDict

from vison.pipe.task import HKKeys
from vison.support import context
#from vison.pipe import lib as pilib
#from vison.point import lib as polib
#from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#import datetime
#from InjTask import InjTask
from vison.image import performance
from vison.datamodel import inputs
from vison.dark.DarkTask import DarkTask
# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

STRAY00_commvalues = dict(program='CALCAMP', test='STRAY00',
                          flushes=7, siflsh=1, siflsh_p=500,
                          inisweep=1,
                          vstart=0, vend=2086,
                          toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                          chinj=0,
                          s_tpump=0,
                          v_tpump=0,
                          exptime=0., shuttr=0, e_shuttr=0,
                          mirr_on=0,
                          wave=4,
                          motr_on=0,
                          source='flat',
                          comments='')


class STRAY00_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    # manifesto.update(OrderedDict(sorted([
    #        ('IDH',([float],'Injection Drain High Voltage.')),
    #        ('toi_chinj',([int],'TOI Charge Injection.')),
    #        ('chinj_on',([int],'Number of lines injected per cycle.')),
    #        ('chinj_of',([int],'Number of lines NON injected per cycle.'))
    #        ])))


class STRAY00(DarkTask):
    """ """

    inputsclass = STRAY00_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data)]
        super(STRAY00, self).__init__(inputs=inputs, log=log, drill=drill,
                                      debug=debug, cleanafter=cleanafter)
        self.name = 'STRAY00'
        self.type = 'Simple'

        self.HKKeys = HKKeys
        self.figdict = dict()
        self.inputs['subpaths'] = dict(figs='figs')

    def set_inpdefaults(self, **kwargs):
        """ """
        self.inpdefaults = dict()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """
        Builds STRAY00 script structure dictionary.
        :param diffvalues: dict, opt, differential values.

        """

        STRAY00_sdict = dict()

        # start with all lights on in LAB
        STRAY00_sdict['col001'] = dict(frames=1, exptime=0, shuttr=0, wave=1,
                                       source='flat',
                                       comments='LABLIT')
        STRAY00_sdict['col002'] = dict(frames=1, exptime=100, shuttr=0, wave=1,
                                       source='flat',
                                       comments='LABLIT')
        # switch off all lights in LAB
        STRAY00_sdict['col003'] = dict(frames=1, exptime=0, shuttr=0, wave=1,
                                       source='flat',
                                       comments='LABOUTW1')
        STRAY00_sdict['col004'] = dict(frames=1, exptime=100, shuttr=0, wave=1,
                                       source='flat',
                                       comments='LABOUTW1')
        # change light source wavelength
        STRAY00_sdict['col005'] = dict(frames=1, exptime=0, shuttr=0, wave=4,
                                       source='flat',
                                       comments='LABOUTW4')
        STRAY00_sdict['col006'] = dict(frames=1, exptime=100, shuttr=0, wave=4,
                                       source='flat',
                                       comments='LABOUTW4')
        # switch off Light Source
        STRAY00_sdict['col007'] = dict(frames=1, exptime=0, shuttr=0, wave=1,
                                       comments='TUNGSOFF')
        STRAY00_sdict['col008'] = dict(frames=1, exptime=100, shuttr=0, wave=1,
                                       source='flat',
                                       comments='TUNGSOFF')
        # switch off Pressure Gauge - Light Source still Off
        STRAY00_sdict['col009'] = dict(frames=1, exptime=0, shuttr=0, wave=1,
                                       comments='GAUGEOFF')
        STRAY00_sdict['col010'] = dict(frames=1, exptime=100, shuttr=0, wave=1,
                                       source='flat',
                                       comments='GAUGEOFF')
        # Switch on Pressure Gauge and Light Source
        STRAY00_sdict['col011'] = dict(frames=1, exptime=0, shuttr=0, wave=1,
                                       source='flat',
                                       comments='TUNGGAUON')
        STRAY00_sdict['col012'] = dict(frames=1, exptime=100, shuttr=0, wave=1,
                                       source='flat',
                                       comments='TUNGGAUON')

        # NON-supervised, Lights Off, Gauge On, Tungsten On

        # DARK FF-source wave = 1, 4

        STRAY00_sdict['col013'] = dict(frames=1, exptime=300, shuttr=0, wave=1,
                                       source='flat',
                                       comments='DARK-FF')

        STRAY00_sdict['col014'] = dict(frames=1, exptime=300, shuttr=0, wave=4,
                                       source='flat',
                                       comments='DARK-FF')

        # DARK PSF-source wave = 1, mirror at near end

        STRAY00_sdict['col015'] = dict(frames=1, exptime=300, shuttr=0, wave=1,
                                       source='point',
                                       mirr_on=1, mirr_pos=1.,
                                       comments='DARK-PSF')

        # DARK PSF-source wave = 1, mirror at far end

        STRAY00_sdict['col016'] = dict(frames=1, exptime=300, shuttr=0, wave=1,
                                       source='point',
                                       mirr_on=1, mirr_pos=99.,
                                       comments='DARK-PSF')

        # DARK PSF-source wave = 4, mirror at near end

        STRAY00_sdict['col017'] = dict(frames=1, exptime=300, shuttr=0, wave=4,
                                       source='point',
                                       mirr_on=1, mirr_pos=1.,
                                       comments='DARK-PSF')

        # DARK PSF-source wave = 4, mirror at far end

        STRAY00_sdict['col018'] = dict(frames=1, exptime=300, shuttr=0, wave=4,
                                       source='point',
                                       mirr_on=1, mirr_pos=99.,
                                       comments='DARK-PSF')

        # PSF - Short, wave = 1

        STRAY00_sdict['col019'] = dict(frames=1, exptime=5., shuttr=1, wave=1,
                                       source='point',
                                       mirr_on=1, mirr_pos=50.,
                                       comments='PSF-S')

        STRAY00_sdict['col020'] = dict(frames=1, exptime=50., shuttr=1, wave=1,
                                       source='point',
                                       mirr_on=1, mirr_pos=50.,
                                       comments='PSF-L')

        # PSF - Long, wave = 4

        STRAY00_sdict['col021'] = dict(frames=1, exptime=5., shuttr=1, wave=4,
                                       source='point',
                                       mirr_on=1, mirr_pos=50.,
                                       comments='PSF-S')

        STRAY00_sdict['col022'] = dict(frames=1, exptime=50., shuttr=1, wave=4,
                                       source='point',
                                       mirr_on=1, mirr_pos=50.,
                                       comments='PSF-L')

        Ncols = len(list(STRAY00_sdict.keys()))
        STRAY00_sdict['Ncols'] = Ncols

        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(STRAY00_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
                diffvalues = diffvalues = dict()

        STRAY00_sdict = sc.update_structdict(
            STRAY00_sdict, commvalues, diffvalues)

        return STRAY00_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(STRAY00, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                    wavedkeys=wavedkeys)
