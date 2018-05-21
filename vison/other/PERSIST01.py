#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PERSIST01

CCD Persistence test

Created on Tue Aug 29 17:39:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
#import numpy as np
from pdb import set_trace as stop
#import os
#import datetime
from copy import deepcopy
from collections import OrderedDict

from vison.support import context
#from vison.pipe import lib as pilib
#from vison.point import lib as polib
from vison.support import ogse
from vison.datamodel import scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
#from vison.datamodel import ccd
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
from vison.pipe.task import Task
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

HKKeys = []

PER01_commvalues = dict(program='CALCAMP', test='PERSIST01',
                        rdmode='fwd_bas',
                        flushes=7,  # exptime=0.,
                        vstart=0, vend=2086,
                        shuttr=1,
                        siflsh=1, siflsh_p=500,
                        chinj=0,
                        wave=4,
                        mirr_pos=ogse.mirror_nom['F4'],
                        mirr_on=1,
                        comments='')


class PERSIST01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptSATUR', ([float], 'Exposure times to produce latent.')),
        ('exptLATEN', ([float], 'Exposure times to quantify latent.')),
    ])))


class PERSIST01(Task):
    """ """

    inputsclass = PERSIST01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(PERSIST01, self).__init__(inputs, log, drill, debug)
        self.name = 'PERSIST01'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = dict()
        self.inputs['subpaths'] = dict(figs='figs')

    def set_inpdefaults(self, **kwargs):
        """ """
        self.inpdefaults = dict(
            exptSATUR=15.,
            exptLATEN=565.
        )

    def set_perfdefaults(self, **kwargs):
        """ """
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ 
        Builds PERSISTENCE01 script structure dictionary.

        :param exptSATUR: int, saturation exposure time.
        :param exptLATEN: int, latency exposure time.
        :param diffvalues: dict, opt, differential values.


        """
        exptSATUR = self.inputs['exptSATUR']
        exptLATEN = self.inputs['exptLATEN']

        PER01_sdict = dict(col1=dict(frames=5, exptime=0, shuttr=0, comments='REFER.'),
                           col2=dict(frames=1, exptime=exptSATUR,
                                     comments='EXPOSE'),
                           col3=dict(frames=3, exptime=exptLATEN, shuttr=0, comments='LATENT'))

        Ncols = len(PER01_sdict.keys())
        PER01_sdict['Ncols'] = Ncols

        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(PER01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        PER01_sdict = sc.update_structdict(PER01_sdict, commvalues, diffvalues)

        return PER01_sdict

    def filterexposures(self, structure, explogf, datapath, OBSID_lims):
        """ """
        wavedkeys = []
        return super(PERSIST01, self).filterexposures(structure, explogf, datapath, OBSID_lims, colorblind=True,
                                                      wavedkeys=wavedkeys)

    def check_data(self):
        """ 

        PERSIST01: Checks quality of ingested data.


        **METACODE**

        ::

            check common HK values are within safe / nominal margins
            check voltages in HK match commanded voltages, within margins

            f.e.ObsID:
                f.e.CCD:
                    f.e.Q.:
                        measure offsets in pre-, over-
                        measure std in pre-, over-
                        measure fluence in apertures around Point Sources

            assess std in pre- (~RON) is within allocated margins
            assess offsets in pre-, and over- are equal, within allocated  margins
            assess fluence is ~expected within apertures (PS) for each frame (pre-satur, satur, post-satur)


            plot point source fluence vs. OBSID, all sources
            [plot std vs. time]

            issue any warnings to log
            issue update to report          


        """
        raise NotImplementedError

    def prep_data(self):
        """

        **METACODE**

        ::

            Preparation of data for further analysis:

            f.e. ObsID [images with TPing only]:
                f.e.CCD:
                    f.e.Q:
                        subtract offset

        """

        raise NotImplementedError

    def basic_analysis(self):
        """

        Basic analysis of data.

        **METACODE**

        ::

            f.e.CCD:
                f.e.Q:
                    use SATURATED frame to generate pixel saturation MASK
                    measure stats in pix satur MASK across OBSIDs 
                     (pre-satur, satur, post-satur)


        """

        raise NotImplementedError

    def meta_analysis(self):
        """

        Meta-analysis of data.

        **METACODE**

        ::


            f.e.CCD:
                f.e.Q:
                    estimate delta-charge_0 and decay tau from time-series

            report:  
                persistence level (delta-charge_0) and time constant


        """
        raise NotImplementedError
