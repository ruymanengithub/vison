# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BF01

Brighter-Fatter Analysis
   Using data from test PTC01 or PTC02

Created on Wed Mar 7 10:57:00 2018

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

from vison.support import context
#from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import HKtools
from vison.datamodel import core, ccd
from vison.image import calibration
import ptc as ptclib
from vison.image import performance
#from FlatTask import FlatTask
from PTC0X import PTC0X
from vison.datamodel import inputs
#import BF01aux
from vison.analysis import Guyonnet15 as G15
from vison.image import covariance as covlib

from vison.support.files import cPickleRead, cPickleDumpDictionary
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


class BF01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
        ('Npix', ([int], 'Number of Pixels (linear) to consider for Covariance Matrix')),
    ])))


class BF01(PTC0X):
    """ """

    inputsclass = BF01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(BF01, self).__init__(inputs, log, drill, debug)
        self.name = 'BF01'
        #self.type = 'Simple'
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract_COV', self.extract_COV),
                         ('extract_BF', self.extract_BF),
                         ('meta', self.meta_analysis)]
        #self.HKKeys = HKKeys
        self.figdict = dict()  # BF01aux.gt_BF01figs(self.inputs['surrogate'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles')

    def set_inpdefaults(self, **kwargs):
        """ """

        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_inpdefaults(**_kwargs)
        self.inpdefaults['test'] = kwargs['test']

    def set_perfdefaults(self, **kwargs):

        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_perfdefaults(**_kwargs)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds PTC0X script structure dictionary.

        #:param exptimes: list of ints [ms], exposure times.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 800 nm.
        :param diffvalues: dict, opt, differential values.   

        """

        raise NotImplementedError(
            "%s: This Task does not build a script, it uses data from another test" % self.name)

    def filterexposures(self, structure, explogf, datapath, OBSID_lims):
        """

        """
        wavedkeys = ['motr_siz']
        return super(BF01, self).filterexposures(structure, explogf, datapath, OBSID_lims, colorblind=False,
                                                 wavedkeys=wavedkeys, surrogate=self.inputs['surrogate'])

    def prepare_images(self):
        super(BF01, self).prepare_images(doExtract=True, doMask=True,
                                         doOffset=True, doBias=False, doFF=False)

    def extract_COV(self):
        """

        Performs basic analysis of images:
            - extracts COVARIANCE matrix for each fluence

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extractCOV', Title='Covariance-Matrices Extraction', level=0)

        Npix = self.inputs['Npix']

        # labels should be the same accross CCDs. PATCH.
        label = self.dd.mx['label'][:, 0].copy()

        indices = copy.deepcopy(self.dd.indices)
        #nObs, nCCD, nQuad = indices.shape
        CCDs = indices[indices.names.index('CCD')].vals
        Quads = indices[indices.names.index('Quad')].vals

        if not self.drill:

            self.dd.products['COV'] = OrderedDict()

            dpath = self.inputs['subpaths']['ccdpickles']

            def fullinpath_adder(path): return os.path.join(dpath, path)
            vfullinpath_adder = np.vectorize(fullinpath_adder)

            ulabels = np.unique(label)

            for jCCD, CCD in enumerate(CCDs):

                self.dd.products['COV'][CCD] = OrderedDict()
                for Q in Quads:
                    self.dd.products['COV'][CCD][Q] = OrderedDict()

                for ulabel in ulabels:

                    six = np.where(label == ulabel)

                    ccdobjList = vfullinpath_adder(
                        self.dd.mx['ccdobj_name'][six[0], jCCD])

                    icovdict = covlib.get_cov_maps(
                        ccdobjList, Npix=Npix, doTest=False)

                    for Q in Quads:
                        self.dd.products['COV'][CCD][Q][ulabel] = copy.deepcopy(
                            icovdict[Q])

        # Plots
        # PENDING

        # REPORTING
        # PENDING

        return None

    def extract_BF(self):
        """ 

        Performs basic analysis of images:
            - extracts BF matrix for each COV matrix

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extractBF', Title='BF-Matrices Extraction', level=0)

        # labels should be the same accross CCDs. PATCH.
        label = self.dd.mx['label'][:, 0].copy()

        indices = copy.deepcopy(self.dd.indices)
        #nObs, nCCD, nQuad = indices.shape
        CCDs = np.array(indices[indices.names.index('CCD')].vals)
        Quads = np.array(indices[indices.names.index('Quad')].vals)

        ulabels = np.unique(label)

        # INITIALISATIONS

        kernel_FWHMx = np.zeros(
            (len(ulabels), len(CCDs), len(Quads)), dtype='float32')
        kernel_FWHMy = np.zeros(
            (len(ulabels), len(CCDs), len(Quads)), dtype='float32')
        kernel_e = np.zeros(
            (len(ulabels), len(CCDs), len(Quads)), dtype='float32')

        fluence = np.zeros(len(ulabels), dtype='float32') + np.nan

        self.dd.products['BF'] = OrderedDict()
        self.dd.products['BF']['ulabels'] = ulabels.copy()
        self.dd.products['BF']['CCDs'] = CCDs.copy()
        self.dd.products['BF']['Quads'] = Quads.copy()

        for jCCD, CCD in enumerate(CCDs):

            self.dd.products['BF'][CCD] = OrderedDict()

            for kQ, Q in enumerate(Quads):

                self.dd.products['BF'][CCD][Q] = OrderedDict()

            self.dd.products['BF']['kernel_FWHMx'] = kernel_FWHMx.copy()
            self.dd.products['BF']['kernel_FWHMy'] = kernel_FWHMy.copy()
            self.dd.products['BF']['kernel_e'] = kernel_e.copy()
            self.dd.products['BF']['fluence'] = fluence.copy()

        if not self.drill:

            singlepixmap = np.zeros((101, 101), dtype='float32') + 0.01
            singlepixmap[50, 50] = 1.

            for jCCD, CCD in enumerate(CCDs):

                for kQ, Q in enumerate(Quads):

                    for ix, ulabel in enumerate(ulabels):

                        COV_dict = self.dd.products['COV'][CCD][Q][ulabel]

                        COV_mx = COV_dict['av_covmap'].copy()

                        fluence[ix] = COV_dict['av_mu'].copy()

                        Asol_Q, psmooth_Q = G15.solve_for_A_linalg(
                            COV_mx, var=1., mu=1., returnAll=True, doplot=False)

                        kernel_Q = G15.degrade_estatic(singlepixmap, Asol_Q)

                        kerQshape = G15.get_cross_shape_rough(
                            kernel_Q, pitch=12.)

                        self.dd.products['BF'][CCD][Q][ulabel] = OrderedDict(Asol=Asol_Q.copy(),
                                                                             psmooth=psmooth_Q.copy(),
                                                                             kernel=kernel_Q.copy())

                        kernel_FWHMx[ix, jCCD, kQ] = kerQshape['FWHMx']
                        kernel_FWHMy[ix, jCCD, kQ] = kerQshape['FWHMy']
                        kernel_e[ix, jCCD, kQ] = kerQshape['e']

        self.dd.products['BF']['kernel_FWHMx'] = kernel_FWHMx.copy()
        self.dd.products['BF']['kernel_FWHMy'] = kernel_FWHMy.copy()
        self.dd.products['BF']['kernel_e'] = kernel_e.copy()
        self.dd.products['BF']['fluence'] = fluence.copy()

        # Plots
        # PENDING

        # REPORTING
        # PENDING

    def meta_analysis(self):
        """

        Analyzes the BF results across fluences.


        """

        raise NotImplementedError
