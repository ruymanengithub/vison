#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FLAT0X

Flat-fields acquisition / analysis script

Created on Tue Aug 29 17:32:52 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
import copy
from collections import OrderedDict

from vison.support import context
from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
#from vison.pipe.task import Task
from vison.flat.FlatTask import FlatTask
from vison.image import performance
from vison.datamodel import inputs
import FLAT0Xaux as FL0Xaux
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

FLAT0X_commvalues = dict(program='CALCAMP',
                         IPHI1=1, IPHI2=1, IPHI3=1, IPHI4=0,
                         rdmode='fwd_bas',
                         flushes=7, vstart=0, vend=2086,
                         exptime=0., shuttr=1,
                         siflush=0,
                         wave=4,
                         comments='')

FLU_lims = dict(CCD1=dict(
    col1=0.25 * 2**16 * (1.+np.array([-0.10, 0.10])),
    col2=0.50 * 2**16 * (1.+np.array([-0.10, 0.10])),
    col3=0.75 * 2**16 * (1.+np.array([-0.10, 0.10]))))
for i in [2, 3]:
    FLU_lims['CCD%i' % i] = copy.deepcopy(FLU_lims['CCD1'])


class FLATS0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
    ])))


class FLAT0X(FlatTask):
    """ """

    inputsclass = FLATS0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(FLAT0X, self).__init__(inputs, log, drill, debug)
        self.name = 'FLAT0X'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('indivflats', self.do_indiv_flats),
                         ('masterflat', self.do_master_flat),
                         ('prmask', self.do_prdef_mask)]
        self.HKKeys = HKKeys
        self.figdict = FL0Xaux.gt_FL0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       ccdflats='ccdflats')

    def set_inpdefaults(self, **kwargs):
        """ """

        try:
            wavelength = kwargs['wavelength']
        except KeyError:
            wavelength = 800
        try:
            test = kwargs['test']
        except KeyError:
            test = 'FLAT0X'

        t_dummy_F0X = np.array([25., 50., 75])/100.
        exptimesF0X = (ogse.tFWC_flat['nm%i' %
                                      wavelength] * t_dummy_F0X).tolist()  # s
        framesF0X = [80, 60, 30]

        self.inpdefaults = dict(exptimes=exptimesF0X,
                                frames=framesF0X,
                                wavelength=wavelength,
                                test=test)

    def set_perfdefaults(self, **kwargs):
        #wavelength = self.inputs['wavelength']
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        self.perfdefaults['FLU_lims'] = FLU_lims.copy()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds FLAT0X script structure dictionary.

        :param diffvalues: dict, opt, differential values.

        """

        exptimes = self.inputs['exptimes']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        test = self.inputs['test']

        FW_ID = ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])

        FLAT0X_commvalues['wave'] = FW_IDX
        FLAT0X_commvalues['test'] = test

        assert len(exptimes) == len(frames)
        #assert len(exptimes) == len(flags)

        FLAT0X_sdict = dict()
        for i, exptime in enumerate(exptimes):
            FLAT0X_sdict['col%i' % (i+1,)] = dict(frames=frames[i], exptime=exptimes[i],
                                                  comments='EXP%.1e' % exptime)  # ,comments=flags[i])

        Ncols = len(FLAT0X_sdict.keys())
        FLAT0X_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(FLAT0X_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        FLAT0X_sdict = sc.update_structdict(
            FLAT0X_sdict, commvalues, diffvalues)

        return FLAT0X_sdict

    def filterexposures(self, structure, explogf, datapath, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(FLAT0X, self).filterexposures(structure, explogf, datapath, OBSID_lims, colorblind=True,
                                                   wavedkeys=wavedkeys)

    def prepare_images(self):
        """

        FLAT0X: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            [bias structure subtraction, if available]
            cosmetics masking

        """
        super(FLAT0X, self).prepare_images(doExtract=True,
                                           doMask=True, doOffset=True, doBias=True)

    def do_indiv_flats(self):
        """

        **METACODE**

        ::

            Preparation of data for further analysis and 
            produce flat-field for each OBSID.

            f.e. ObsID:
                f.e.CCD:

                    load ccdobj

                    f.e.Q:

                        model 2D fluence distro in image area
                        produce average profile along rows
                        produce average profile along cols

                    save 2D model and profiles in a pick file for each OBSID-CCD
                    divide by 2D model to produce indiv-flat
                    save indiv-Flat to FITS(?), update add filename

            plot average profiles f. each CCD and Q (color coded by time)

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='indivFF', Title='Individual Flat-Fields', level=0)

        # INITIALISATION

        indices = copy.deepcopy(self.dd.indices)

        ObsIDs = indices[indices.names.index('ObsID')].vals
        CCDs = indices[indices.names.index('CCD')].vals

        ccdindices = copy.deepcopy(self.dd.mx['CCD'].indices)

        self.dd.initColumn('indiv_flats', ccdindices,
                           dtype='S100', valini='None')

        # The heavy-lifting

        if not self.drill:

            dpath = self.inputs['subpaths']['ccdpickles']
            fpath = self.inputs['subpaths']['ccdflats']

            def fullinpath_adder(path): return os.path.join(dpath, path)
            vfullinpath_adder = np.vectorize(fullinpath_adder)
            fccdobj_names = vfullinpath_adder(self.dd.mx['ccdobj_name'])

            wavelength = self.inputs['wavelength']

            for iObs, ObsID in enumerate(ObsIDs):

                for jCCD, CCD in enumerate(CCDs):
                    ilabel = self.dd.mx['label'][iObs, jCCD]

                    self.dd.mx['indiv_flats'][iObs, jCCD] = 'EUC_FF_%inm_%s_ROE1_%s.fits' %\
                        (wavelength, ilabel, ObsID, CCD)

            def fulloutpath_adder(path): return os.path.join(fpath, path)
            vfulloutpath_adder = np.vectorize(fulloutpath_adder)
            findiv_flats = vfulloutpath_adder(self.dd.mx['indiv_flats'])

            ffsettings = dict()

            FFing.produce_IndivFlats(fccdobj_names.flatten(), findiv_flats.flatten(),
                                     settings=ffsettings,
                                     runonTests=False,
                                     processes=self.processes)

            # MISSING: 1D profiles. Do as part of produce_IndivFlats, or separate?

            # Show 1D Profiles for each fluence across CCDs:
            #   2 profiles x N-fluences = 2N plots (each with 4Qx3CCD subplots)

    def do_master_flat(self):
        """ 

        **METACODE**

        ::

            Produces Master Flat-Field

            f.e.CCD:
                f.e.Q:
                    stack individual flat-fields by chosen estimator
            save Master FF to FITS
            measure PRNU and 
            report PRNU figures

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='MasterFF', Title='Master Flat-Fields', level=0)

        wavelength = self.inputs['wavelength']
        settings = dict()

        indices = copy.deepcopy(self.dd['indiv_flats'].indices)
        #ObsIDs = indices[indices.names.index('ObsID')].vals
        CCDs = indices[indices.names.index('CCD')].vals

        dpath = self.inputs['subpaths']['ccdflats']
        cdppath = self.inputs['subpaths']['cdps']

        labels = self.dd.mx['label'][:].copy()
        ulabels = np.unique(labels)

        if not self.drill:

            PRNU = OrderedDict()

            self.dd.products['MasterFFs'] = OrderedDict()

            for ulabel in ulabels:

                PRNU[ulabel] = OrderedDict()

                self.dd.products['MasterFFs'][ulabel] = OrderedDict()

                for jCCD, CCD in enumerate(CCDs):

                    PRNU[ulabel][CCD] = OrderedDict()

                    FFname = 'EUC_FF_%inm_%s_ROE1_%s.fits' % \
                        (wavelength, ulabel, CCD)

                    FFpath = os.path.join(cdppath, FFname)

                    selix = np.where((labels == ulabel) & (CCDs == CCD))

                    FFlist = self.dd.mx['indiv_flats'][selix].flatten().copy()

                    def fullinpath_adder(
                        path): return os.path.join(dpath, path)
                    vfullinpath_adder = np.vectorize(fullinpath_adder)
                    FFlist = vfullinpath_adder(FFlist)

                    # MISSING: proper defects and useful area masking
                    #   (mask-out pre/over scans)

                    FFing.produce_MasterFlat(
                        FFlist, FFpath, mask=None, settings=settings)

                    self.dd.products['MasterFFs'][ulabel][CCD] = FFpath

                    # Measure Flat-Field (PRNU): per-Q

                    FF = FFing.FlatField(FFpath)

                    iFextension = FF.extnames.index('FLAT')

                    Quads = FF.Quads

                    for Q in Quads:

                        iQ_PRNU = FF.get_stats(Q, sector='img', statkeys=['std'],
                                               ignore_pover=True,
                                               extension=iFextension)[0]

                        PRNU[ulabel][CCD][Q] = iQ_PRNU

        # REPORT PRNU results

        # SHOW FFs

    def do_prdef_mask(self):
        """
        **METACODE**

        ::

            Produces mask of defects in Photo-Response
            Could use master FF, or a stack of a subset of images (in order
            to produce mask, needed by other tasks, quicker).

            f.e.CCD:
                f.e.Q:
                    produce mask of PR defects
                    save mask of PR defects
                    count dead pixels / columns 

            report PR-defects stats

        """

        raise NotImplementedError
