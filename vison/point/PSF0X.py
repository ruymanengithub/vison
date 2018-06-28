# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PSF0X

PSF vs. Fluence, and Wavelength
   PSF01 - nominal temperature
   PSF02 - alternative temperatures

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract offset level.
    - Divide by Flat-field.
    - Crop stamps of the sources on each CCD/Quadrant.
       - save snapshot figures of sources.
    - for each source:
       - measure shape using weighted moments
       - measure shape using Gaussian Fit
       - Bayesian Forward Modelling the optomechanic+detector PSF
    - Produce synoptic figures.
    - Save results.

Created on Thu Dec 29 15:01:07 2016

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

from vison.support import context
from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
#from vison.datamodel import HKtools
from vison.datamodel import scriptic as sc
#from vison.point import lib as polib
#from vison.support.report import Report
from vison.support import files
#from vison.pipe.task import Task
from PointTask import PointTask
import PSF0Xaux
from vison.image import performance
from vison.datamodel import inputs
from vison.support.files import cPickleRead
# END IMPORT

isthere = os.path.exists

stampw = polib.stampw

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

PSF0X_commvalues = dict(program='CALCAMP',
                        IPHI1=1, IPHI2=1, IPHI3=1, IPHI4=0,
                        rdmode='fwd_bas',
                        flushes=7,
                        exptime=0.,
                        vstart=0,
                        vend=2086,
                        shuttr=1,
                        siflsh=1,
                        siflush_p=500,
                        motr_on=1,
                        motr_cnt=2,
                        motr_siz=120,
                        source='point',
                        mirr_on=1,
                        wave=4,
                        mirr_pos=ogsemod.mirror_nom['F4'],
                        comments='')

PSF0X_relfluences = np.array([5., 25., 50., 75., 90.])


def get_testdefaults(ogseobj=None):

    if ogseobj is None:
        ogseobj = ogsemod.Ogse()

    testdefaults = dict(waves=[590, 640, 730, 800, 880, 0],
                        exptimes=dict(),
                        frames=[20, 14, 10, 4, 4])

    for w in testdefaults['waves']:
        tFWC_pointw = ogseobj.profile['tFWC_point']['nm%i' % w]
        testdefaults['exptimes']['nm%i' % w] = \
            (PSF0X_relfluences/100.*tFWC_pointw).tolist()

    return testdefaults


def get_Flu_lims(relfluences):

    plusminus10pc = 1.+np.array([-0.1, 0.1])
    satur_fluence = 2.*2.E5 / performance.gains['CCD1']

    Flu_lims = OrderedDict(
        CCD1=OrderedDict(
            E=OrderedDict(
                ALPHA=OrderedDict())))  # +/-10%

    Nfluences = len(relfluences)
    for i in range(1, Nfluences+1):
        Flu_lims['CCD1']['E']['ALPHA']['col%i' % i] = \
            2. * relfluences[i-1]*plusminus10pc*satur_fluence

    for Spot in ['BRAVO', 'CHARLIE', 'DELTA', 'ECHO']:
        Flu_lims['CCD1']['E'][Spot] = Flu_lims['CCD1']['E']['ALPHA']
    for Q in ['F', 'G', 'H']:
        Flu_lims['CCD1'][Q] = copy.deepcopy(Flu_lims['CCD1']['E'])
    for CCD in [2, 3]:
        Flu_lims['CCD%i' % CCD] = copy.deepcopy(Flu_lims['CCD1'])

    return Flu_lims


FWHM_lims = OrderedDict(CCD1=OrderedDict(
    E=OrderedDict(
        ALPHA=[1.2, 2.5])))
for Spot in ['BRAVO', 'CHARLIE', 'DELTA', 'ECHO']:
    FWHM_lims['CCD1']['E'][Spot] = FWHM_lims['CCD1']['E']['ALPHA']
for Q in ['F', 'G', 'H']:
    FWHM_lims['CCD1'][Q] = copy.deepcopy(FWHM_lims['CCD1']['E'])
for CCD in [2, 3]:
    FWHM_lims['CCD%i' % CCD] = copy.deepcopy(FWHM_lims['CCD1'])

BGD_lims = OrderedDict(CCD1=OrderedDict(
    E=[-1., 5.]))
for Q in ['F', 'G', 'H']:
    BGD_lims['CCD1'][Q] = copy.deepcopy(BGD_lims['CCD1']['E'])
for CCD in [2, 3]:
    BGD_lims['CCD%i' % CCD] = copy.deepcopy(BGD_lims['CCD1'])


class PSF0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
    ])))


class PSF0X(PointTask):

    inputsclass = PSF0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(PSF0X, self).__init__(inputs, log, drill, debug)
        self.name = 'PSF0X'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis), 
                         ('bayes', self.bayes_analysis),
                         ('meta', self.meta_analysis)]
        for v in self.subtasks: self.inputs['todo_flags'][v[0]] = False
        self.HKKeys = HKKeys
        self.figdict = PSF0Xaux.gt_PSF0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles')
        self.init_todo_flags()
        

    def set_inpdefaults(self, **kwargs):

        testkey = kwargs['test']
        testdefaults = get_testdefaults(self.ogse)

        if 'PSF01' in testkey or 'PSFLUX00' in testkey:
            wavelength = kwargs['wavelength']
        elif 'PSF02' in testkey:
            wavelength = 800

        exptimes = testdefaults['exptimes']['nm%i' % wavelength]
        frames = testdefaults['frames']

        self.inpdefaults = dict(wavelength=wavelength,
                                frames=frames,
                                exptimes=exptimes)

    def set_perfdefaults(self, **kwargs):
        super(PSF0X, self).set_perfdefaults(**kwargs)

        self.perfdefaults['BGD_lims'] = BGD_lims  # ADUs
        self.perfdefaults['Flu_lims'] = get_Flu_lims(PSF0X_relfluences)  # ADUs
        self.perfdefaults['FWHM_lims'] = FWHM_lims  # Pixels

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ 

        Builds PSF0X script structure dictionary.

        #:param exptimes: list of ints, [ms], exposure times.
        #:param frames: list of frame numbers. Same length as exptimes.
        #:param wavelength: int, [nm], wavelength.
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """
        exptimes = self.inputs['exptimes']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        test = self.inputs['test']

        assert len(exptimes) == len(frames)

        FW_ID = self.ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])

        PSF0X_commvalues['wave'] = FW_IDX
        PSF0X_commvalues['mirr_pos'] = self.ogse.profile['mirror_nom']['F%i' % FW_IDX]

        ncols = len(exptimes)

        PSF0X_sdict = dict()

        for ic in range(ncols):
            colid = 'col%i' % (ic+1,)
            PSF0X_sdict[colid] = dict(frames=frames[ic], exptime=exptimes[ic],
                                      test=test)

        Ncols = len(PSF0X_sdict.keys())
        PSF0X_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(PSF0X_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        PSF0X_sdict = sc.update_structdict(PSF0X_sdict, commvalues, diffvalues)

        return PSF0X_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """
        """
        wavedkeys = []
        return super(PSF0X, self).filterexposures(structure, explog, OBSID_lims, colorblind=False,
                                                  wavedkeys=wavedkeys)

    def prep_data(self):
        """

        PSF0X: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies (per CCD  & Q):
            cosmetics masking 
            offset subtraction
            bias subtraction
            FF division

            per spot:               
                cuts-out and save stamps of pre-processed spots for further analysis.


        """
        super(PSF0X, self).prepare_images(
            doExtract=True, doMask=True, doOffset=True, doBias=True, doFF=True)

        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices[dIndices.names.index('CCD')].vals
        #nCCD = len(CCDs)
        Quads = dIndices[dIndices.names.index('Quad')].vals
        nQuads = len(Quads)
        nObs = dIndices[dIndices.names.index('ix')].len
        SpotNames = dIndices[dIndices.names.index('Spot')].vals
        nSpots = len(SpotNames)

        CIndices = copy.deepcopy(dIndices)
        CIndices.pop(CIndices.names.index('Quad'))
        CIndices.pop(CIndices.names.index('Spot'))

        # INITIALIZATIONS

        self.dd.initColumn('spots_name', CIndices, dtype='S100', valini='None')

        if not self.drill:

            #rpath = self.inputs['resultspath']
            picklespath = self.inputs['subpaths']['ccdpickles']
            spotspath = self.inputs['subpaths']['spots']

            for iObs in range(nObs):

                for jCCD, CCD in enumerate(CCDs):

                    CCDkey = 'CCD%i' % CCD

                    fullccdobj_name = os.path.join(
                        picklespath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    # Cut-out "spots"

                    spots_array = np.zeros((nQuads, nSpots), dtype=object)

                    for kQ, Quad in enumerate(Quads):

                        for lS, SpotName in enumerate(SpotNames):

                            coo = polib.Point_CooNom[CCDkey][Quad][SpotName]
                            lSpot = polib.extract_spot(
                                ccdobj, coo, Quad, stampw=stampw)

                            spots_array[kQ, lS] = copy.deepcopy(lSpot)

                    # save "spots" to a separate file and keep name in dd

                    self.dd.mx['spots_name'][iObs,
                                             jCCD] = '%s_spots' % self.dd.mx['File_name'][iObs, jCCD]

                    fullspots_name = os.path.join(
                        spotspath, '%s.pick' % self.dd.mx['spots_name'][iObs, jCCD])

                    spotsdict = dict(spots=spots_array)

                    files.cPickleDumpDictionary(spotsdict, fullspots_name)

        if self.log is not None:
            self.log.info('Saved spot "bag" files to %s' % spotspath)

    def basic_analysis(self):
        """Performs basic analysis on spots:
             - shape from moments
             - Gaussian fitting: peak intensity, position, width_x, width_y

        """

        raise NotImplementedError

    # def bayes_analysis(dd,report,inputs,log=None):

    def bayes_analysis(self):
        """ 
        Performs bayesian decomposition of the spot images:
            - optomechanic PSF and detector PSF.
        Also measures R2, fwhm and ellipticity of "extracted" detector PSF.
        """
        raise NotImplementedError

    # def meta_analysis(dd,report,inputs,log=None):
    def meta_analysis(self):
        """

        Analyzes the relation between detector PSF and fluence.

        """
        raise NotImplementedError
