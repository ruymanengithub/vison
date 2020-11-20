# -*- coding: utf-8 -*-
"""

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

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import string as st
import warnings
import copy
from collections import OrderedDict
from skimage import exposure
from sklearn import linear_model
from astropy.nddata import Cutout2D
from astropy.modeling import models
from scipy import stats

from vison.support import utils
from vison.pipe.task import HKKeys
from vison.support import context
from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
from vison.datamodel import scriptic as sc
from vison.support import files
from vison.point import PointTask as PT
from . import PSF0Xaux
from vison.support.files import cPickleRead, cPickleDumpDictionary
from vison.datamodel import ccd
from vison.xtalk import opt_xtalk as oxt
from vison.xtalk import xtalk as xt
from vison.analysis import Guyonnet15 as G15
# END IMPORT

isthere = os.path.exists

stampw = polib.stampw

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

PSF0X_commvalues = dict(program='CALCAMP',
                        flushes=7, siflsh=1, siflsh_p=500,
                        inisweep=1,
                        vstart=0, vend=2086,
                        toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                        chinj=0,
                        s_tpump=0,
                        v_tpump=0,
                        shuttr=1, 
                        e_shuttr=0,
                        mirr_on=1,
                        wave=4,
                        mirr_pos=ogsemod.mirror_nom['F4'],
                        motr_on=1,
                        motr_cnt=2,
                        motr_siz=120,
                        source='point',
                        comments='')

PSF0X_relfluences = np.array([5., 25., 50., 75., 90.])


def get_testdefaults(ogseobj=None):

    if ogseobj is None:
        ogseobj = ogsemod.Ogse()

    testdefaults = dict(waves=[590, 730, 800, 880, 0],
                        exptimes=dict(),
                        frames=[20, 15, 10, 4, 4])

    for w in testdefaults['waves']:
        tFWC_pointw = ogseobj.profile['tFWC_point']['nm%i' % w]
        testdefaults['exptimes']['nm%i' % w] = \
            (PSF0X_relfluences / 100. * tFWC_pointw).tolist()

    return testdefaults


def get_PeakFlu_lims(relfluences):

    plusminus30pc = 1. + np.array([-0.3, 0.3])
    satur_fluence = context.eff_satur

    # assuming a best fwhm~2.5pix (default), and a gaussian profile
    # F = 2pi(fwhm/2.355)**2*I0

    #I02F = 2*np.pi*(fwhm/2.355)**2.

    PeakFlu_lims = OrderedDict(
        CCD1=OrderedDict(
            E=OrderedDict(
                ALPHA=OrderedDict())))  # +/-10%

    Nfluences = len(relfluences)
    for i in range(1, Nfluences + 1):
        relflu = min(relfluences[i - 1] / 100., 1)

        PeakFlu_lims['CCD1']['E']['ALPHA']['col%03i' % i] = \
            relflu * plusminus30pc * satur_fluence
#            I02F * relflu * plusminus30pc * satur_fluence

    for Spot in ['BRAVO', 'CHARLIE', 'DELTA', 'ECHO']:
        PeakFlu_lims['CCD1']['E'][Spot] = PeakFlu_lims['CCD1']['E']['ALPHA']
    for Q in ['F', 'G', 'H']:
        PeakFlu_lims['CCD1'][Q] = copy.deepcopy(PeakFlu_lims['CCD1']['E'])
    for CCD in [2, 3]:
        PeakFlu_lims['CCD%i' % CCD] = copy.deepcopy(PeakFlu_lims['CCD1'])

    return PeakFlu_lims


FWHM_lims = OrderedDict(CCD1=OrderedDict(
    E=OrderedDict(
        ALPHA=[1.2, 3.])))
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


class PSF0X_inputs(PT.Point_inputs):
    manifesto = PT.Point_inputs.manifesto.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
    ])))


class PSF0X(PT.PointTask):

    inputsclass = PSF0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('lock', self.lock_on_stars),
                         ('relock', self.relock),
                         ('check', self.check_data),
                         ('prep', self.prep_data),
                         ('corrBFE_G15', self.corrBFE_G15),
                         ('basic', self.basic_analysis),
                         ('basic_nobfe', self.basic_analysis_noBFE),
                         ('bayes', self.bayes_analysis),
                         ('meta', self.meta_analysis),
                         ('xtalk_sex', self.opt_xtalk_sextract),
                         ('xtalk_build', self.opt_xtalk_build),
                         ('xtalk_meta', self.opt_xtalk_meta),
                         ('debugtask', self.debugtask)]

        super(PSF0X, self).__init__(inputs=inputs, log=log, drill=drill,
                                    debug=debug, cleanafter=cleanafter)
        self.name = 'PSF0X'
        self.type = 'Simple'

        self.HKKeys = HKKeys
        self.CDP_lib = PSF0Xaux.get_CDP_lib(self.inputs['test'])
        self.figdict = PSF0Xaux.get_PSF0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       products='products', spots='spots',
                                       xtalk='xtalk')
        
        if 'inCDPs' in self.inputs:
            if isinstance(self.inputs['inCDPs'],dict):

                try:
                    selectFF = ('wavelength' in self.inputs) and \
                    ('FF' in self.inputs['inCDPs'])
                except TypeError:
                    selectFF = False

                if selectFF:
                    wavelength = self.inputs['wavelength']
                    nmkey = 'nm%i' % wavelength
                    nmFFs = self.inputs['inCDPs']['FF'][nmkey]
                    self.inputs['inCDPs']['FF'] = nmFFs.copy()


    def set_inpdefaults(self, **kwargs):

        testkey = kwargs['test']
        testdefaults = get_testdefaults(self.ogse)

        if 'PSF01' in testkey or 'PSFLUX00' in testkey:
            wavelength = kwargs['wavelength']
        elif 'PSF02' in testkey:
            wavelength = 800

        exptimes = testdefaults['exptimes']['nm%i' % wavelength]
        frames = testdefaults['frames']

        self.inpdefaults = dict(
            offsetxy=[0., 0.],
            wavelength=wavelength,
            frames=frames,
            exptimes=exptimes)

    def set_perfdefaults(self, **kwargs):
        super(PSF0X, self).set_perfdefaults(**kwargs)

        self.perfdefaults['BGD_lims'] = BGD_lims  # ADUs
        self.perfdefaults['PeakFlu_lims'] = get_PeakFlu_lims(PSF0X_relfluences)  # ADUs
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
            colid = 'col%03i' % (ic + 1,)
            PSF0X_sdict[colid] = dict(frames=frames[ic], exptime=exptimes[ic],
                                      test=test)

        Ncols = len(list(PSF0X_sdict.keys()))
        PSF0X_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(PSF0X_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
                diffvalues = diffvalues = dict()

        PSF0X_sdict = sc.update_structdict(PSF0X_sdict, commvalues, diffvalues)

        return PSF0X_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """
        """
        wavedkeys = []
        return super(PSF0X, self).filterexposures(structure, explog, OBSID_lims, colorblind=False,
                                                  wavedkeys=wavedkeys)

    def lock_on_stars(self):
        """ """

        FW_ID = self.dd.mx['wave'][0, 0]
        wave = self.ogse.get_wavelength(FW_ID)
        tFWC_point = self.ogse.profile['tFWC_point']['nm%i' % wave]
        exptime = self.dd.mx['exptime'][:, 0]

        sexconfig = dict(MINAREA=2., DET_THRESH=15.,
                         MAG_ZEROPOINT=20.)

        # single Obsid-locking
        #iObs = np.abs(exptime-tFWC_point/2.).argmin()
        # PT.PointTask.lock_on_stars(self,iObs=iObs,
        #                           sexconfig=sexconfig)

        ulabels = np.unique(self.dd.mx['label'][:, 0]).tolist()

        # ulabels = [ulabels[0]] # TEST

        ObsList = []

        for ulabel in ulabels:

            ixsel = np.where(self.dd.mx['label'][:, 0] == ulabel)[0][0]
            ObsList.append(ixsel)

        PT.PointTask.lock_on_stars(self, 
            iObs=ObsList,
            labels=ulabels,
            sexconfig=sexconfig)

        self._aprox_missing_locks()

    def relock(self):
        """ """
        PT.PointTask.relock(self, check_ok=False)
        self._aprox_missing_locks()


    def _aprox_missing_locks(self):
        """ """
        lock_tb_cdp = files.cPickleRead(self.dd.products['LOCK_TB_CDP'])
        CCDs = list(self.ogse.startrackers.keys())
        ulabels = list(self.ogse.startrackers[CCDs[0]].keys())

        enough_stars = 10

        lock_tb = lock_tb_cdp['data']['LOCK_TB'].copy()


        for iCCD,CCDk in enumerate(CCDs):
            for k, ulabel in enumerate(ulabels):

                jj = np.where((lock_tb.CCD==iCCD+1) & (lock_tb.LABEL==k))[0][0]

                if lock_tb.NMATCH[jj] < enough_stars:
                    ixgood = np.where((lock_tb.CCD==iCCD+1) &\
                     (lock_tb.NMATCH>=enough_stars))[0]
                    jjgood = ixgood[np.abs(ixgood-jj).argmin()]

                    kgood = lock_tb.LABEL[jjgood]
                    ulabelgood = ulabels[kgood]

                    self.ogse.startrackers[CCDk][ulabel] = \
                        copy.deepcopy(self.ogse.startrackers[CCDk][ulabelgood])

                    if self.log is not None:
                        msg = 'WARNING: lock of %s/%s replaced with %s/%s' % \
                            (CCDk, ulabel, CCDk, ulabelgood)
                        self.log.info(msg)

                else:
                    pass

    def debugtask(self):
        """ """

        if self.report is not None:
            self.report.add_Section(
                keyword='debug', Title='Debugging', level=0)


        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices.get_vals('CCD')
        nCCDs = len(CCDs)
        Quads = dIndices.get_vals('Quad')
        nQuads = len(Quads)
        nObs = dIndices[dIndices.names.index('ix')].len
        SpotNames = dIndices.get_vals('Spot')
        nSpots = len(SpotNames)

        CIndices = copy.deepcopy(dIndices)
        CIndices.pop(CIndices.names.index('Quad'))
        CIndices.pop(CIndices.names.index('Spot'))

        spotspath = self.inputs['subpaths']['spots']
        picklespath = self.inputs['subpaths']['ccdpickles']

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        spots_cdp = self.CDP_lib['SPOTS']

        spots_cdp.rootname = '{}_{}_{}'.format(spots_cdp.rootname,
                self.inputs['test'],
                self.inputs['BLOCKID'])

        spots_cdp.path = spotspath

        spots_cdp.header = CDP_header.copy()
        spots_cdp.meta = dict()
        spots_cdp.data = OrderedDict(spots=np.zeros((nObs, nCCDs, nQuads, nSpots), dtype=object))

        strackers = self.ogse.startrackers
        stlabels = self.ogse.labels

        psCCDcoodicts = OrderedDict(names=strackers['CCD1']['col001'].starnames,
            CCDs=CCDs,
            labels=stlabels)

        for jCCD, CCDk in enumerate(CCDs):
            psCCDcoodicts[CCDk] = OrderedDict()
            for ilabel, label in enumerate(stlabels):
                psCCDcoodicts[CCDk][label] = strackers[CCDk][label].get_allCCDcoos(
                    nested=True)

        for iObs in range(nObs):
                # for iObs in range(3): # TESTS

                for jCCD, CCDk in enumerate(CCDs):

                    ilabel = self.dd.mx['label'][iObs, jCCD]

                    fullccdobj_name = os.path.join(
                        picklespath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    # Cut-out "spots"

                
                    for kQ, Quad in enumerate(Quads):

                        for lS, SpotName in enumerate(SpotNames):

                            #coo = polib.Point_CooNom[CCDk][Quad][SpotName]
                            coo = psCCDcoodicts[CCDk][ilabel][Quad][SpotName]
                            lSpot = polib.extract_spot(
                                ccdobj, coo, Quad, stampw=stampw)

                            spots_cdp.data['spots'][iObs, jCCD, kQ, lS] = copy.deepcopy(lSpot)


        spots_cdp.meta.update(dict(nObs=nObs,
            nCCDs=nCCDs,
            nQuads=nQuads,
            nSpots=nSpots,
            SpotNames=SpotNames,
            structure='spots:np.zeros((nObs, nCCDs, nQuads, nSpots), dtype=object)'))

        self.save_CDP(spots_cdp)
        self.pack_CDP_to_dd(spots_cdp, 'SPOTS')


        if self.log is not None:
            self.log.info('Saved spot "bag" files to %s' % spotspath)

        

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

        
        onTests = False

        if onTests:

            for colname in self.dd.colnames:
                if colname != 'ObsID':
                    self.dd.mx[colname].array = np.expand_dims(
                        self.dd.mx[colname].array[37, ...], 0).copy()
                else:
                    self.dd.mx[colname] = self.dd.mx[colname].array[37].copy()
            self.dd.indices[0].vals = [0]
            self.dd.indices[0].len = 1


        bypass = False
        if not bypass: #TEST

            super(PSF0X, self).prepare_images(
                doExtract=True, 
                doBadPixels=True,
                doMask=True, 
                doOffset=True, 
                doBias=False, 
                doFF=True)

        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices.get_vals('CCD')
        nCCDs = len(CCDs)
        Quads = dIndices.get_vals('Quad')
        nQuads = len(Quads)
        nObs = dIndices[dIndices.names.index('ix')].len
        SpotNames = dIndices.get_vals('Spot')
        nSpots = len(SpotNames)

        CIndices = copy.deepcopy(dIndices)
        CIndices.pop(CIndices.names.index('Quad'))
        CIndices.pop(CIndices.names.index('Spot'))

        # INITIALIZATIONS

        #self.dd.initColumn('spots_name', CIndices, dtype='U100', valini='None')

        if not self.drill:

            function, module = utils.get_function_module()
            CDP_header = self.CDP_header.copy()
            CDP_header.update(dict(function=function, module=module))
            CDP_header['DATE'] = self.get_time_tag()

            #rpath = self.inputs['resultspath']
            picklespath = self.inputs['subpaths']['ccdpickles']
            spotspath = self.inputs['subpaths']['spots']

            spots_cdp = self.CDP_lib['SPOTS']

            spots_cdp.rootname = '{}_{}_{}'.format(spots_cdp.rootname,
                self.inputs['test'],
                self.inputs['BLOCKID'])

            spots_cdp.path = spotspath

            spots_cdp.header = CDP_header.copy()
            spots_cdp.meta = dict()
            spots_cdp.data = OrderedDict(spots=np.zeros((nObs, nCCDs, nQuads, nSpots), dtype=object))


            strackers = self.ogse.startrackers
            stlabels = self.ogse.labels

            psCCDcoodicts = OrderedDict(names=strackers['CCD1']['col001'].starnames,
                CCDs=CCDs,
                labels=stlabels)

            for jCCD, CCDk in enumerate(CCDs):
                psCCDcoodicts[CCDk] = OrderedDict()
                for ilabel, label in enumerate(stlabels):

                    psCCDcoodicts[CCDk][label] = strackers[CCDk][label].get_allCCDcoos(
                        nested=True)

            for iObs in range(nObs):
                # for iObs in range(3): # TESTS

                for jCCD, CCDk in enumerate(CCDs):

                    ilabel = self.dd.mx['label'][iObs, jCCD]

                    fullccdobj_name = os.path.join(
                        picklespath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    # Cut-out "spots"

                    for kQ, Quad in enumerate(Quads):

                        for lS, SpotName in enumerate(SpotNames):

                            #coo = polib.Point_CooNom[CCDk][Quad][SpotName]
                            coo = psCCDcoodicts[CCDk][ilabel][Quad][SpotName]
                            lSpot = polib.extract_spot(
                                ccdobj, coo, Quad, stampw=stampw)

                            spots_cdp.data['spots'][iObs, jCCD, kQ, lS] = copy.deepcopy(lSpot)


                    if onTests:
                    #doPlot=True
                    #if doPlot and iObs==3 and CCDk=='CCD1':

                        from matplotlib import pyplot as plt

                        for kQ, Quad in enumerate(Quads):

                            for lS, SpotName in enumerate(SpotNames):

                                img = spots_cdp.data['spots'][iObs, jCCD, kQ, lS].data
                                fig = plt.figure()
                                ax = fig.add_subplot(111)
                                ax.imshow(img.T, origin='lower left')
                                ax.set_title('Obs%i, %s-%s, %s' % \
                                    (iObs+1,CCDk,Quad,SpotName))
                                plt.show()

        spots_cdp.meta.update(dict(nObs=nObs,
            nCCDs=nCCDs,
            nQuads=nQuads,
            nSpots=nSpots,
            SpotNames=SpotNames,
            structure='spots:np.zeros((nObs, nCCDs, nQuads, nSpots), dtype=object)'))

        self.save_CDP(spots_cdp)
        self.pack_CDP_to_dd(spots_cdp, 'SPOTS')


        if self.log is not None:
            self.log.info('Saved spot "bag" files to %s' % spotspath)

    def _get_psCCDcoodicts(self):

        CCDs = self.dd.indices.get_vals('CCD')
        strackers = self.ogse.startrackers
        stlabels = self.ogse.labels

        psCCDcoodicts = OrderedDict(names=strackers['CCD1']['col001'].starnames,
            CCDs=CCDs,
            labels=stlabels)

        for jCCD, CCDk in enumerate(CCDs):
            psCCDcoodicts[CCDk] = OrderedDict()
            for ilabel, label in enumerate(stlabels):
                psCCDcoodicts[CCDk][label] = \
                    strackers[CCDk][label].get_allCCDcoos(
                    nested=True)

        return psCCDcoodicts


    def corrBFE_G15(self):
        """Corrects BFE using the matrices from Guyonnet et al. 2015."""

        if self.report is not None:
            self.report.add_Section(
                keyword='corrG15', Title='Correction of BFE using G15', level=0)


        if 'BFE' not in self.inputs['inCDPs'].keys():

            if self.log is not None:
                self.log.info('Cannot correct BFE without BFE CDP!')

            if self.reports is not None:
                self.report.add_Text('Cannot correct BFE without BFE CDP!')

            return None

        # Loading the BFE model


        wavelength = self.inputs['wavelength']


        BFEddpick = self.inputs['inCDPs']['BFE']['nm%i' % wavelength]
        BFEall = files.cPickleRead(BFEddpick).products['BF'].copy()
        ulabels = BFEall['ulabels']
        midrangekey = ulabels[len(ulabels)//2-1]
        inCCDs = BFEall['CCDs']
        inQuads = BFEall['Quads']


        Asols = dict()
        for inCCD in inCCDs:
            Asols[inCCD] = dict()
            for Q in inQuads:
                Asols[inCCD][Q] = BFEall[inCCD][Q][midrangekey]['Asol'].copy()

        
        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices.get_vals('CCD')
        nCCDs = len(CCDs)
        Quads = dIndices.get_vals('Quad')
        nQuads = len(Quads)
        nObs = dIndices[dIndices.names.index('ix')].len
        SpotNames = dIndices.get_vals('Spot')
        nSpots = len(SpotNames)

        CIndices = copy.deepcopy(dIndices)
        CIndices.pop(CIndices.names.index('Quad'))
        CIndices.pop(CIndices.names.index('Spot'))

        # INITIALIZATIONS

        spotspath = self.inputs['subpaths']['spots']
        picklespath = self.inputs['subpaths']['ccdpickles']


        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        spotsOUT_cdp = self.CDP_lib['SPOTS_NOBFE']

        spotsOUT_cdp.rootname = '{}_{}_{}'.format(spotsOUT_cdp.rootname,
                self.inputs['test'],
                self.inputs['BLOCKID'])

        spotsOUT_cdp.path = spotspath

        spotsOUT_cdp.header = CDP_header.copy()
        spotsOUT_cdp.meta = dict()
        spotsOUT_cdp.data = OrderedDict(spots=np.zeros((nObs, nCCDs, nQuads, nSpots), dtype=object))


        if not self.drill:
            
            spotsIN_array = files.cPickleRead(self.dd.products['SPOTS'])['data']['spots'].copy()

            for iObs in range(nObs):

                for jCCD, CCDk in enumerate(CCDs):

                    for kQ, Quad in enumerate(Quads):

                        _Asol = Asols[CCDk][Quad]

                        for lS, SpotName in enumerate(SpotNames):


                            inSpot = copy.deepcopy(spotsIN_array[iObs, jCCD, kQ, lS])
                            img = inSpot.data.copy()
                            try:
                                oimg = G15.correct_estatic(img, _Asol)
                                outSpot = copy.deepcopy(inSpot)
                                outSpot.data = oimg.copy()
                                spotsOUT_cdp.data['spots'][iObs, jCCD, kQ, lS] = copy.deepcopy(outSpot)
                            except:
                                spotsOUT_cdp.data['spots'][iObs, jCCD, kQ, lS] = None

        spotsOUT_cdp.meta.update(dict(nObs=nObs,
            nCCDs=nCCDs,
            nQuads=nQuads,
            nSpots=nSpots,
            SpotNames=SpotNames,
            structure='spots:np.zeros((nObs, nCCDs, nQuads, nSpots), dtype=object)'))

        self.save_CDP(spotsOUT_cdp)
        self.pack_CDP_to_dd(spotsOUT_cdp, 'SPOTS_NOBFE')

        if self.log is not None:
            self.log.info('Saved BFE-corrected spot "bag" file to %s' % spotspath)

    def _extract_basic(self, spotsbag='SPOTS', prefix=''):
        """ """

        CQSindices = copy.deepcopy(self.dd.indices)
        Npxres = 10 # Npx of of residal image of gaussian fit (on a side)

        assert 'Spot' in CQSindices.names

        colnames = ['bas_bgd', 'bas_peak', 'bas_fluence', 'bas_efluence',
            'bas_x', 'bas_y', 'bas_x_ccd', 'bas_y_ccd',
            'bas_fwhmx', 'bas_fwhmy',
            'sh_x', 'sh_y', 'sh_e1', 'sh_e2',
            'sh_ell', 'sh_R2', 'sh_R2arcsec',
            'sh_a', 'sh_b', 
            'gau_bgd', 'gau_ebgd',
            'gau_i00', 'gau_ei00',
            'gau_xcen', 'gau_excen',
            'gau_ycen', 'gau_eycen',
            'gau_sigmax', 'gau_esigmax',
            'gau_fwhmx', 'gau_efwhmx',
            'gau_sigmay', 'gau_esigmay',
            'gau_fwhmy', 'gau_efwhmy', 
            'gau_didfit',
            'gau_res_xskew', 'gau_res_yskew']

        ncols = len(colnames)
        for i in range(ncols):
            colnames[i] = '%s%s' % (prefix, colnames[i])

        for i in range(ncols):
            colname = colnames[i]
            if 'didfit' in colname:
                idtype = 'int32'
                valini = 0
            else:
                idtype = 'float32'
                valini = 0.0
            self.dd.initColumn(colname, CQSindices, dtype=idtype,
                 valini=valini)
            

        nObs = CQSindices[0].len
        CCDs = CQSindices.get_vals('CCD')
        Quads = CQSindices.get_vals('Quad')
        SpotNames = CQSindices.get_vals('Spot')

        if self.drill:
            return

        spots_array = files.cPickleRead(self.dd.products[spotsbag])['data']['spots'].copy()

        #psCCDcoodicts = self._get_psCCDcoodicts()

        bas_keys = ['bgd','peak','fluence','efluence',
                        'x','y','x_ccd','y_ccd','fwhmx','fwhmy']
        bas_res_def = dict(zip(bas_keys, np.zeros(
                                    len(bas_keys), dtype='float32')))

        for iObs in range(nObs):

            for jCCD, CCDk in enumerate(CCDs):

                for kQ, Quad in enumerate(Quads):

                    for lS, SpotName in enumerate(SpotNames):

                        ixtup = (iObs, jCCD, kQ, lS)

                        inSpot = copy.deepcopy(spots_array[ixtup])

                        try: 
                            bas_res = inSpot.measure_basic(rap=10, rin=15, rout=-1,
                            gain=3.5, debug=False)
                        except BaseException:
                            continue
                        # bas_res = dict(bgd=bgd, peak=peak, fluence=flu, efluence=eflu,
                        #    x=x, y=y, x_ccd=x_ccd, y_ccd=y_ccd,
                        #    fwhmx=fwhmx, fwhmy=fwhmy)

                        self.dd.mx['%sbas_bgd' % prefix][ixtup] = bas_res['bgd']
                        self.dd.mx['%sbas_peak' % prefix][ixtup] = bas_res['peak']
                        self.dd.mx['%sbas_fluence' % prefix][ixtup] = bas_res['fluence']
                        self.dd.mx['%sbas_efluence' % prefix][ixtup] = bas_res['efluence']
                        self.dd.mx['%sbas_x' % prefix][ixtup] = bas_res['x']
                        self.dd.mx['%sbas_y' % prefix][ixtup] = bas_res['y']
                        self.dd.mx['%sbas_x_ccd' % prefix][ixtup] = bas_res['x_ccd']
                        self.dd.mx['%sbas_y_ccd' % prefix][ixtup] = bas_res['y_ccd']
                        self.dd.mx['%sbas_fwhmx' % prefix][ixtup] = bas_res['fwhmx']
                        self.dd.mx['%sbas_fwhmy' % prefix][ixtup] = bas_res['fwhmy']


                        #inSpot.data -= loc_bgd # UNNECESSARY? background subtraction

                        inSpot.shsettings = dict(iterations=4,
                               sampling=1.0,  # oversampling if > 1
                               platescale = 120.0,  # microns / arcsecond
                               pixelSize=12.0,  # um/pix
                               sigma=0.75,  # arcseconds
                               weighted=True,  # use weighted moments
                               conservePeak=True,
                               debug=False,
                               fixedPosition=True,
                               fixedX=bas_res['x'],
                               fixedY=bas_res['y'])

                        ref_quad_res = inSpot.measureRefinedEllipticity()

                        #ref_quad_res = dict(centreX=quad['centreX'] + 1, 
                        #    centreY=quad['centreY'] + 1,
                        #    e1=quad['e1'], e2=quad['e2'],
                        #    ellipticity=quad['ellipticity'],
                        #    R2=R2,
                        #    R2arcsec=R2arcsec,
                        #    GaussianWeighted=GaussianWeighted,
                        #    a=quad['a'], b=quad['b'])

                        self.dd.mx['%ssh_x' % prefix][ixtup] = ref_quad_res['centreX']
                        self.dd.mx['%ssh_y' % prefix][ixtup] = ref_quad_res['centreY']
                        self.dd.mx['%ssh_e1' % prefix][ixtup] = ref_quad_res['e1']
                        self.dd.mx['%ssh_e2' % prefix][ixtup] = ref_quad_res['e2']
                        self.dd.mx['%ssh_ell' % prefix][ixtup] = ref_quad_res['ellipticity']
                        self.dd.mx['%ssh_R2' % prefix][ixtup] = ref_quad_res['R2']
                        self.dd.mx['%ssh_R2arcsec' % prefix][ixtup] = ref_quad_res['R2arcsec']
                        self.dd.mx['%ssh_a' % prefix][ixtup] = ref_quad_res['a']
                        self.dd.mx['%ssh_b' % prefix][ixtup] = ref_quad_res['b']


                        #gauss_res = inSpot.fit_Gauss()
                        # tuple: (vals, evals)
                        # i00, xcen, ycen, xsigma, ysigma

                        gauss_res = inSpot.get_shape_Gauss()

                        self.dd.mx['%sgau_bgd' % prefix][ixtup] = gauss_res['bgd']
                        self.dd.mx['%sgau_ebgd' % prefix][ixtup] = gauss_res['ebgd']
                        self.dd.mx['%sgau_i00' % prefix][ixtup] = gauss_res['i0']
                        self.dd.mx['%sgau_ei00' % prefix][ixtup] = gauss_res['ei0']
                        self.dd.mx['%sgau_xcen' % prefix][ixtup] = gauss_res['x']
                        self.dd.mx['%sgau_xcen' % prefix][ixtup] = gauss_res['ex']
                        self.dd.mx['%sgau_ycen' % prefix][ixtup] = gauss_res['y']
                        self.dd.mx['%sgau_eycen' % prefix][ixtup] = gauss_res['ey']
                        self.dd.mx['%sgau_sigmax' % prefix][ixtup] = gauss_res['sigma_x']
                        self.dd.mx['%sgau_esigmax' % prefix][ixtup] = gauss_res['esigma_x']
                        self.dd.mx['%sgau_fwhmx' % prefix][ixtup] = gauss_res['fwhm_x']
                        self.dd.mx['%sgau_efwhmx' % prefix][ixtup] = gauss_res['efwhm_x']
                        self.dd.mx['%sgau_sigmay' % prefix][ixtup] = gauss_res['sigma_y']
                        self.dd.mx['%sgau_esigmay' % prefix][ixtup] = gauss_res['esigma_y']
                        self.dd.mx['%sgau_fwhmy' % prefix][ixtup] = gauss_res['fwhm_y']
                        self.dd.mx['%sgau_efwhmy' % prefix][ixtup] = gauss_res['efwhm_y']
                        self.dd.mx['%sgau_didfit' % prefix][ixtup] = gauss_res['didFit']

                        if gauss_res['didFit']:

                            ibgd = gauss_res['bgd']
                            i0 = gauss_res['i0']
                            ixcen = gauss_res['x']
                            iycen = gauss_res['y']
                            ix_stddev = gauss_res['sigma_x']
                            iy_stddev = gauss_res['sigma_y']

                            gaussmod = models.Const2D(ibgd) + models.Gaussian2D(
                                i0, ixcen, iycen, ix_stddev, iy_stddev, theta=0.)

                            XX, YY = np.meshgrid(np.arange(inSpot.NX),
                                 np.arange(inSpot.NY), indexing='xy')

                            gaussEval = gaussmod(XX,YY)
            
                            iresiduals = inSpot.data - gaussEval

                            ipos = (ixcen, iycen)
                            isize = (Npxres, Npxres)     # pixels
                            irescutout = Cutout2D(iresiduals, ipos, isize).data.copy()
                            irescutout /= gauss_res['fluence'] # normalization

                            self.dd.mx['%sgau_res_xskew' % prefix][ixtup] = \
                                stats.skew(np.mean(irescutout,axis=1))
                            self.dd.mx['%sgau_res_yskew' % prefix][ixtup] = \
                                stats.skew(np.mean(irescutout,axis=0))


    def basic_analysis(self):
        """Performs basic analysis on spots:
             - shape from moments
             - Gaussian fitting: peak intensity, position, width_x, width_y

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='basic', Title='Basic Extraction of Spots Shapes', level=0)

        self._extract_basic(spotsbag='SPOTS', prefix='')

    def basic_analysis_noBFE(self):
        """Performs basic analysis on [BFE-corrected] spots:
             - shape from moments
             - Gaussian fitting: peak intensity, position, width_x, width_y

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='nobfe_basic', Title='BFE corrected: Basic Extraction of Spots Shapes', level=0)

        self._extract_basic(spotsbag='SPOTS_NOBFE', prefix='nobfe_')        

    def _build_SpotsPoster(self, spotsbag='SPOTS', histequ=True):
        """ """

        CQSindices = copy.deepcopy(self.dd.indices)

        assert 'Spot' in CQSindices.names

        nObs = CQSindices[0].len
        CCDs = CQSindices.get_vals('CCD')
        nCCDs = len(CCDs)
        Quads = CQSindices.get_vals('Quad')
        nQuads = len(Quads)
        SpotNames = CQSindices.get_vals('Spot')
        nSpots = len(SpotNames)

        if self.drill:
            return

        spots_array = files.cPickleRead(self.dd.products[spotsbag])['data']['spots'].copy()

        psCCDcoodicts = self._get_psCCDcoodicts()

        NY = stampw * nObs
        NX = stampw * nCCDs * nQuads * nSpots

        SpotsPoster = np.zeros((NY,NX),dtype='float32')


        for iObs in range(nObs):

            for jCCD, CCDk in enumerate(CCDs):


                for kQ, Quad in enumerate(Quads):

                    for lS, SpotName in enumerate(SpotNames):

                        ixtup = (iObs, jCCD, kQ, lS)

                        inSpot = copy.deepcopy(spots_array[ixtup])

                        lly = stampw * iObs
                        llx = stampw * (jCCD*(nQuads+nSpots)+kQ*(nSpots)+lS)

                        SpotsPoster[lly:lly+stampw, llx:llx+stampw] = inSpot.data.copy()

        SpotsPoster[SpotsPoster<=0.] = 1.
        SpotsPoster = np.log10(SpotsPoster)

        if histequ:
            SpotsPoster = exposure.equalize_hist(SpotsPoster, nbins=256)

        return SpotsPoster

    def bayes_analysis(self):
        """
        Performs bayesian decomposition of the spot images:
            - optomechanic PSF and detector PSF.
        Also measures R2, fwhm and ellipticity of "extracted" detector PSF.
        """

        # TESTS

        from pylab import plot, show
        from scipy.optimize import curve_fit

        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices.get_vals('CCD')
        #nCCD = len(CCDs)
        Quads = dIndices.get_vals('Quad')
        nQuads = len(Quads)
        nObs = dIndices[dIndices.names.index('ix')].len
        SpotNames = dIndices.get_vals('Spot')
        nSpots = len(SpotNames)

        jCCD = 0
        lS = 0
        iQ = Quads.index('E')

        fwhmx = self.dd.mx['chk_fwhmx'][:, jCCD, iQ, lS].copy()
        fwhmy = self.dd.mx['chk_fwhmy'][:, jCCD, iQ, lS].copy()
        peak = self.dd.mx['chk_peak'][:, jCCD, iQ, lS].copy()

        def f_gauss_ccd_fit(peakflu, *p):
            hwc = 2.**15.
            return np.sqrt(p[0]**2. + (p[1] * peakflu / hwc)**2.)

        selx = np.where((peak < 5.e4) & (fwhmx > 1.5) & (fwhmx < 2.5))
        fitcoefsx, pcov = curve_fit(f_gauss_ccd_fit, peak[selx], fwhmx[selx], p0=[2.0, 0.1])

        sely = np.where((peak < 5.e4) & (fwhmy > 1.5) & (fwhmy < 2.5))
        fitcoefsy, pcov = curve_fit(f_gauss_ccd_fit, peak[sely], fwhmy[sely], p0=[2.0, 0.1])

        spotspath = self.inputs['subpaths']['spots']

        spotslist = []

        for iObs in range(nObs):

            spotspick = os.path.join(spotspath, '%s.pick' % self.dd.mx['spots_name'][iObs, jCCD])
            spotsobj = cPickleRead(spotspick)
            spotslist.append(spotsobj['spots'][iQ, lS].data.copy())

        stop()

        for iObs in range(nObs):

            img = spotslist[iObs]
            centroid = np.unravel_index(img.argmax(), img.shape)
            plot(img[centroid[0], :] / img.max())
        show()

        stop()

    def _get_fwhm_flu_bfit(self, iCCD, kQ, fwhmkey, bfecorr):
        """ """

        xkey = 'gau_i00'
        ykey = 'gau_%s' % fwhmkey
        fitkey = 'gau_didfit'
        if bfecorr:
            xkey = 'nobfe_%s' % xkey
            ykey = 'nobfe_%s' % ykey
            fitkey = 'nobfe_%s' % fitkey

        xdata = self.dd.mx[xkey][:,iCCD,kQ,:].copy()
        ydata = self.dd.mx[ykey][:,iCCD,kQ,:].copy()
        didfit = self.dd.mx[fitkey][:,iCCD,kQ,:].copy()

        nSpots = xdata.shape[1]

        slopes = []
        intercepts = []

        for i in range(nSpots):
            _x = xdata[:,i]
            _y = ydata[:,i]
            _df = didfit[:,i]

            ixsel = np.where((_x>1.e3) & (_df==1))

            if len(ixsel[0])>5:
                ransac = linear_model.RANSACRegressor()

                x2fit = (_x[ixsel]-2.**15)/1.e4

                x2fit = np.expand_dims(x2fit, 1)

                ransac.fit(x2fit, np.expand_dims(_y[ixsel], 1))

                slopes.append(ransac.estimator_.coef_[0][0])
                intercepts.append(ransac.estimator_.intercept_[0])

        slope = np.nanmean(slopes)
        intercept = np.nanmean(intercepts)
        coeffs = [slope,intercept]

        p = np.poly1d(coeffs)
        xfit = (np.linspace(1.,2.**16,5)-2.**15)/1.e4
        ybest = np.polyval(p,xfit)

        xbest = xfit+2.**15/1.e4

        return xbest, ybest, coeffs


    def meta_analysis(self):
        """

        Analyzes the relation between detector PSF and fluence, etc.

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='Meta Analysis: Shapes', level=0)

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        doBuildPosters = False

        if doBuildPosters:

            SpotsPoster = self._build_SpotsPoster(spotscol='SPOTS', 
                histequ=False)
            fdict = self.figdict['SpotsPoster'][1]
            fdict['data'] = SpotsPoster.copy()

            SpotsPosterNOBFE = self._build_SpotsPoster(spotscol='SPOTS_NOBFE',
                histequ=False)
            fdictNOBFE = self.figdict['SpotsPosterNOBFE'][1]
            fdictNOBFE['data'] = SpotsPosterNOBFE.copy()        


            #normfunction = Normalize(vmin=0.,vmax=np.log10(2.**16),clip=False)
            #fdict['meta']['corekwargs']['norm'] = normfunction

            if self.report is not None:
                self.addFigures_ST(figkeys=['SpotsPoster',
                    'SpotsPosterNOBFE'],
                    dobuilddata=False)

        

        # INITIALISATIONS

        indices = copy.deepcopy(self.dd.indices)
        nObs, nC, nQ, nSpots = indices.shape
        CCDs = np.array(indices.get_vals('CCD'))
        Quads = np.array(indices.get_vals('Quad'))
        SpotNames = np.array(indices.get_vals('Spot'))


        # PLOT OF FWHMj vs. peak fluence (w/w-o BFE corr)

        plot_FWHM_dict = OrderedDict()

        for tag in ['fwhmx', 'fwhmy']:
            plot_FWHM_dict[tag] = OrderedDict(labelkeys=['BFE', 'noBFE', 'ideal'])
            for CCDk in CCDs:
                plot_FWHM_dict[tag][CCDk] = OrderedDict()
                for Q in Quads:
                    plot_FWHM_dict[tag][CCDk][Q] = OrderedDict()
                    plot_FWHM_dict[tag][CCDk][Q]['x'] = OrderedDict()
                    plot_FWHM_dict[tag][CCDk][Q]['y'] = OrderedDict()



        for iCCD, CCDk in enumerate(CCDs):

            for kQ, Q in enumerate(Quads):

                for bfecorr in [True, False]:

                    if bfecorr:
                        BFEtag = 'noBFE'
                    else:
                        BFEtag = 'BFE'

                    x_fwhmx,y_fwhmx, px = self._get_fwhm_flu_bfit(iCCD, kQ, 
                            'fwhmx',bfecorr=bfecorr)


                    plot_FWHM_dict['fwhmx'][CCDk][Q]['x'][BFEtag] = x_fwhmx
                    plot_FWHM_dict['fwhmx'][CCDk][Q]['y'][BFEtag] = y_fwhmx

                    if bfecorr:
                        xideal = np.array([1.,2.**16])
                        pideal = np.poly1d([0.,px[1]])
                        yideal = np.polyval(pideal,xideal)
                        plot_FWHM_dict['fwhmx'][CCDk][Q]['y']['ideal'] = yideal
                        plot_FWHM_dict['fwhmx'][CCDk][Q]['x']['ideal'] = xideal-2.**15/1.e4

                    x_fwhmy,y_fwhmy, py = self._get_fwhm_flu_bfit(iCCD, kQ, 
                            'fwhmy',bfecorr=bfecorr)

                    plot_FWHM_dict['fwhmy'][CCDk][Q]['x'][BFEtag] = x_fwhmy 
                    plot_FWHM_dict['fwhmy'][CCDk][Q]['y'][BFEtag] = y_fwhmy

                    if bfecorr:
                        xideal = np.array([1.,2.**16])
                        pideal = np.poly1d([0.,py[1]])
                        yideal = np.polyval(pideal,xideal)
                        plot_FWHM_dict['fwhmy'][CCDk][Q]['y']['ideal'] = yideal
                        plot_FWHM_dict['fwhmy'][CCDk][Q]['x']['ideal'] = xideal-2.**15/1.e4


        for tag in ['fwhmx', 'fwhmy']:

            fdict_fwhm = self.figdict['PSF0X_%s_v_flu' % tag][1]
            fdict_fwhm['data'] = plot_FWHM_dict[tag].copy()

            if self.report is not None:
                self.addFigures_ST(figkeys=['PSF0X_%s_v_flu' % tag],
                                   dobuilddata=False)

        # Plot of skweness of gauss-fit residuals vs. peak fluence and x/y directions 
        plot_skew_dict = OrderedDict()


        for tag1 in ['vs_pos', 'vs_fluence']:
            plot_skew_dict[tag1] = OrderedDict()
            for tag2 in ['dirx', 'diry']:
                plot_skew_dict[tag1][tag2] = OrderedDict(labelkeys=['BFE', 'noBFE'])
                for CCDk in CCDs:
                    plot_skew_dict[tag1][tag2][CCDk] = OrderedDict()
                    for Q in Quads:
                        plot_skew_dict[tag1][tag2][CCDk][Q] = OrderedDict()
                        plot_skew_dict[tag1][tag2][CCDk][Q]['x'] = OrderedDict()
                        plot_skew_dict[tag1][tag2][CCDk][Q]['y'] = OrderedDict()
                        plot_skew_dict[tag1][tag2][CCDk][Q]['ey'] = OrderedDict()


        for iCCD, CCDk in enumerate(CCDs):

            for kQ, Q in enumerate(Quads):

                for bfecorr in [True, False]:

                    if bfecorr:
                        BFEtag = 'noBFE'
                        colbfetag = 'nobfe_'
                    else:
                        BFEtag = 'BFE'
                        colbfetag = ''


                    # x-skew vs. x-ccd

                    _x_vpos_kx = self.dd.mx['chk_x_ccd'][:,iCCD, kQ,:].flatten()

                    if Q in ['F','G']:
                        _x_vpos_kx = self.ccdcalc.NAXIS1 - _x_vpos_kx

                    _y_vpos_kx = self.dd.mx['{}gau_res_xskew'.format(colbfetag)][:,iCCD, kQ,:].flatten()

                    bin_vpos_kx = PSF0Xaux._f_xy_bin(_x_vpos_kx,_y_vpos_kx,Nbins=3)
                    x_vpos_kx, y_vpos_kx, sigy_vpos_kx, n_vpos_kx = bin_vpos_kx

                    plot_skew_dict['vs_pos']['dirx'][CCDk][Q]['x'][BFEtag] = x_vpos_kx  
                    plot_skew_dict['vs_pos']['dirx'][CCDk][Q]['y'][BFEtag] = y_vpos_kx
                    plot_skew_dict['vs_pos']['dirx'][CCDk][Q]['ey'][BFEtag] = sigy_vpos_kx / n_vpos_kx**0.5

                    # y-skew vs. y-ccd

                    _x_vpos_ky = self.dd.mx['chk_y_ccd'][:,iCCD, kQ,:].flatten()

                    stop()

                    if Q in ['E','F']:
                        _x_vpos_ky = self.ccdcalc.NAXIS2 - _x_vpos_ky


                    _y_vpos_ky = self.dd.mx['{}gau_res_yskew'.format(colbfetag)][:,iCCD, kQ,:].flatten()

                    bin_vpos_ky = PSF0Xaux._f_xy_bin(_x_vpos_ky,_y_vpos_ky,Nbins=3)
                    x_vpos_ky, y_vpos_ky, sigy_vpos_ky, n_vpos_ky = bin_vpos_ky

                    plot_skew_dict['vs_pos']['diry'][CCDk][Q]['x'][BFEtag] = x_vpos_ky  
                    plot_skew_dict['vs_pos']['diry'][CCDk][Q]['y'][BFEtag] = y_vpos_ky
                    plot_skew_dict['vs_pos']['diry'][CCDk][Q]['ey'][BFEtag] = sigy_vpos_ky / n_vpos_ky**0.5

                    # x-skew vs. fluence

                    _x_vflu_kx = self.dd.mx['bas_fluence'][:,iCCD, kQ,:].flatten() / 1.e3
                    _y_vflu_kx = self.dd.mx['{}gau_res_xskew'.format(colbfetag)][:,iCCD, kQ,:].flatten()

                    Nbinsflu = len(self.inputs['exptimes'])

                    bin_vflu_kx = PSF0Xaux._f_xy_bin(_x_vflu_kx,_y_vflu_kx,Nbins=Nbinsflu)
                    x_vflu_kx, y_vflu_kx, sigy_vflu_kx, n_vflu_kx = bin_vflu_kx

                    plot_skew_dict['vs_fluence']['dirx'][CCDk][Q]['x'][BFEtag] = x_vflu_kx
                    plot_skew_dict['vs_fluence']['dirx'][CCDk][Q]['y'][BFEtag] = y_vflu_kx
                    plot_skew_dict['vs_fluence']['dirx'][CCDk][Q]['ey'][BFEtag] = sigy_vflu_kx / n_vflu_kx**0.5

                    # y-skew vs. fluence

                    _x_vflu_ky = self.dd.mx['bas_fluence'][:,iCCD, kQ,:].flatten() / 1.e3
                    _y_vflu_ky = self.dd.mx['{}gau_res_yskew'.format(colbfetag)][:,iCCD, kQ,:].flatten()

                    bin_vflu_ky = PSF0Xaux._f_xy_bin(_x_vflu_ky,_y_vflu_ky,Nbins=Nbinsflu)
                    x_vflu_ky, y_vflu_ky, sigy_vflu_ky, n_vflu_ky = bin_vflu_ky

                    plot_skew_dict['vs_fluence']['diry'][CCDk][Q]['x'][BFEtag] = x_vflu_ky
                    plot_skew_dict['vs_fluence']['diry'][CCDk][Q]['y'][BFEtag] = y_vflu_ky
                    plot_skew_dict['vs_fluence']['diry'][CCDk][Q]['ey'][BFEtag] = sigy_vflu_ky / n_vflu_ky**0.5


        for tag1 in ['vs_pos', 'vs_fluence']:
            for tag2 in ['dirx', 'diry']:

                figtag = 'PSF0X_skew_%s_%s' % (tag2,tag1)

                fdict_skew = self.figdict[figtag][1]
                fdict_skew['data'] = plot_skew_dict[tag1][tag2].copy()

                if self.report is not None:
                    self.addFigures_ST(figkeys=[figtag],
                                   dobuilddata=False)


    def opt_xtalk_sextract(self):
        """Runs sextractor on images for optical-crosstalk measurements."""

        if self.report is not None:
            self.report.add_Section(
                keyword='xtalksex', Title='Cross-Talk SExtraction', level=0)

        DDindices = copy.deepcopy(self.dd.indices)
        nObs = DDindices.get_len('ix')
        CCDs = DDindices.get_vals('CCD')

        if not self.drill:

            for iObs in range(nObs):

                print('\nopt_xtalk: sextracting ObsID %i/%i' % (iObs + 1, nObs))

                ObsID = self.dd.mx['ObsID'][iObs]
                timestamp = self.dd.mx['date'][iObs, 0]

                iobs_pick = os.path.join(self.inputs['subpaths']['xtalk'],
                                         'EUC_%i_%s_ROE1_sex.pick' % (ObsID, timestamp))

                iobs_sexdata = OrderedDict()

                for jCCD, CCDkey in enumerate(CCDs):

                    dpath = self.dd.mx['datapath'][iObs, jCCD]

                    infits = os.path.join(dpath, '%s.fits' %
                                          self.dd.mx['File_name'][iObs, jCCD])

                    ccdobj = ccd.CCD(infits)

                    sextag = 'EUC_%i_%s_%s_sex' % (ObsID, timestamp, CCDkey)
                    iresSEx = oxt.exe_SEx(ccdobj, tag=sextag)

                    iobs_sexdata[CCDkey] = copy.deepcopy(iresSEx)

                cPickleDumpDictionary(iobs_sexdata, iobs_pick)

            if self.report is not None:
                self.report.add_Text('%i ObsIDs SExtracted!' % nObs)

    def opt_xtalk_build(self):
        """ """

        if self.report is not None:
            self.report.add_Section(
                keyword='xtalkbuild', Title='Cross-Talk Matrix Building', level=0)

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        xtalkpath = self.inputs['subpaths']['xtalk']

        if not self.drill:

            rawcrosstalks = oxt.get_rawcrosstalk_mx(self.dd, xtalkpath)

            rawct_cdp = self.CDP_lib['RAW_CTALK']

            rawct_cdp.rootname = '%s_%s_%s' % (rawct_cdp.rootname,
                                               self.inputs['test'],
                                               self.inputs['BLOCKID'])
            rawct_cdp.path = self.inputs['subpaths']['xtalk']

            rawct_cdp.header = CDP_header.copy()
            rawct_cdp.meta = dict()
            rawct_cdp.data = rawcrosstalks.copy()

            self.save_CDP(rawct_cdp)
            self.pack_CDP_to_dd(rawct_cdp, 'RAW_CTALK')

            if self.report is not None:
                self.report.add_Text('Raw Crosstalk Ramps extracted.')

    def opt_xtalk_meta(self):
        """ """

        if self.report is not None:
            self.report.add_Section(
                keyword='xtalkmeta', Title='Cross-Talk Matrix Meta-Analysis', level=0)

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        Quads = self.ccdcalc.Quads

        figspath = self.inputs['subpaths']['figs']

        if not self.drill:

            raw_ctalk_pick = self.dd.products['RAW_CTALK']
            raw_ctalk = cPickleRead(raw_ctalk_pick)['data'].copy()

            crosstalks = oxt.get_crosstalk_mx(raw_ctalk, CCDs, Quads)

            ct_cdp = self.CDP_lib['CTALK']

            ct_cdp.rootname = '%s_%s_%s' % (ct_cdp.rootname,
                                            self.inputs['test'],
                                            self.inputs['BLOCKID'])
            ct_cdp.path = self.inputs['subpaths']['xtalk']

            ct_cdp.header = CDP_header.copy()
            ct_cdp.meta = dict()
            ct_cdp.data = crosstalks.copy()

            self.save_CDP(ct_cdp)
            self.pack_CDP_to_dd(ct_cdp, 'CTALK')

            fignameratio = self.figdict['PSF0X_crosstalk_RATIO'][1]['figname']
            fullfignameratio = os.path.join(figspath, fignameratio)

            xt.PlotSummaryFig(crosstalks, suptitle='', figname=fullfignameratio,
                              scale='RATIO', showvalues=True)

            fignameadu = self.figdict['PSF0X_crosstalk_ADU'][1]['figname']
            fullfignameadu = os.path.join(figspath, fignameadu)

            xt.PlotSummaryFig(crosstalks, suptitle='', figname=fullfignameadu,
                              scale='ADU', showvalues=True)

            if self.report is not None:
                crosstalkfigs = ['PSF0X_crosstalk_RATIO',
                                 'PSF0X_crosstalk_ADU']
                self.addFigures_ST(figkeys=crosstalkfigs, dobuilddata=False)

        self.canbecleaned = True 
        
