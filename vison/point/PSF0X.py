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

from vison.support import utils
from vison.pipe.task import HKKeys
from vison.support import context
from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
from vison.datamodel import scriptic as sc
from vison.support import files
from vison.point import PointTask as PT
import PSF0Xaux
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
                        shuttr=1, e_shuttr=0,
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
                         ('corrBFE_G15',self.corrBFE_G15),
                         ('basic', self.basic_analysis),
                         ('bayes', self.bayes_analysis),
                         ('meta', self.meta_analysis),
                         ('xtalk_sex', self.opt_xtalk_sextract),
                         ('xtalk_build', self.opt_xtalk_build),
                         ('xtalk_meta', self.opt_xtalk_meta)]

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

        Ncols = len(PSF0X_sdict.keys())
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
        CCDs = self.ogse.startrackers.keys()
        ulabels = self.ogse.startrackers[CCDs[0]].keys()

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

                    spots_array = np.zeros((nQuads, nSpots), dtype=object)

                    for kQ, Quad in enumerate(Quads):

                        for lS, SpotName in enumerate(SpotNames):

                            #coo = polib.Point_CooNom[CCDk][Quad][SpotName]
                            coo = psCCDcoodicts[CCDk][ilabel][Quad][SpotName]
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

                    if onTests:
                    #doPlot=True
                    #if doPlot and iObs==3 and CCDk=='CCD1':

                        from matplotlib import pyplot as plt

                        for kQ, Quad in enumerate(Quads):

                            for lS, SpotName in enumerate(SpotNames):

                                img = spotsdict['spots'][kQ, lS].data
                                fig = plt.figure()
                                ax = fig.add_subplot(111)
                                ax.imshow(img.T, origin='lower left')
                                ax.set_title('Obs%i, %s-%s, %s' % \
                                    (iObs+1,CCDk,Quad,SpotName))
                                plt.show()


        if self.log is not None:
            self.log.info('Saved spot "bag" files to %s' % spotspath)


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
        midrangekey = ulabels[len(ulabels)/2]
        inCCDs = BFEall['CCDs']
        inQuads = BFEall['Quads']

        Asols = dict()
        for inCCD in inCCDs:
            Asols[inCCD] = dict()
            for Q in inQuads:
                Asols[inCCD][Q] = BFEall[inCCD][Q][midrangekey]['Asol'].copy()


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

        outcol = 'spots_name_nobfe'
        self.dd.initColumn(outcol, CIndices, dtype='S100', valini='None')

        if not self.drill:

            spotspath = self.inputs['subpaths']['spots']
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
            #for iObs in range(3): # TESTS
            #for iObs in [20]: # TESTs

                for jCCD, CCDk in enumerate(CCDs):

                    fullINspots_name = os.path.join(
                        spotspath, '%s.pick' % self.dd.mx['spots_name'][iObs, jCCD])

                    spotsIN_array = files.cPickleRead(fullINspots_name)['spots']
                    

                    # Process individual "spots"

                    spotsOUT_array = np.zeros((nQuads, nSpots), dtype=object)

                    for kQ, Quad in enumerate(Quads):

                        _Asol = Asols[CCDk][Quad]

                        for lS, SpotName in enumerate(SpotNames):


                            inSpot = copy.deepcopy(spotsIN_array[kQ, lS])
                            img = inSpot.data.copy()
                            try:
                                oimg = G15.correct_estatic(img, _Asol)
                                outSpot = copy.deepcopy(inSpot)
                                outSpot.data = oimg.copy()
                                spotsOUT_array[kQ, lS] = copy.deepcopy(outSpot)
                            except:
                                spotsOUT_array[kQ, lS] = None

                    # save "spots" to a separate file and keep name in dd

                    self.dd.mx[outcol][iObs,
                                             jCCD] = '%s_spots_nobfe' % self.dd.mx['File_name'][iObs, jCCD]

                    fullspots_name = os.path.join(
                        spotspath, '%s.pick' % self.dd.mx[outcol][iObs, jCCD])

                    spotsdict = dict(spots=spotsOUT_array)

                    files.cPickleDumpDictionary(spotsdict, fullspots_name)

        if self.log is not None:
            self.log.info('Saved BFE-corrected spot "bag" files to %s' % spotspath)

    def basic_analysis(self):
        """Performs basic analysis on spots:
             - shape from moments
             - Gaussian fitting: peak intensity, position, width_x, width_y

        """

        raise NotImplementedError

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

        CCDs = dIndices[dIndices.names.index('CCD')].vals
        #nCCD = len(CCDs)
        Quads = dIndices[dIndices.names.index('Quad')].vals
        nQuads = len(Quads)
        nObs = dIndices[dIndices.names.index('ix')].len
        SpotNames = dIndices[dIndices.names.index('Spot')].vals
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

    # def meta_analysis(dd,report,inputs,log=None):
    def meta_analysis(self):
        """

        Analyzes the relation between detector PSF and fluence.

        """
        raise NotImplementedError

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
