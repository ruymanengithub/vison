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

from vison.pipe.task import HKKeys
from vison.support import context
from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
from vison.datamodel import scriptic as sc
from vison.support import files
from vison.point import PointTask as PT
import PSF0Xaux
from vison.support.files import cPickleRead
from vison.datamodel import ccd
# END IMPORT

isthere = os.path.exists

stampw = polib.stampw

#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
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
                        mirr_pos = ogsemod.mirror_nom['F4'],
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
            (PSF0X_relfluences/100.*tFWC_pointw).tolist()

    return testdefaults


def get_Flu_lims(relfluences):

    plusminus10pc = 1.+np.array([-0.1, 0.1])
    satur_fluence = 2.**16
    
    # assuming a best fwhm~2pix, and a gaussian profile
    # F = 2pi(fwhm/2.355)**2*I0
    
    I02F = 2*np.pi*(2./2.355)**2.

    Flu_lims = OrderedDict(
        CCD1=OrderedDict(
            E=OrderedDict(
                ALPHA=OrderedDict())))  # +/-10%

    Nfluences = len(relfluences)
    for i in range(1, Nfluences+1):
        relflu = min( relfluences[i-1]/100.,1)
        
        Flu_lims['CCD1']['E']['ALPHA']['col%03i' % i] = \
            I02F * relflu * plusminus10pc*satur_fluence

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


class PSF0X_inputs(PT.Point_inputs):
    manifesto = PT.Point_inputs.manifesto.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
    ])))


class PSF0X(PT.PointTask):

    inputsclass = PSF0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('lock', self.lock_on_stars),
                         ('check', self.check_data), 
                         ('prep', self.prep_data),
                         ('basic', self.basic_analysis), 
                         ('bayes', self.bayes_analysis),
                         ('meta', self.meta_analysis),
                         ('xtalk', self.opt_xtalk)]
        super(PSF0X, self).__init__(inputs, log, drill, debug)
        self.name = 'PSF0X'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.CDP_lib = PSF0Xaux.get_CDP_lib(self.inputs['test'])
        self.figdict = PSF0Xaux.get_PSF0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                   products='products', spots='spots',
                   xtalk='xtalk')
        
        

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
                                offsetxy=[0.,0.],
                                wavelength=wavelength,
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
            colid = 'col%03i' % (ic+1,)
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
        
    def lock_on_stars(self):
        """ """
        
        FW_ID = self.dd.mx['wave'][0,0]
        wave = self.ogse.get_wavelength(FW_ID)
        tFWC_point = self.ogse.profile['tFWC_point']['nm%i' % wave]
        exptime = self.dd.mx['exptime'][:,0]
        
        iObs = np.abs(exptime-tFWC_point/2.).argmin()
        
        
        PT.PointTask.lock_on_stars(self,iObs=iObs,
                                   sexconfig=dict(
                                           MINAREA=3.,
                                           DET_THRESH=15.,
                                           MAG_ZEROPOINT=20.)) 

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
            doExtract=True, doBadPixels=True,            
            doMask=True, doOffset=True, doBias=True, doFF=True)

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
            #for iObs in range(3): # TESTS

                for jCCD, CCDk in enumerate(CCDs):
                    

                    fullccdobj_name = os.path.join(
                        picklespath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    # Cut-out "spots"

                    spots_array = np.zeros((nQuads, nSpots), dtype=object)

                    for kQ, Quad in enumerate(Quads):

                        for lS, SpotName in enumerate(SpotNames):

                            coo = polib.Point_CooNom[CCDk][Quad][SpotName]
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
        
    def opt_xtalk(self):
        """Does analysis of cross-talk using point sources as estimulators."""
        
        if self.report is not None:
            self.report.add_Section(
                keyword='xtalk', Title='Cross-Talk Analysis', level=1)
        
        DDindices = copy.deepcopy(self.dd.indices)
        nObs = DDindices.get_len('ix')
        CCDs = DDindices.get_vals('CCD')
        
        if not self.drill:
            
            for iObs in range(nObs):
                
                for jCCD, CCDkey in enumerate(CCDs):
                    
                    dpath = self.dd.mx['datapath'][iObs, jCCD]
                    
                    infits = os.path.join(dpath, '%s.fits' %
                                          self.dd.mx['File_name'][iObs, jCCD])
                    
                    ccdobj = ccd.CCD(infits)
                    
                    
                    
