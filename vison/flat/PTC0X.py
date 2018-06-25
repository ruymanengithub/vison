# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PTC_0X

Photon-Transfer-Curve Analysis
   PTC01 - nominal temperature and wavelength
   PTC02 - alternative temperatures / wavelengths

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract pairs of exposures with equal fluence
    - Synoptic analysis:
        variance vs. fluence
        variance(binned difference-frames) vs. fluence
    - extract: RON, gain, gain(fluence)
    - produce synoptic figures
    - Save results.



Created on Mon Apr  3 17:00:24 2017

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
from vison.ogse import ogse as ogsemod
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import HKtools
from vison.datamodel import core, ccd
from vison.image import calibration
import ptc as ptclib
from vison.image import performance
from FlatTask import FlatTask
from vison.datamodel import inputs
import PTC0Xaux
from vison.support.files import cPickleRead, cPickleDumpDictionary
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

PTC0X_commvalues = dict(program='CALCAMP',
                        flushes=7,
                        exptime=0.,
                        vstart=0,
                        vend=2086,
                        shuttr=1,
                        siflsh=1,
                        siflsh_p=500,
                        wave=4,
                        source='flat',
                        comments='')


PTC01_relfluences = np.array(
    [5., 10., 20., 30., 50., 70., 80., 90., 100., 110., 120.])

PTC02_relfluences = np.array([10., 30., 50., 70., 80., 90.])


def get_testdefaults_PTC0X(ogseobj=None):

    if ogseobj is None:
        ogseobj = ogsemod.Ogse()

    tFWC800 = ogseobj.profile['tFWC_flat']['nm800']

    PTC01_exptimes = (PTC01_relfluences / 100. *
                      tFWC800).tolist()  # ms
    PTC02waves = [590, 640, 730, 800, 880, 0]

    PTC02TEMP_exptimes = (PTC02_relfluences / 100. *
                          tFWC800).tolist()

    testdefaults = dict(PTC01=dict(exptimes=PTC01_exptimes,
                                   frames=[10, 10, 10, 10, 10,
                                           10, 10, 10, 4, 4, 4],
                                   wavelength=800),
                        PTC02WAVE=dict(waves=PTC02waves,
                                       frames=[4, 4, 4, 4, 4, 4],
                                       exptimes=dict()),
                        PTC02TEMP=dict(frames=[4, 4, 4, 4, 4, 4],
                                       exptimes=PTC02TEMP_exptimes,
                                       wavelength=800))

    for w in testdefaults['PTC02WAVE']['waves']:
        tFWCw = ogseobj.profile['tFWC_flat']['nm%i' % w]
        testdefaults['PTC02WAVE']['exptimes']['nm%i' % w] = (
            PTC02_relfluences/100.*tFWCw).tolist()

    return testdefaults


plusminus10pcent = 1.+np.array([-0.10, 0.10])

FLU_lims_PTC01 = dict(CCD1=dict())
for iflu, rflu in enumerate(PTC01_relfluences):
    _cenval = min(rflu / 100., 1.) * 2.**16
    _lims = _cenval * plusminus10pcent
    FLU_lims_PTC01['CCD1']['col%i' % (iflu+1)] = _lims

for i in [2, 3]:
    FLU_lims_PTC01['CCD%i' % i] = copy.deepcopy(FLU_lims_PTC01['CCD1'])

FLU_lims_PTC02 = dict(CCD1=dict())
for iflu, rflu in enumerate(PTC02_relfluences):
    _cenval = min(rflu / 100., 1.) * 2.**16
    _lims = _cenval * plusminus10pcent
    FLU_lims_PTC02['CCD1']['col%i' % (iflu+1)] = _lims

for i in [2, 3]:
    FLU_lims_PTC02['CCD%i' % i] = copy.deepcopy(FLU_lims_PTC02['CCD1'])


class PTC0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
    ])))


class PTC0X(FlatTask):
    """ """

    inputsclass = PTC0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(PTC0X, self).__init__(inputs, log, drill, debug)
        self.name = 'PTC0X'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract', self.extract_PTC),
                         ('meta', self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = PTC0Xaux.gt_PTC0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles')

    def set_inpdefaults(self, **kwargs):
        """ """

        testdefaults = get_testdefaults_PTC0X(self.ogse)

        try:
            testkey = kwargs['test']
        except KeyError:
            testkey = 'PTC01'

        if testkey == 'PTC01':
            _testkey = 'PTC01'
        if 'PTC02' in testkey:
            if testkey[-1] == 'K':
                _testkey = 'PTC02TEMP'
            else:
                _testkey = 'PTC02WAVE'

        if _testkey == 'PTC02WAVE':
            try:
                wavelength = kwargs['wavelength']
            except KeyError:
                wavelength = 800
            exptimes = testdefaults[_testkey]['exptimes']['nm%i' % wavelength]
        else:
            exptimes = testdefaults[_testkey]['exptimes']
            wavelength = testdefaults[_testkey]['wavelength']

        frames = testdefaults[_testkey]['frames']

        self.inpdefaults = dict(test=testkey, wavelength=wavelength,
                                frames=frames, exptimes=exptimes)

    def set_perfdefaults(self, **kwargs):
        super(PTC0X, self).set_perfdefaults(**kwargs)

        try:
            testkey = kwargs['test']
        except KeyError:
            testkey = 'PTC01'

        if 'PTC01' in testkey:
            FLU_lims = FLU_lims_PTC01
        elif 'PTC02' in testkey:
            FLU_lims = FLU_lims_PTC02

        self.perfdefaults['FLU_lims'] = FLU_lims  # dict

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds PTC0X script structure dictionary.

        #:param exptimes: list of ints [ms], exposure times.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 800 nm.
        :param diffvalues: dict, opt, differential values.   

        """

        testkey = self.inputs['test']
        exptimes = self.inputs['exptimes']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']

        assert len(exptimes) == len(frames)

        FW_ID = self.ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])

        PTC0X_commvalues['test'] = testkey
        PTC0X_commvalues['wave'] = FW_IDX

        PTC0X_sdict = dict()

        for ix, ifra in enumerate(frames):
            iexp = exptimes[ix]

            colkey = 'col%i' % (ix+1,)

            PTC0X_sdict[colkey] = dict(frames=ifra, exptime=iexp)

        Ncols = len(PTC0X_sdict.keys())
        PTC0X_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(PTC0X_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        PTC0X_sdict = sc.update_structdict(PTC0X_sdict, commvalues, diffvalues)

        return PTC0X_sdict

    def filterexposures(self, structure, explogf, datapath, OBSID_lims):
        """

        """
        wavedkeys = ['motr_siz']
        return super(PTC0X, self).filterexposures(structure, explogf, datapath, OBSID_lims, colorblind=False,
                                                  wavedkeys=wavedkeys)

    def prepare_images(self):
        super(PTC0X, self).prepare_images(doExtract=True, doMask=True,
                                          doOffset=True, doBias=False, doFF=False)

    def extract_PTC(self):
        """

        Performs basic analysis of images:
            - builds PTC curves: both on non-binned and binned images

        **METACODE**

        ::

            create list of OBSID pairs

            create segmentation map given grid parameters

            f.e. OBSID pair:
                CCD:
                    Q:
                        subtract CCD images
                        f.e. segment:
                            measure central value
                            measure variance

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='PTC Extraction', level=0)

        # HARDWIRED VALUES
        wpx = 300
        hpx = 300

        # labels should be the same accross CCDs. PATCH.
        label = self.dd.mx['label'][:, 0].copy()
        ObsIDs = self.dd.mx['ObsID'][:].copy()

        indices = copy.deepcopy(self.dd.indices)

        nObs, nCCD, nQuad = indices.shape

        Quads = indices.get_vals('Quad')
        CCDs = indices.get_vals('CCD')

        emptyccdobj = ccd.CCD()
        tile_coos = dict()
        for Quad in Quads:
            tile_coos[Quad] = emptyccdobj.get_tile_coos(Quad, wpx, hpx)
        Nsectors = tile_coos[Quads[0]]['Nsamps']
        sectornames = np.arange(Nsectors)

        Sindices = copy.deepcopy(self.dd.indices)
        if 'Sector' not in Sindices.names:
            Sindices.append(core.vIndex('Sector', vals=sectornames))

        # Initializing new columns

        valini = 0.

        self.dd.initColumn('sec_med', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('sec_var', Sindices, dtype='float32', valini=valini)

        # Pairing ObsIDs

        self.dd.initColumn(
            'ObsID_pair', self.dd.mx['ObsID'].indices, dtype='int64', valini=np.nan)

        ulabels = np.unique(label)

        for ulabel in ulabels:
            six = np.where(label == ulabel)
            nsix = len(six[0])
            ixeven = np.arange(0, nsix, 2)
            ixodd = np.arange(1, nsix, 2)

            self.dd.mx['ObsID_pair'][six[0][ixeven]] = ObsIDs[six[0][ixodd]]

        if not self.drill:

            # if self.proc_histo['Masked']:
            #    estimators = dict(median=np.ma.median,std=np.ma.std)
            # else:
            #    estimators = dict(median=np.median,std=np.std)

            dpath = self.inputs['subpaths']['ccdpickles']

            for iObs in range(nObs):

                _ObsID_pair = self.dd.mx['ObsID_pair'][iObs]
                if np.isnan(_ObsID_pair):
                    continue
                iObs_pair = np.where(ObsIDs == _ObsID_pair)[0][0]

                for jCCD, CCDk in enumerate(CCDs):

                    ccdobj_odd_f = os.path.join(
                        dpath, self.dd.mx['ccdobj_name'][iObs, jCCD])
                    ccdobj_eve_f = os.path.join(
                        dpath, self.dd.mx['ccdobj_name'][iObs_pair, jCCD])

                    ccdobj_odd = copy.deepcopy(
                        cPickleRead(ccdobj_odd_f))  # ['ccdobj'])
                    ccdobj_eve = copy.deepcopy(
                        cPickleRead(ccdobj_eve_f))  # ['ccdobj'])

                    evedata = ccdobj_eve.extensions[-1].data.copy()

                    # easy way to subtract one image from the other
                    ccdobj_odd.sub_bias(evedata, extension=-1)

                    for kQ in range(nQuad):

                        Quad = Quads[kQ]

                        _tile_coos = tile_coos[Quad]

                        _meds = ccdobj_odd.get_tile_stats(
                            Quad, _tile_coos, 'median', extension=-1)
                        _vars = ccdobj_odd.get_tile_stats(
                            Quad, _tile_coos, 'std', extension=-1)**2.

                        self.dd.mx['sec_med'][iObs, jCCD, kQ, :] = _meds.copy()
                        self.dd.mx['sec_var'][iObs, jCCD, kQ, :] = _vars.copy()

        return None

    def meta_analysis(self):
        """

        Analyzes the variance and fluence:
        gain, and gain(fluence)

        METACODE

        ::

            f.e. CCD:
                Q:
                    (using stats across segments:)
                    fit PTC to quadratic model
                    solve for gain
                    solve for alpha (pixel-correls, Guyonnet+15)
                    solve for blooming limit (ADU)
                        convert bloom limit to electrons, using gain

            plot PTC curves with best-fit f.e. CCD, Q
            report on gain estimates f. e. CCD, Q (table)
            report on blooming limits (table)

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='PTC Analysis', level=0)

        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices.get_vals('CCD')
        Quads = dIndices.get_vals('Quad')

        # Initializations of output data-products

        gain_mx = OrderedDict()

        g_tmp_keys = ['a0', 'ea0', 'a1', 'ea1', 'a2',
                      'ea2', 'gain', 'egain', 'alpha', 'rn']

        for CCDkey in CCDs:
            gain_mx[CCDkey] = dict()
            for Quad in Quads:
                gain_mx[CCDkey][Quad] = dict()
                for key in g_tmp_keys:
                    gain_mx[CCDkey][Quad][key] = np.nan

        b_tmp_keys = ['bloom_ADU', 'bloom_e']

        bloom_mx = OrderedDict()

        for CCDkey in CCDs:
            bloom_mx[CCDkey] = dict()
            for Quad in Quads:
                bloom_mx[CCDkey][Quad] = dict()
                for key in b_tmp_keys:
                    bloom_mx[CCDkey][Quad][key] = np.nan

        # fitting the PTCs

        for iCCD, CCDkey in enumerate(CCDs):

            for jQ, Q in enumerate(Quads):

                raw_var = self.dd.mx['sec_var'][:, iCCD, jQ, :]
                raw_med = self.dd.mx['sec_med'][:, iCCD, jQ, :]

                ixnonan = np.where(~np.isnan(raw_var) & ~np.isnan(raw_med))
                var = raw_var[ixnonan]
                med = raw_med[ixnonan]

                #fitresults = dict(fit=p,efit=ep,gain=g,cuadterm=cuadterm,rn=rn,badresult=badresult)

                _fitresults = ptclib.fitPTC(med, var)

                for zx in range(2):
                    gain_mx[CCDkey][Q]['a%i' % zx] = _fitresults['fit'][2-zx]
                    gain_mx[CCDkey][Q]['ea%i' % zx] = _fitresults['efit'][2-zx]

                gain_mx[CCDkey][Q]['gain'] = _fitresults['gain']
                gain_mx[CCDkey][Q]['rn'] = _fitresults['rn']
                gain_mx[CCDkey][Q]['quality'] = _fitresults['quality']

                _bloom = ptclib.foo_bloom(med, var)

                bloom_mx[CCDkey][Q]['bloom_ADU'] = _bloom['bloom']
                bloom_mx[CCDkey][Q]['bloom_ADU'] = _bloom['bloom']

        self.dd.products['gain'] = copy.deepcopy(gain_mx)
        self.dd.products['bloom'] = copy.deepcopy(bloom_mx)

        # Build Tables

        # Do plots

        # Add reports
