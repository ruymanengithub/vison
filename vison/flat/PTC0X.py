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

"""

# IMPORT STUFF
import numpy as np
import os
import copy
from collections import OrderedDict
import pandas as pd

from vison.pipe.task import Task
from vison.support import context
from vison.datamodel import core
#from vison.pipe import lib as pilib
from vison.ogse import ogse as ogsemod
from vison.datamodel import scriptic as sc
import ptc as ptclib
from FlatTask import FlatTask
from vison.datamodel import inputs
import PTC0Xaux
from vison.support import utils
from vison.datamodel import cdp
from vison.support.files import cPickleRead
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']



PTC0X_commvalues = dict(program='CALCAMP',
                        flushes=7,siflsh=1,siflsh_p=500,                        
                        inisweep=1,
                        vstart=0,vend=2086,
                        exptime=0.,shuttr=1,e_shuttr=0,
                        mirr_on=0,
                        wave=4,
                        motr_on=0,
                        source='flat',
                        comments='')


PTC01_relfluences = np.array(
    [5., 10., 20., 30., 50., 70., 80., 90., 100., 110., 120.])

PTC02_relfluences = np.array([10., 30., 50., 70., 80., 90.])

FLATFLUX00_relfluences = np.array([5., 20., 50., 80.])

def get_testdefaults_PTC0X(ogseobj=None):

    if ogseobj is None:
        ogseobj = ogsemod.Ogse()

    tFWC800 = ogseobj.profile['tFWC_flat']['nm800']

    PTC01_exptimes = (PTC01_relfluences / 100. *
                      tFWC800).tolist()  # ms
    PTC02waves = [590, 730, 800, 880, 0]
    FLATFLUX00waves = [590, 730, 800, 880, 0]

    PTC02TEMP_exptimes = (PTC02_relfluences / 100. *
                          tFWC800).tolist()

    testdefaults = dict(PTC01=dict(exptimes=PTC01_exptimes,
                                   frames=[10, 10, 10, 10, 10,
                                           10, 10, 10, 4, 4, 4],
                                   wavelength=800),
                        PTC02WAVE=dict(waves=PTC02waves,
                                       frames=[4, 4, 4, 4, 4, 4],
                                       exptimes=dict()),
                        FLATFLUX00=dict(waves=FLATFLUX00waves,
                                       frames=[1, 1, 1, 1],
                                       exptimes=dict()),
                        PTC02TEMP=dict(frames=[4, 4, 4, 4, 4, 4],
                                       exptimes=PTC02TEMP_exptimes,
                                       wavelength=800))

    for w in testdefaults['PTC02WAVE']['waves']:
        tFWCw = ogseobj.profile['tFWC_flat']['nm%i' % w]
        testdefaults['PTC02WAVE']['exptimes']['nm%i' % w] = (
            PTC02_relfluences/100.*tFWCw).tolist()
        
    for w in testdefaults['FLATFLUX00']['waves']:
        tFWCw = ogseobj.profile['tFWC_flat']['nm%i' % w]
        testdefaults['FLATFLUX00']['exptimes']['nm%i' % w] = (
            FLATFLUX00_relfluences/100.*tFWCw).tolist()

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

FLU_lims_FLATFLUX00 = dict(CCD1=dict())
for iflu, rflu in enumerate(FLATFLUX00_relfluences):
    _cenval = min(rflu / 100., 1.) * 2.**16
    _lims = _cenval * plusminus10pcent
    FLU_lims_FLATFLUX00['CCD1']['col%i' % (iflu+1)] = _lims

for i in [2, 3]:
    FLU_lims_FLATFLUX00['CCD%i' % i] = copy.deepcopy(FLU_lims_FLATFLUX00['CCD1'])


FLU_lims_FLATFLUX00 = dict(CCD1=dict())
for iflu, rflu in enumerate(FLATFLUX00_relfluences):
    _cenval = min(rflu / 100., 1.) * 2.**16
    _lims = _cenval * plusminus10pcent
    FLU_lims_FLATFLUX00['CCD1']['col%i' % (iflu+1)] = _lims

for i in [2, 3]:
    FLU_lims_FLATFLUX00['CCD%i' % i] = copy.deepcopy(FLU_lims_FLATFLUX00['CCD1'])


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
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract', self.extract_PTC),
                         ('meta', self.meta_analysis)]

        super(PTC0X, self).__init__(inputs, log, drill, debug)
        self.name = 'PTC0X'
        self.type = 'Simple'
        self.HKKeys = HKKeys
        self.CDP_lib = PTC0Xaux.CDP_lib.copy()
        self.figdict = PTC0Xaux.gt_PTC0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                   products='products')
        self.window = dict(wpx=300,hpx=300)

    def set_inpdefaults(self, **kwargs):
        """ """

        testdefaults = get_testdefaults_PTC0X(self.ogse)

        try:
            testkey = kwargs['test']
        except KeyError:
            testkey = 'PTC01'

        if testkey == 'PTC01':
            _testkey = 'PTC01' 
        if 'FLATFLUX00' in testkey:
            _testkey = 'FLATFLUX00'
        
        if 'PTC02' in testkey:
            if testkey[-1] == 'K':
                _testkey = 'PTC02TEMP'
            else:
                _testkey = 'PTC02WAVE'
        if 'FLATFLUX00' in testkey:
            _testkey = 'FLATFLUX00'

        if _testkey in ['PTC02WAVE','FLATFLUX00']:
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
        elif 'FLATFLUX00' in testkey:
            FLU_lims = FLU_lims_FLATFLUX00

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

    def filterexposures(self, structure, explog, OBSID_lims):
        """

        """
        
        return super(PTC0X, self).filterexposures(structure, explog, OBSID_lims, 
                    wavedkeys=['motr_siz'],colorblind=False)

    def prepare_images(self):
        Task.prepare_images(self, doExtract=True, doMask=True,
                                         doOffset=True, doBias=False, doFF=False)
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
        
        # HARDWIRED VALUES
        wpx = self.window['wpx']
        hpx = self.window['hpx']

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='PTC Extraction', level=0)
            self.report.add_Text('Segmenting on %i x %i windows...' % (wpx,hpx))

        # labels should be the same accross CCDs. PATCH.
        label = self.dd.mx['label'][:, 0].copy()
        ObsIDs = self.dd.mx['ObsID'][:].copy()

        indices = copy.deepcopy(self.dd.indices)

        nObs, nCCD, nQuad = indices.shape[0:3]

        Quads = indices.get_vals('Quad')
        CCDs = indices.get_vals('CCD')

        tile_coos = dict()
        for Quad in Quads:
            tile_coos[Quad] = self.ccdcalc.get_tile_coos(Quad, wpx, hpx)
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
                        dpath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])
                    ccdobj_eve_f = os.path.join(
                        dpath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs_pair, jCCD])

                    ccdobj_odd = copy.deepcopy(
                        cPickleRead(ccdobj_odd_f))  
                    ccdobj_eve = copy.deepcopy(
                        cPickleRead(ccdobj_eve_f)) 
                    
                    evedata = ccdobj_eve.extensions[-1].data.copy()

                    # easy way to subtract one image from the other
                    ccdobj_sub = copy.deepcopy(ccdobj_odd)
                    ccdobj_sub.sub_bias(evedata, extension=-1)

                    for kQ in range(nQuad):

                        Quad = Quads[kQ]

                        _tile_coos = tile_coos[Quad]

                        _meds = ccdobj_odd.get_tiles_stats(
                            Quad, _tile_coos, 'median', extension=-1)
                        
                        # IT'S A SUBTRACTION, SO WE HAVE TO DIVIDE BY 2 THE VARIANCE!
                        _vars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, 'std', extension=-1)**2./2. 

                        self.dd.mx['sec_med'][iObs, jCCD, kQ, :] = _meds.copy()
                        self.dd.mx['sec_var'][iObs, jCCD, kQ, :] = _vars.copy()
                        
        # MISSING: any figure?
        # MISSING: any Table?
        

        return 

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
        nC = len(CCDs)
        Quads = dIndices.get_vals('Quad')
        nQ = len(Quads)
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        prodspath = self.inputs['subpaths']['products']

        # Initializations of output data-products
        
        NP = nC * nQ
        
        GAIN_TB = OrderedDict()
        
        GAIN_TB['CCD'] = np.zeros(NP,dtype='int32')
        GAIN_TB['Q'] = np.zeros(NP,dtype='int32')
        GAIN_TB['gain'] = np.zeros(NP,dtype='float32')
        GAIN_TB['egain'] = np.zeros(NP,dtype='float32')
        GAIN_TB['alpha'] = np.zeros(NP,dtype='float32')
        GAIN_TB['fqual'] = np.zeros(NP,dtype='int32')
        GAIN_TB['bloom'] = np.zeros(NP,dtype='float32')
        
        curves_cdp = cdp.CDP()
        curves_cdp.header = CDP_header.copy()
        curves_cdp.path = prodspath
        curves_cdp.data = OrderedDict()
        
        
        for CCDk in CCDs:
            curves_cdp.data[CCDk] = OrderedDict()            
            for Q in Quads:
                curves_cdp.data[CCDk][Q] = OrderedDict()
                curves_cdp.data[CCDk][Q]['x'] = OrderedDict()
                curves_cdp.data[CCDk][Q]['y'] = OrderedDict()
        
        curves_cdp.data['labelkeys'] = ['data','fit']
        
        gain_mx = OrderedDict()

        g_tmp_keys = ['a0', 'ea0', 'a1', 'ea1', 'a2',
                      'ea2', 'gain', 'egain', 'alpha', 'rn']

        for CCDk in CCDs:
            gain_mx[CCDk] = dict()
            for Q in Quads:
                gain_mx[CCDk][Q] = dict()
                for key in g_tmp_keys:
                    gain_mx[CCDk][Q][key] = np.nan

        b_tmp_keys = ['bloom_ADU', 'bloom_e']

        bloom_mx = OrderedDict()

        for CCDkey in CCDs:
            bloom_mx[CCDkey] = dict()
            for Quad in Quads:
                bloom_mx[CCDkey][Quad] = dict()
                for key in b_tmp_keys:
                    bloom_mx[CCDkey][Quad][key] = np.nan

        # fitting the PTCs

        for iCCD, CCDk in enumerate(CCDs):

            for jQ, Q in enumerate(Quads):
                
                
                ixsel = np.where(~np.isnan(self.dd.mx['ObsID_pair'][:]))

                raw_var = self.dd.mx['sec_var'][ixsel, iCCD, jQ, :]
                raw_med = self.dd.mx['sec_med'][ixsel, iCCD, jQ, :]
                
                ixnonan = np.where(~np.isnan(raw_var) & ~np.isnan(raw_med))
                var = raw_var[ixnonan]
                med = raw_med[ixnonan]

                #fitresults = dict(fit=p,efit=ep,gain=g,cuadterm=cuadterm,rn=rn,badresult=badresult)

                _fitresults = ptclib.fitPTC(med, var)

                for zx in range(2+1):
                    gain_mx[CCDk][Q]['a%i' % zx] = _fitresults['fit'][2-zx]
                    gain_mx[CCDk][Q]['ea%i' % zx] = _fitresults['efit'][2-zx]

                gain_mx[CCDk][Q]['gain'] = _fitresults['gain']
                gain_mx[CCDk][Q]['egain'] = _fitresults['efit'][1]/_fitresults['gain']**2.
                gain_mx[CCDk][Q]['alpha'] = -_fitresults['quadterm']
                gain_mx[CCDk][Q]['rn'] = _fitresults['rn']
                gain_mx[CCDk][Q]['quality'] = _fitresults['quality']
                
                _bloom = ptclib.foo_bloom(med, var)

                bloom_mx[CCDk][Q]['bloom_ADU'] = _bloom['bloom']
                
                kk = iCCD * nQ + jQ
                
                GAIN_TB['CCD'][kk] = iCCD+1
                GAIN_TB['Q'][kk] = jQ+1
                GAIN_TB['gain'][kk] = gain_mx[CCDk][Q]['gain']
                GAIN_TB['egain'][kk] = gain_mx[CCDk][Q]['egain']
                GAIN_TB['alpha'][kk] = gain_mx[CCDk][Q]['alpha']
                GAIN_TB['fqual'][kk] = gain_mx[CCDk][Q]['quality']
                GAIN_TB['bloom'][kk] = bloom_mx[CCDk][Q]['bloom_ADU']
                
                curves_cdp.data[CCDk][Q]['x']['data'] = med.copy()
                curves_cdp.data[CCDk][Q]['y']['data'] = var.copy()
                
                fkmed = np.linspace(0.,2.**16,100)
                bfvar = np.polyval(_fitresults['fit'],fkmed)
                
                curves_cdp.data[CCDk][Q]['x']['fit'] = fkmed.copy()
                curves_cdp.data[CCDk][Q]['y']['fit'] = bfvar.copy()

        self.dd.products['gain_mx'] = copy.deepcopy(gain_mx)
        self.dd.products['bloom_mx'] = copy.deepcopy(bloom_mx)

        # Build Tables
        
        GAIN_TB_dddf = OrderedDict(GAIN_TB = pd.DataFrame.from_dict(GAIN_TB))
        
        gain_tb_cdp = self.CDP_lib['GAIN_TB']
        gain_tb_cdp.rootname = gain_tb_cdp.rootname % \
          (self.inputs['test'],self.inputs['wavelength'])
        gain_tb_cdp.path = self.inputs['subpaths']['products']
        gain_tb_cdp.ingest_inputs(
                data = GAIN_TB_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )

        gain_tb_cdp.init_wb_and_fillAll(header_title='%s (%inm): PTC TABLE' % \
                (self.inputs['test'],self.inputs['wavelength']))
        self.save_CDP(gain_tb_cdp)
        self.pack_CDP_to_dd(gain_tb_cdp, 'GAIN_TB_CDP')

        if self.report is not None:
            
            fccd = lambda x: CCDs[x-1]
            fq = lambda x: Quads[x-1]
            fi = lambda x: '%i' % x
            ff = lambda x: '%.2f' % x
            fE = lambda x: '%.2E' % x
            
            cov_formatters=[fccd,fq,ff,fE,fE,fi,fE]
            
            caption = '%s (%inm): PTC TABLE' % \
                (self.inputs['test'],self.inputs['wavelength'])
            Gtex = gain_tb_cdp.get_textable(sheet='GAIN_TB', caption=caption,
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=cov_formatters)
            
            
            self.report.add_Text(Gtex)        
        
        
        # Do plots
        
        fdict_PTC = self.figdict['PTC0X_PTC_curves'][1]
        fdict_PTC['data'] = curves_cdp.data.copy()
        if self.report is not None:
            self.addFigures_ST(figkeys=['PTC0X_PTC_curves'], 
                               dobuilddata=False)
        

       
