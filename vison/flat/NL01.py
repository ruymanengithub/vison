# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: NL01

End-To-End Non-Linearity Curve


Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract offset level.
    - Divide by Flat-field.
    - Synoptic analysis:
        fluence ratios vs. extime ratios >> non-linearity curve
    - extract: Non-Linearity curve for each CCD and quadrant
    - produce synoptic figures
    - Save results.


Created on Mon Apr  3 17:38:00 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import pandas as pd

from vison.support import context
from vison.datamodel import cdp
from vison.datamodel import scriptic as sc
from vison.support import files
#from vison.pipe.task import Task
from FlatTask import FlatTask
from vison.datamodel import inputs, core
from vison.support import utils
import NL01aux
import nl as nllib
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


NL01_commvalues = dict(program='CALCAMP',
                       test='NL01',
                       IPHI1=1, IPHI2=1, IPHI3=1, IPHI4=0,
                       rdmode='fwd_bas',
                       flushes=7, vstart=0, vend=2086,
                       exptime=0., shuttr=1,
                       siflsh=1, siflsh_p=500,
                       wave=5,
                       source='flat',
                       comments='')


plusminus10pcent = 1.+np.array([-0.10, 0.10])

NL01_relfluences = np.array(
    [1.,2.,3.,5., 10., 20., 30., 50., 70., 80., 85., 90., 95., 100., 110.])

FLU_lims = dict(CCD1=dict())
for iflu, rflu in enumerate(NL01_relfluences):
    _cenval = min(rflu / 100., 1.) * 2.**16
    _lims = _cenval * plusminus10pcent
    FLU_lims['CCD1']['col%i' % (iflu+1)] = _lims

for i in [2, 3]:
    FLU_lims['CCD%i' % i] = copy.deepcopy(FLU_lims['CCD1'])


class NL01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('exptinter', ([
            float], 'Exposure time for interleaved fluence stability control frames.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength'))
    ])))


class NL01(FlatTask):
    """ """

    inputsclass = NL01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('extract', self.extract_stats),
                         ('NL', self.produce_NLCs),
                         ('satCTE', self.do_satCTE)]
        super(NL01, self).__init__(inputs, log, drill, debug)
        self.name = 'NL01'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.CDP_lib = NL01aux.CDP_lib.copy()
        self.figdict = NL01aux.NL01figs.copy()
        # dict(figs='figs',pickles='ccdpickles')
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                   products='products')
        self.window = dict(wpx=300,hpx=300)

    def set_inpdefaults(self, **kwargs):

        wave = 880

        tFWCw = self.ogse.profile['tFWC_flat']['nm%i' % wave]

        expts = (NL01_relfluences/100. *
                 tFWCw).tolist()  # ms
        self.inpdefaults = dict(exptimes=expts,
                                exptinter=0.5 * tFWCw,
                                frames=(np.ones(len(expts), dtype='int32')*5).tolist(),
                                wavelength=wave,
                                )

    def set_perfdefaults(self, **kwargs):
        super(NL01, self).set_perfdefaults(**kwargs)
        self.perfdefaults['FLU_lims'] = FLU_lims  # dict

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds NL01 script structure dictionary.

        #:param expts: list of ints [ms], exposure times.
        #:param exptinter: int, ms, exposure time of interleaved source-stability exposures.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 0 (Neutral Density Filter)
        :param diffvalues: dict, opt, differential values.
        """

        expts = self.inputs['exptimes']
        exptinter = self.inputs['exptinter']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        
        Nbgd=3
        Nstab0 = 1

        assert len(expts) == len(frames)

        FW_ID = self.ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])

        NL01_commvalues['wave'] = FW_IDX

        NL01_sdict = dict()

        NL01_sdict['col1'] = dict(frames=Nbgd, exptime=0, comment='BGD')
        NL01_sdict['col2'] = dict(frames=Nstab0, exptime=exptinter, comment='STAB')

        for ix, ifra in enumerate(frames):

            iexp = expts[ix]

            colkeyFlu = 'col%i' % (ix*2+3,)

            NL01_sdict[colkeyFlu] = dict(
                frames=ifra, exptime=iexp, comment='Fluence%i' % (ix+1,))

            colkeySta = 'col%i' % (ix*2+3+1,)

            NL01_sdict[colkeySta] = dict(
                frames=1, exptime=exptinter, comment='STAB')

        Ncols = len(NL01_sdict.keys())
        NL01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(NL01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        NL01_sdict = sc.update_structdict(NL01_sdict, commvalues, diffvalues)

        return NL01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """Loads a list of Exposure Logs and selects exposures from test PSF0X.

        The filtering takes into account an expected structure for the 
        acquisition script.

        The datapath becomes another column in DataDict. This helps dealing
        with tests that run overnight and for which the input data is in several
        date-folders.

        """
        wavedkeys = []
        return super(NL01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                 wavedkeys=wavedkeys)

    def prep_data(self):
        """

        Takes Raw Data and prepares it for further analysis. 

        **METACODE**

        ::

            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        subtract offset
                    opt: [sub bias frame]
                    opt: [divide by FF]
                    opt: [mask-out defects]

        """
        super(NL01, self).prepare_images(doExtract=True, doMask=True,
                                         doOffset=True, doBias=True, 
                                         doFF=True)

    def extract_stats(self):
        """

        Performs basic analysis: extracts statistics from 
        image regions to later build NLC.

        **METACODE**

        ::

            create segmentation map given grid parameters    

            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        f.e. "img-segment": (done elsewhere)
                            measure central value
                            measure variance

        """
        
        # HARDWIRED VALUES
        wpx = self.window['wpx']
        hpx = self.window['hpx']

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='Image Extraction', level=0)
            self.report.add_Text('Segmenting on %i x %i windows...' % (wpx,hpx))    
        
        
        ObsIDs = self.dd.mx['ObsID'][:].copy()

        indices = copy.deepcopy(self.dd.indices)

        nObs, nCCD, nQuad = indices.shape

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

        if not self.drill:

            dpath = self.inputs['subpaths']['ccdpickles']

            for iObs, ObsID in enumerate(ObsIDs):

                for jCCD, CCDk in enumerate(CCDs):

                    ccdobj_f = os.path.join(
                        dpath, self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = copy.deepcopy(files.cPickleRead(ccdobj_f))

                    for kQ in range(nQuad):

                        Quad = Quads[kQ]

                        _tile_coos = tile_coos[Quad]

                        _meds = ccdobj.get_tile_stats(
                            Quad, _tile_coos, 'median', extension=-1)
                        _vars = ccdobj.get_tile_stats(
                            Quad, _tile_coos, 'std', extension=-1)**2.

                        self.dd.mx['sec_med'][iObs, jCCD, kQ, :] = _meds.copy()
                        self.dd.mx['sec_var'][iObs, jCCD, kQ, :] = _vars.copy()

    def produce_NLCs(self):
        """ 

        **METACODE**

        ::

            Obtains Best-Fit Non-Linearity Curve

            f.e. CCD:
                f.e. Q:

                    [opt] apply correction for source variability (interspersed exposure 
                      with constant exptime)
                    Build NL Curve (NLC) - use stats and exptimes
                    fit poly. shape to NL curve

            plot NL curves for each CCD, Q
            report max. values of NL (table)

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='NL', Title='Non-Linearity Analysis', level=0)

        dIndices = copy.deepcopy(self.dd.indices)
        CCDs = dIndices.get_vals('CCDs')
        Quads = dIndices.get_vals('Quads')
        
        nC = len(CCDs)
        nQ = len(Quads)        
        NP = nC * nQ
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        prodspath = self.inputs['subpaths']['products']
        
        # INITIALISATIONS
        
        
        # NON-LINEARITY TABLE
        
        NL_TB = OrderedDict()
        
        NL_TB['CCD'] = np.zeros(NP,dtype='int32')
        NL_TB['Q'] = np.zeros(NP,dtype='int32')
        NL_TB['MAXNLPC'] = np.zeros(NP,dtype='float32')
        NL_TB['FLU_MAXNLPC'] = np.zeros(NP,dtype='float32')
        
        
        # NON LINEARITY RESULTS

        NLall_mx = OrderedDict()

        for CCDkey in CCDs:
            NLall_mx[CCDkey] = OrderedDict()
            for Quad in Quads:
                NLall_mx[CCDkey][Quad] = OrderedDict()
                
        # NL CURVES
        
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

        # Fitting the NL curves

        for iCCD, CCDkey in enumerate(CCDs):

            for jQ, Q in enumerate(Quads):
                
                
                kk = iCCD * nQ + jQ

                raw_med = self.dd.mx['sec_med'][:, iCCD, jQ, :].copy()
                col_labels = self.dd.mx['label'][:, iCCD].copy()
                exptimes = self.dd.mx['exptime'][:, iCCD].copy()
                dtobjs = self.dd.mx['time'][:, iCCD].copy()
                
                # fitresults = OrderedDict(coeffs, NLdeg, maxNLpc,flu_maxNLpc, bgd)
                _fitresults = nllib.wrap_fitNL(raw_med, exptimes, col_labels, dtobjs, TrackFlux=True,
                                               subBgd=True)
                
                NLall_mx[CCDkey][Quad].update(_fitresults)
                
                NL_TB['CCD'][kk] = iCCD+1
                NL_TB['Q'][kk] = jQ+1
                NL_TB['MAXNLPC'] = _fitresults['maxNLpc']
                NL_TB['FLU_MAXNLPC'] = _fitresults['flu_maxNLpc']
                
                curves_cdp.data[CCDkey][Q]['x']['data'] = _fitresults['inputcurve']['YL'].copy()
                curves_cdp.data[CCDkey][Q]['y']['data'] = _fitresults['inputcurve']['Z'].copy()
                curves_cdp.data[CCDkey][Q]['x']['fit'] = _fitresults['outputcurve']['YL'].copy()
                curves_cdp.data[CCDkey][Q]['y']['fit'] = _fitresults['outputcurve']['Z'].copy()
                

        self.dd.products['NL'] = copy.deepcopy(NLall_mx)

        # Build Tables
        
        NL_TB_dddf = OrderedDict(NL_TB = pd.DataFrame.from_dict(NL_TB))
        
        nl_tb_cdp = self.CDP_lib['NL_TB']
        nl_tb_cdp.path = prodspath
        nl_tb_cdp.ingest_inputs(
                data = NL_TB_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )

        nl_tb_cdp.init_wb_and_fillAll(header_title='NL01: RESULTS TABLE')
        self.save_CDP(nl_tb_cdp)
        self.pack_CDP_to_dd(nl_tb_cdp, 'NL_TB_CDP')
        
        if self.report is not None:
            
            fccd = lambda x: CCDs[x-1]
            fq = lambda x: Quads[x-1]
            ff = lambda x: '%.2f' % x
            
            formatters=[fccd,fq,ff,ff]
            
            caption = 'NL01 results TABLE' 
            Ntex = nl_tb_cdp.get_textable(sheet='NL_TB', caption=caption,
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=formatters)
            
            
            self.report.add_Text(Ntex)        

        # Do plots
        
        

        fdict_NL = self.figdict['NL01_fit_curves'][1]
        fdict_NL['data'] = curves_cdp.data.copy()
        if self.report is not None:
            self.addFigures_ST(figkeys=['NL01_fit_curves'], 
                               dobuilddata=False)
        
        
        
        

    def do_satCTE(self):
        """

        **METACODE**

        ::

            select ObsIDs with fluence(exptime) >~ 0.5 FWC

            f.e. ObsID: 
                CCD: 
                    Q:
                        measure CTE from amount of charge in over-scan relative to fluence

            f.e. CCD: 
                Q:
                    get curve of CTE vs. fluence
                    measure FWC from curve in ADU

            report FWCs in electrons [via gain in inputs] f.e. CCD, Q (table)

        """
        raise NotImplementedError
