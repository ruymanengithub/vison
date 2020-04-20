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

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import pandas as pd


from vison.pipe.task import HKKeys
from vison.support import context
from vison.datamodel import cdp
from vison.datamodel import scriptic as sc
from vison.support import files
#from vison.pipe.task import Task
from .FlatTask import FlatTask
from vison.datamodel import inputs, core
from vison.support import utils
from . import NL01aux
from . import nl as nllib
from vison import ogse_profiles

from pylab import plot, show
# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


NL01_commvalues = dict(program='CALCAMP', test='NL01',
                       flushes=7, siflsh=1, siflsh_p=500,
                       inisweep=1,
                       vstart=0, vend=2086,
                       exptime=0., shuttr=1, e_shuttr=0,
                       mirr_on=0,
                       wave=6,
                       motr_on=0,
                       source='flat',
                       comments='')


plusminus10pcent = 1. + np.array([-0.10, 0.10])

NL01_relfluences = np.array(
    [0.5, 0.7, 1., 2., 3., 5., 10., 20., 30., 55., 70., 80., 85., 90., 95., 100., 110.])


def get_Flu_lims(relfluences):

    reflims = 0.5 * 2**16 * plusminus10pcent

    FLU_lims = OrderedDict(CCD1=OrderedDict())
    FLU_lims['CCD1']['col001'] = np.array([-10., 20.])  # BGD. ADU
    FLU_lims['CCD1']['col002'] = reflims

    for iflu, rflu in enumerate(relfluences):
        _cenval = min(rflu / 100., 1.) * 2.**16
        _lims = _cenval * plusminus10pcent
        FLU_lims['CCD1']['col%03i' % (2 * iflu + 3,)] = _lims
        FLU_lims['CCD1']['col%03i' % (2 * iflu + 3 + 1,)] = reflims

    for i in [2, 3]:
        FLU_lims['CCD%i' % i] = copy.deepcopy(FLU_lims['CCD1'])

    return FLU_lims


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

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.Nbgd = 3
        self.Nstab0 = 1
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('extract', self.extract_stats),
                         ('NL', self.produce_NLCs),
                         ('satCTE', self.do_satCTE)]
        super(NL01, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                                   cleanafter=cleanafter)
        self.name = 'NL01'
        self.type = 'Simple'

        self.HKKeys = HKKeys
        self.CDP_lib = NL01aux.get_CDP_lib()
        self.figdict = NL01aux.get_NL01figs()
        # dict(figs='figs',pickles='ccdpickles')
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       products='products')
        self.window = dict(wpx=300, hpx=300)

    def set_inpdefaults(self, **kwargs):
        wave = 0
        tFWCw = self.ogse.profile['tFWC_flat']['nm%i' % wave]
        expts = (NL01_relfluences / 100. *
                 tFWCw).tolist()  # ms
        self.inpdefaults = dict(exptimes=expts,
                                exptinter=0.5 * tFWCw,
                                frames=(np.ones(len(expts), dtype='int32') * 4).tolist(),
                                wavelength=wave,
                                )

    def set_perfdefaults(self, **kwargs):
        super(NL01, self).set_perfdefaults(**kwargs)
        #wave = self.inputs['wavelength']
        #exptimes = np.array(self.inputs['exptimes'])
        #tsatur = self.ogse.profile['tFWC_flat']['nm%i' % wave]
        #relfluences = exptimes / tsatur * 100.
        self.perfdefaults['FLU_lims'] = get_Flu_lims(NL01_relfluences)  # dict

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

        assert len(expts) == len(frames)

        FW_ID = self.ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])

        NL01_commvalues['wave'] = FW_IDX

        NL01_sdict = dict()

        NL01_sdict['col001'] = dict(frames=self.Nbgd, exptime=0, comments='BGD')
        NL01_sdict['col002'] = dict(frames=self.Nstab0, exptime=exptinter, comments='STAB')

        for ix, ifra in enumerate(frames):

            iexp = expts[ix]

            colkeyFlu = 'col%03i' % (ix * 2 + 3,)

            NL01_sdict[colkeyFlu] = dict(
                frames=ifra, exptime=iexp, comments='Fluence%i' % (ix + 1,))

            colkeySta = 'col%03i' % (ix * 2 + 3 + 1,)

            NL01_sdict[colkeySta] = dict(
                frames=1, exptime=exptinter, comments='STAB')

        Ncols = len(list(NL01_sdict.keys()))
        NL01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(NL01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
                diffvalues = diffvalues = dict()

        NL01_sdict = sc.update_structdict(NL01_sdict, commvalues, diffvalues)

        return NL01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """Loads a list of Exposure Logs and selects exposures from test NL01.

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

        super(NL01, self).prepare_images(doExtract=True,
                                         doBadPixels=True,
                                         doMask=True,
                                         doOffset=True,
                                         doBias=True,
                                         doFF=False)

    def recalibrate_exptimes(self, exptimes):
        """Corrects exposure times given independent calibration of the shutter."""

        _path = os.path.abspath(ogse_profiles.__file__)
        ogsepath = os.path.dirname(_path)
        fpath = os.path.join(ogsepath, self.ogse.profile['SHUTTER_CALIB'])

        if not os.path.exists(fpath):

            if self.log is not None:
                self.log.info('Exposure time calibration file not found! [%s]' % fpath)

            raise RuntimeError

        newexptimes = nllib.recalibrate_exptimes(exptimes, calibrationfile=fpath)

        return newexptimes

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
            self.report.add_Text('Segmenting on %i x %i windows...' % (wpx, hpx))

        ObsIDs = self.dd.mx['ObsID'][:].copy()

        indices = copy.deepcopy(self.dd.indices)

        ishape = indices.shape

        nObs = ishape[0]
        nCCD = ishape[1]
        nQuad = ishape[2]

        Quads = indices.get_vals('Quad')
        CCDs = indices.get_vals('CCD')

        tile_coos = dict()
        for Q in Quads:
            tile_coos[Q] = self.ccdcalc.get_tile_coos(Q, wpx, hpx, noedges=True)

        Nsectors = tile_coos[Quads[0]]['Nsamps']
        sectornames = np.arange(Nsectors)

        Sindices = copy.deepcopy(self.dd.indices)
        if 'Sector' not in Sindices.names:
            Sindices.append(core.vIndex('Sector', vals=sectornames))

        # Initializing new columns

        valini = 0.

        self.dd.initColumn('sec_med', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('sec_var', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('sec_X', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('sec_Y', Sindices, dtype='float32', valini=valini)

        if not self.drill:

            dpath = self.inputs['subpaths']['ccdpickles']

            for iObs, ObsID in enumerate(ObsIDs):

                print(('Extracting from OBSID %i/%i' % (iObs + 1, nObs)))

                for jCCD, CCDk in enumerate(CCDs):

                    ccdobj_f = os.path.join(
                        dpath, '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = copy.deepcopy(files.cPickleRead(ccdobj_f))

                    for kQ in range(nQuad):

                        Q = Quads[kQ]

                        _tile_coos = tile_coos[Q]

                        _meds = ccdobj.get_tiles_stats(
                            Q, _tile_coos, 'median', extension=-1)
                        _vars = ccdobj.get_tiles_stats(
                            Q, _tile_coos, 'std', extension=-1)**2.

                        self.dd.mx['sec_med'][iObs, jCCD, kQ, :] = _meds.copy()
                        self.dd.mx['sec_var'][iObs, jCCD, kQ, :] = _vars.copy()

                        Xtiles = np.array([item[0] for item in _tile_coos['ccpix']])
                        Ytiles = np.array([item[1] for item in _tile_coos['ccpix']])

                        # tiles coordinates are NOT in canonical reference system,
                        # with Output Node in lower-left corner, so we address that

                        if Q in ['F', 'G']:
                            Xtiles = 2048. - Xtiles
                        if Q in ['E', 'F']:
                            Ytiles = 2066. - Ytiles

                        self.dd.mx['sec_X'][iObs, jCCD, kQ, :] = Xtiles.copy()
                        self.dd.mx['sec_Y'][iObs, jCCD, kQ, :] = Ytiles.copy()

                        # tests
                        # if self.dd.mx['exptime'][iObs,jCCD]>5.:
                        # if _meds.mean()>5E4:
                        #    print('%i %s %s' % (ObsID,CCDk,Q))
                        #    plot(_meds.flatten(),_vars.flatten(),'k.')
                        #    show()
                        #    stop()

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

        CCDs = dIndices.get_vals('CCD')
        Quads = dIndices.get_vals('Quad')

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

        NL_TB['CCD'] = np.zeros(NP, dtype='int32')
        NL_TB['Q'] = np.zeros(NP, dtype='int32')
        NL_TB['MAXNLPC'] = np.zeros(NP, dtype='float32')
        NL_TB['FLU_MAXNLPC'] = np.zeros(NP, dtype='float32')

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

        curves_cdp.data['labelkeys'] = ['data', 'fit']

        # Fitting the NL curves

        for iCCD, CCDkey in enumerate(CCDs):

            for jQ, Q in enumerate(Quads):

                kk = iCCD * nQ + jQ

                raw_med = self.dd.mx['sec_med'][:, iCCD, jQ, :].copy()
                raw_var = self.dd.mx['sec_var'][:, iCCD, jQ, :].copy()
                #col_labels = self.dd.mx['label'][:, iCCD].copy()
                exptimes = self.dd.mx['exptime'][:, iCCD].copy()
                dtobjs = self.dd.mx['time'][:, iCCD].copy()

                # fitresults = OrderedDict(coeffs, NLdeg, maxNLpc,flu_maxNLpc, bgd)
                _fitresults = nllib.wrap_fitNL_SingleFilter(raw_med, raw_var,
                                                            exptimes,
                                                            dtobjs,
                                                            TrackFlux=True,
                                                            subBgd=True)

                NLall_mx[CCDkey][Q].update(_fitresults)

                NL_TB['CCD'][kk] = iCCD + 1
                NL_TB['Q'][kk] = jQ + 1
                NL_TB['MAXNLPC'][kk] = _fitresults['maxNLpc']
                NL_TB['FLU_MAXNLPC'][kk] = _fitresults['flu_maxNLpc']

                curves_cdp.data[CCDkey][Q]['x']['data'] = _fitresults['inputcurve']['X'].copy()
                curves_cdp.data[CCDkey][Q]['y']['data'] = _fitresults['inputcurve']['Y'].copy()
                curves_cdp.data[CCDkey][Q]['x']['fit'] = _fitresults['outputcurve']['X'].copy()
                curves_cdp.data[CCDkey][Q]['y']['fit'] = _fitresults['outputcurve']['Y'].copy()

        self.dd.products['NL'] = copy.deepcopy(NLall_mx)

        # Build Tables

        NL_TB_dddf = OrderedDict(NL_TB=pd.DataFrame.from_dict(NL_TB))

        nl_tb_cdp = self.CDP_lib['NL_TB']
        nl_tb_cdp.path = prodspath
        nl_tb_cdp.ingest_inputs(
            data=NL_TB_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy()
        )

        nl_tb_cdp.init_wb_and_fillAll(header_title='NL01: RESULTS TABLE')
        self.save_CDP(nl_tb_cdp)
        self.pack_CDP_to_dd(nl_tb_cdp, 'NL_TB_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def ff(x): return '%.2f' % x

            formatters = [fccd, fq, ff, ff]

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

        self.canbecleaned = True

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
