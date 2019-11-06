#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:01:36 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
from collections import OrderedDict
import string as st
import os
import matplotlib.cm as cm

from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files


from matplotlib import pyplot as plt
plt.switch_backend('TkAgg')
from matplotlib.colors import Normalize
# END IMPORT

cols2keep = [
    'test',
    'sn_ccd1',
    'sn_ccd2',
    'sn_ccd3',
    'sn_roe',
    'sn_rpsu',
    'exptime',
    'vstart',
    'vend',
    'rdmode',
    'flushes',
    'siflsh',
    'siflsh_p',
    'swellw',
    'swelldly',
    'inisweep',
    'cdpu_clk',
    'chinj',
    'chinj_on',
    'chinj_of',
    'id_wid',
    'id_dly',
    'chin_dly',
    'v_tpump',
    's_tpump',
    'v_tp_mod',
    's_tp_mod',
    'v_tp_cnt',
    's_tp_cnt',
    'dwell_v',
    'dwell_s',
    'toi_fl',
    'toi_tp',
    'toi_ro',
    'toi_ch',
    'motr',
    'motr_cnt',
    'motr_siz',
    'source',
    'wave',
    'mirr_on',
    'mirr_pos',
    'R1C1_TT',
    'R1C1_TB',
    'R1C2_TT',
    'R1C2_TB',
    'R1C3_TT',
    'R1C3_TB',
    'IDL',
    'IDH',
    'IG1_1_T',
    'IG1_2_T',
    'IG1_3_T',
    'IG1_1_B',
    'IG1_2_B',
    'IG1_3_B',
    'IG2_T',
    'IG2_B',
    'OD_1_T',
    'OD_2_T',
    'OD_3_T',
    'OD_1_B',
    'OD_2_B',
    'OD_3_B',
    'RD_T',
    'RD_B',
    'time',
    'chk_NPIXOFF',
    'chk_NPIXSAT',
    'offset_pre',
    'offset_ove',
    'deltaoff_pre',
    'deltaoff_ove',
    'flu_med_img',
    'flu_std_img',
    'std_pre',
    'std_ove']


class MetaBF(MetaCal):
    """ """

    def __init__(self, *args, **kwargs):
        """ """

        super(MetaBF, self).__init__(*args, **kwargs)

        self.testnames = ['BF01', 'BF01_590', 'BF01_730', 'BF01_880']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.products['BF'] = OrderedDict()
        self.products['BFFIT'] = OrderedDict()
        self.products['BFCOV'] = OrderedDict()

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.init_fignames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """

        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']

        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0]),
                                     vcore.vIndex('CCD', vals=self.CCDs),
                                     vcore.vIndex('Quad', vals=self.Quads)])

        #idd = copy.deepcopy(inventoryitem['dd'])
        sidd = self.parse_single_test_gen(jrep, block, testname, inventoryitem)

        # TEST SCPECIFIC
        # TO BE ADDED:
        #   BLOCK, TEST, REPEAT
        #   wavenm, calibrated HK (voltages),
        #   GAIN, EGAIN, ALPHA, BLOOM_ADU, BLOOM_E
        #   REFERENCES TO CURVES

        CHAMBER = sidd.meta['inputs']['CHAMBER']

        CHAMBER_key = CHAMBER[0]
        chamber_v = np.array([CHAMBER_key])
        sidd.addColumn(chamber_v, 'CHAMBERKEY', IndexS, ix=0)

        ogseobj = ogse.Ogse(CHAMBER=CHAMBER)

        wave = sidd.mx['wave'][0, 0]

        wave_v = np.array([ogseobj.get_wavelength(wave)])
        sidd.addColumn(wave_v, 'WAVENM', IndexS, ix=0)

        block_v = np.array([block])
        sidd.addColumn(block_v, 'BLOCK', IndexS, ix=0)

        test_v = np.array([jrep + 1])
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)

        test_v = np.array([session])
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)

        test_v = np.array([testname])
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)

        productspath = os.path.join(inventoryitem['resroot'], 'products')
        BFfitCDP_pick = os.path.join(productspath,
                                     os.path.split(sidd.products['BFfitTABLE_CDP'])[-1])

        BFCDP_pick = os.path.join(productspath, os.path.split(sidd.products['BFTABLE_CDP'])[-1])

        BFfitCDP = files.cPickleRead(BFfitCDP_pick)

        BFCDP = files.cPickleRead(BFCDP_pick)

        bfcdpkeys_v = np.zeros((1), dtype='S50')
        bffitcdpkeys_v = np.zeros((1), dtype='S50')
        bfcovkeys_v = np.zeros((1), dtype='S50')

        BFCOV = sidd.products['COV'].copy()

        bfcdpkey = '%s_%s_%s_%s' % (testname, block, session, jrep + 1)
        bffitcdpkey = '%s_%s_%s_%s' % (testname, block, session, jrep + 1)
        bfcovkey = '%s_%s_%s_%s' % (testname, block, session, jrep + 1)

        self.products['BF'][bfcdpkey] = copy.deepcopy(BFCDP)
        self.products['BFFIT'][bffitcdpkey] = copy.deepcopy(BFfitCDP)
        self.products['BFCOV'][bfcovkey] = copy.deepcopy(BFCOV)

        bfcdpkeys_v[0] = bfcdpkey
        bffitcdpkeys_v[0] = bffitcdpkey
        bfcovkeys_v[0] = bfcovkey

        sidd.addColumn(bfcdpkeys_v, 'BFCDP_KEY', IndexS)
        sidd.addColumn(bffitcdpkeys_v, 'BFFITCDP_KEY', IndexS)
        sidd.addColumn(bfcovkeys_v, 'BFCOV_KEY', IndexS)

        BFfit_df = BFfitCDP['data']['BFFIT']

        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        fwhmx_hwc_v = tmp_v_CQ.copy()
        fwhmy_hwc_v = tmp_v_CQ.copy()

        fwhmx_slope_v = tmp_v_CQ.copy()
        fwhmy_slope_v = tmp_v_CQ.copy()

        ell_v = tmp_v_CQ.copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                ixsel = np.where((BFfit_df['CCD'] == iCCD + 1) & (BFfit_df['Q'] == kQ + 1))

                fwhmx_hwc_v[0, iCCD, kQ] = BFfit_df['FWHMx_HWC'].as_matrix()[ixsel][0]
                fwhmy_hwc_v[0, iCCD, kQ] = BFfit_df['FWHMy_HWC'].as_matrix()[ixsel][0]

                fwhmx_slope_v[0, iCCD, kQ] = BFfit_df['FWHMx_Slope'].as_matrix()[ixsel][0]
                fwhmy_slope_v[0, iCCD, kQ] = BFfit_df['FWHMy_Slope'].as_matrix()[ixsel][0]

                ell_v[0, iCCD, kQ] = BFfit_df['ELL_HWC'].as_matrix()[ixsel][0]

        sidd.addColumn(fwhmx_hwc_v, 'FWHMX_HWC', IndexCQ)
        sidd.addColumn(fwhmy_hwc_v, 'FWHMY_HWC', IndexCQ)
        sidd.addColumn(fwhmx_slope_v, 'FWHMX_SLOPE', IndexCQ)
        sidd.addColumn(fwhmy_slope_v, 'FWHMY_SLOPE', IndexCQ)
        sidd.addColumn(ell_v, 'ELL_HWC', IndexCQ)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _get_GEN_extractor_fromPT(self, colkey):
        """ """

        def _extract_fromPT(PT, block, CCDk, Q):
            ixblock = self.get_ixblock(PT, block)
            column = '%s_%s_Quad%s' % (colkey, CCDk, Q)

            VAL = np.nanmedian(PT[column][ixblock])
            return VAL

        return _extract_fromPT

    def _get_BFvalues(self, _BF, CCD, iQ, orientation):
        """ """

        data = _BF['data']['BF']

        ixsel = np.where((data['CCD'] == CCD) & (data['Q'] == iQ))

        x = []
        y = []

        for ix in ixsel[0]:
            x.append(data['fluence'][ix])
            y.append(data['FWHM%s' % orientation.lower()][ix])

        x = np.array(x)
        y = np.array(y)

        order = np.argsort(x)
        x = x[order].copy()
        y = y[order].copy()

        return x, y

    def _get_BFFITvalues(self, _BFfit, CCD, iQ, orientation):
        """ """

        data = _BFfit['data']['BFFIT']

        ixsel = np.where((data['CCD'] == CCD) & (data['Q'] == iQ))

        x = np.array([0., 2.**16])

        slope = data['FWHM%s_Slope' %
                     orientation.lower()][ixsel[0]].as_matrix()[0]
        y_hwc = data['FWHM%s_HWC' %
                     orientation.lower()][ixsel[0]].as_matrix()[0]

        y = (x - 2.**16 / 2.) * slope * 1.E-4 + y_hwc

        return x, y

    def _get_FWHMzfitsMAP_from_PT(self, testname, orientation='x'):
        """ """

        FWHMzMAP = OrderedDict()
        FWHMzMAP['labelkeys'] = ['E', 'Efit', 'F', 'Ffit',
                                 'G', 'Gfit', 'H', 'Hfit']

        PT = self.ParsedTable[testname]

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                FWHMzMAP[Ckey] = OrderedDict()

                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                CCD = int(CCDk[-1])

                ixblock = np.where(PT['BLOCK'] == block)

                _bfkey = PT['BFCDP_KEY'][ixblock][0]
                _bffitkey = PT['BFFITCDP_KEY'][ixblock][0]

                _BF = self.products['BF'][_bfkey]
                _BFfit = self.products['BFFIT'][_bffitkey]

                _ccd_bfdict = OrderedDict(x=OrderedDict(),
                                          y=OrderedDict())

                for iQ, Q in enumerate(self.Quads):

                    _xRAW, _yRAW = self._get_BFvalues(_BF, CCD, iQ + 1, orientation)

                    _ccd_bfdict['x'][Q] = _xRAW.copy() / 1.e3
                    _ccd_bfdict['y'][Q] = _yRAW.copy()

                    _xFIT, _yFIT = self._get_BFFITvalues(_BFfit, CCD, iQ + 1, orientation)

                    _ccd_bfdict['x']['%sfit' % Q] = _xFIT.copy() / 1.e3
                    _ccd_bfdict['y']['%sfit' % Q] = _yFIT.copy()

                FWHMzMAP[Ckey] = _ccd_bfdict.copy()

        return FWHMzMAP

    def _get_FWHMZ_vs_WAVE(self, orientation='X'):
        """ """

        x = dict()
        y = dict()

        labelkeys = []

        tests = []
        wavenms = []
        for test in self.testnames:
            tests.append(test)
            wavenms.append(self.ParsedTable[test]['WAVENM'][0])

        tests = np.array(tests)
        wavenms = np.array(wavenms)
        ixord = np.argsort(wavenms)
        tests = tests[ixord]
        wavenms = wavenms[ixord]

        for block in self.flight_blocks:

            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' % CCD
                for kQ, Q in enumerate(self.Quads):
                    pkey = '%s_%s_%s' % (block, CCDk, Q)

                    x[pkey] = wavenms.copy()
                    y[pkey] = []

                    for test in tests:
                        PT = self.ParsedTable[test]
                        ixsel = np.where(PT['BLOCK'] == block)
                        colkey = 'FWHM%s_HWC_%s_Quad%s' % (orientation,
                                                           CCDk, Q)
                        y[pkey].append(PT[colkey][ixsel][0])

                    y[pkey] = np.array(y[pkey])

                    labelkeys.append(pkey)

        # Niemi et al. 2015 - for comparison

        x['Niemi'] = wavenms.copy()
        if orientation == 'X':
            y['Niemi'] = 41. * wavenms**(-0.24)
        elif orientation == 'Y':
            y['Niemi'] = 105.9 * wavenms**(-0.38)

        labelkeys.append('Niemi')

        FWHMZdict = dict(x=x, y=y, labelkeys=labelkeys)

        return FWHMZdict

    def _get_FWHMZ_vs_FLU(self, test, orientation='X'):
        """ """

        x = dict()
        y = dict()

        labelkeys = []

        PT = self.ParsedTable[test]

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                CCD = int(CCDk[-1])

                ixblock = np.where(PT['BLOCK'] == block)

                _bfkey = PT['BFCDP_KEY'][ixblock][0]
                _BF = self.products['BF'][_bfkey]

                for iQ, Q in enumerate(self.Quads):

                    _xRAW, _yRAW = self._get_BFvalues(_BF, CCD, iQ + 1, orientation)

                    labelkey = '%s_CCD%i_%s' % (block, CCD, Q)

                    gain = self.cdps['GAIN'][block][CCDk][Q][0]

                    x[labelkey] = (_xRAW * gain / 1.e3).copy()
                    y[labelkey] = _yRAW.copy()

                    labelkeys.append(labelkey)

        x['Niemi'] = np.linspace(1.5e1, 2.e2, 10)  # kilo-electrons
        if orientation == 'X':
            y['Niemi'] = 6.66 + 1.14E-2 * x['Niemi']
        elif orientation == 'Y':
            y['Niemi'] = 6.72 + 1.12E-2 * x['Niemi']

        labelkeys.append('Niemi')

        BFdict = dict(x=x, y=y, labelkeys=labelkeys)

        return BFdict

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['FWHMX_vs_WAVE'] = os.path.join(self.figspath,
                                                  'FWHMX_vs_WAVELENGTH.png')
        self.figs['FWHMY_vs_WAVE'] = os.path.join(self.figspath,
                                                  'FWHMY_vs_WAVELENGTH.png')

        self.figs['FWHMX_vs_FLU'] = os.path.join(self.figspath,
                                                 'FWHMX_vs_FLUENCE.png')
        self.figs['FWHMY_vs_FLU'] = os.path.join(self.figspath,
                                                 'FWHMY_vs_FLUENCE.png')

        # self.figs['SLOPEXY_vs_WAVE'] = os.path.join(self.figspath,
        #                 'SLOPEXY_vs_WAVELENGTH.png')

        for testname in self.testnames:

            self.figs['FWHMX_curves_MAP_%s' % testname] = os.path.join(
                self.figspath, 'FWHMX_Curves_MAP_%s.png' % testname)

            self.figs['SLOPEX_MAP_%s' % testname] = os.path.join(self.figspath,
                                                                 'SLOPEX_MAP_%s.png' % testname)

            self.figs['FWHMX_MAP_%s' % testname] = os.path.join(self.figspath,
                                                                'FWHMX_MAP_%s.png' % testname)

            self.figs['FWHMY_curves_MAP_%s' % testname] = os.path.join(
                self.figspath, 'FWHMY_curves_MAP_%s.png' % testname)

            self.figs['SLOPEY_MAP_%s' % testname] = os.path.join(self.figspath,
                                                                 'SLOPEY_MAP_%s.png' % testname)

            self.figs['FWHMY_MAP_%s' % testname] = os.path.join(self.figspath,
                                                                'FWHMY_MAP_%s.png' % testname)
            self.figs['ELL_MAP_%s' % testname] = os.path.join(self.figspath,
                                                              'ELL_MAP_%s.png' % testname)

    def dump_aggregated_results(self):
        """ """

        # FWHMX vs. WAVE

        BLOCKcolors = cm.rainbow(np.linspace(0, 1, len(self.flight_blocks)))

        FWHMxWAVEdict = self._get_FWHMZ_vs_WAVE(orientation='X')

        FWHMxWAVEkwargs = dict(
            title='FWHM-X vs. Wavelength',
            doLegend=False,
            xlabel='Wavelength [nm]',
            ylabel='FWHMx [um]',
            ylim=[7.5, 9.5],
            figname=self.figs['FWHMX_vs_WAVE'])

        corekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            for iCCD in self.CCDs:
                for kQ in self.Quads:
                    corekwargs['%s_CCD%i_%s' % (block, iCCD, kQ)] = dict(linestyle='-',
                                                                         marker='.', color=jcolor)
        corekwargs['Niemi'] = dict(linestyle='--', marker='o', color='k')

        FWHMxWAVEkwargs['corekwargs'] = corekwargs

        self.plot_XY(FWHMxWAVEdict, **FWHMxWAVEkwargs)

        # FWHMY vs. WAVE

        FWHMyWAVEdict = self._get_FWHMZ_vs_WAVE(orientation='Y')

        FWHMyWAVEkwargs = dict(
            title='FWHM-Y vs. Wavelength',
            doLegend=False,
            xlabel='Wavelength [nm]',
            ylabel='FWHMy [um]',
            ylim=[7.5, 9.5],
            figname=self.figs['FWHMY_vs_WAVE'])

        FWHMyWAVEkwargs['corekwargs'] = corekwargs

        self.plot_XY(FWHMyWAVEdict, **FWHMyWAVEkwargs)

        # FWHMX vs. FLUENCE

        FWHMxFLUdict = self._get_FWHMZ_vs_FLU('BF01', orientation='X')

        FWHMxFLUkwargs = dict(
            title='FWHM-X (800nm) vs. Fluence',
            doLegend=False,
            xlabel='Fluence [ke-]',
            ylabel='FWHMx [um]',
            ylim=[6., 10.],
            figname=self.figs['FWHMX_vs_FLU'])

        FWHMxFLUkwargs['corekwargs'] = corekwargs

        self.plot_XY(FWHMxFLUdict, **FWHMxFLUkwargs)

        # FWHMY vs. FLUENCE

        FWHMyFLUdict = self._get_FWHMZ_vs_FLU('BF01', orientation='Y')

        FWHMyFLUkwargs = dict(
            title='FWHM-Y (800nm) vs. Fluence',
            doLegend=False,
            xlabel='Fluence [ke-]',
            ylabel='FWHMy [um]',
            ylim=[6., 10.],
            figname=self.figs['FWHMY_vs_FLU'])

        FWHMyFLUkwargs['corekwargs'] = corekwargs

        self.plot_XY(FWHMyFLUdict, **FWHMyFLUkwargs)

        # MAPS

        for testname in self.testnames:

            stestname = st.replace(testname, '_', '\_')

            # FWHMX CURVES MAP

            corekwargs_fwhmzmap = dict(E=dict(linestyle='', marker='.', color='r'),
                                       Efit=dict(linestyle='--', marker='', color='r'),
                                       F=dict(linestyle='', marker='.', color='g'),
                                       Ffit=dict(linestyle='--', marker='', color='g'),
                                       G=dict(linestyle='', marker='.', color='b'),
                                       Gfit=dict(linestyle='--', marker='', color='b'),
                                       H=dict(linestyle='', marker='.', color='m'),
                                       Hfit=dict(linestyle='--', marker='', color='m'))

            FWHMX_fits_MAP = self._get_FWHMzfitsMAP_from_PT(testname, orientation='x')

            self.plot_XYMAP(FWHMX_fits_MAP, **dict(
                suptitle='%s: FWHMx Fits' % stestname,
                doLegend=True,
                xlabel='kADU',
                ylim=[6., 10.],
                ylabel='FWHMx [um]',
                corekwargs=corekwargs_fwhmzmap,
                figname=self.figs['FWHMX_curves_MAP_%s' % testname]
                #figname = ''
            ))

            # FWHMY CURVES MAP

            FWHMY_fits_MAP = self._get_FWHMzfitsMAP_from_PT(testname, orientation='y')

            self.plot_XYMAP(FWHMY_fits_MAP, **dict(
                suptitle='%s: FWHMy Fits' % stestname,
                doLegend=True,
                xlabel='kADU',
                ylim=[6., 10.],
                ylabel='FWHMy [um]',
                corekwargs=corekwargs_fwhmzmap,
                figname=self.figs['FWHMY_curves_MAP_%s' % testname]
                #figname = ''
            ))

            # SLOPEX MAP

            SLOPEX_MAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_GEN_extractor_fromPT('FWHMX_SLOPE'))

            self.plot_SimpleMAP(SLOPEX_MAP, **dict(
                suptitle='%s: Slope-x, um/10kADU' % stestname,
                figname=self.figs['SLOPEX_MAP_%s' % testname]))

            # FWHMX MAP

            FWHMX_MAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_GEN_extractor_fromPT('FWHMX_HWC'))

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(FWHMX_MAP, **dict(
                suptitle='%s: FWHM-X at HWC' % stestname,
                figname=self.figs['FWHMX_MAP_%s' % testname]))

            # SLOPEY MAP

            SLOPEY_MAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_GEN_extractor_fromPT('FWHMY_SLOPE'))

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(SLOPEY_MAP, **dict(
                suptitle='%s: Slope-y, um/10kADU' % stestname,
                figname=self.figs['SLOPEY_MAP_%s' % testname]))

            # FWHMY MAP

            FWHMY_MAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_GEN_extractor_fromPT('FWHMY_HWC'))

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(FWHMY_MAP, **dict(
                suptitle='%s: FWHM-Y at HWC' % stestname,
                figname=self.figs['FWHMY_MAP_%s' % testname]))

            # ELL MAP

            ELL_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname],
                                              extractor=self._get_GEN_extractor_fromPT('ELL_HWC'))

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(ELL_MAP, **dict(
                suptitle='%s: ellipticity at HWC' % stestname,
                figname=self.figs['ELL_MAP_%s' % testname]))
