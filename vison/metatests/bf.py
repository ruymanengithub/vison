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

from vison.support import vcal, utils
from vison.datamodel import core as vcore
from vison.datamodel import cdp
from vison.ogse import ogse
from vison.support import files
from vison.analysis import Guyonnet15 as G15

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
        self.init_outcdpnames()

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

    def _get_BFvalues(self, data, CCD, iQ, orientation):
        """ """

        #data = _BF['data']['BF']

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


    def _get_BFvalues_fromA(self, _BFmx, _BFraw, CCD, iQ, orientation, relflu=0.5):
        """ """

        Npix = 51

        #data = _BF['data']['BF']
        Q = self.Quads[iQ-1]

        ixsel = np.where((_BFmx['CCD'] == CCD) & (_BFmx['Q'] == iQ))

        x = []
        y = []

        
        fluences = _BFmx['fluence'].values[ixsel]
        maxflu = np.nanmax(fluences)
        ixflu = np.nanargmin(np.abs(fluences/maxflu-relflu))
        refcol = _BFmx['col'].values[ixsel][ixflu]
        refcorrmap = _BFraw['CCD%i' % CCD][refcol]['av_corrmap'][Q]
        #var = _BFraw['CCD%i' % CCD][selcol]['av_corrmap'][Q]
        refmu = _BFraw['CCD%i' % CCD][refcol]['av_mu'][Q]

        singlepixmap = np.zeros((Npix, Npix), dtype='float32') + 0.0
        singlepixmap[(Npix - 1) / 2, (Npix - 1) / 2] = 1.

        
        Asol_Q, psmooth_Q = G15.solve_for_A_linalg(
            refcorrmap, var=1., mu=refmu, returnAll=True, 
            doplot=False,
            verbose=False)

        for ix in ixsel[0]:
            imu = _BFmx['fluence'][ix]
            x.append(imu)

            kernel_Q = G15.degrade_estatic(singlepixmap*imu, Asol_Q)

            cross_Q = kernel_Q[Npix / 2 - 1:Npix / 2 + 2, 
                    Npix / 2 - 1:Npix / 2 + 2].copy()
            kerQshape = G15.get_cross_shape_rough(
                cross_Q, pitch=12.)

            fwhm = kerQshape['fwhm%s' % orientation.lower()]

            y.append(fwhm)
        
        x = np.array(x)
        y = np.array(y)

        order = np.argsort(x)
        x = x[order].copy()
        y = y[order].copy()

        return x, y

    def _get_BFFITvalues(self, data, CCD, iQ, orientation):
        """ """

        #data = _BFfit['data']['BFFIT']

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

                _BF = self.products['BF'][_bfkey]['data']['BF']
                _BFfit = self.products['BFFIT'][_bffitkey]['data']['BFFIT']

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

    def _get_FWHMZ_vs_FLU(self, test, orientation='X',
                fixFluKernel=False):
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
                _BFmx = self.products['BF'][_bfkey]['data']['BF']
                _BFraw = self.products['BFCOV'][_bfkey]


                if fixFluKernel:
                    
                    for iQ, Q in enumerate(self.Quads):
                        _xRAW, _yRAW = self._get_BFvalues_fromA(_BFmx, _BFraw, CCD, iQ+1,
                            orientation, relflu=0.5)
                else:

                    for iQ, Q in enumerate(self.Quads):

                        _xRAW, _yRAW = self._get_BFvalues(_BFmx, CCD, iQ + 1, orientation)

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

        self.figs['FWHMX_vs_FLU_fixkernel'] = os.path.join(self.figspath,
                                                 'FWHMX_vs_FLUENCE_fixkernel.png')
        self.figs['FWHMY_vs_FLU_fixkernel'] = os.path.join(self.figspath,
                                                 'FWHMY_vs_FLUENCE_fixkernel.png')

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

    def init_outcdpnames(self):

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        self.outcdps['COV'] = 'BFE_COV_matrices.json'
        self.outcdps['G15_PARS'] = 'BFE_GUYO15_PARS.json'
        self.outcdps['G15_GAUSS_PSF'] = 'BFE_GUYO15_GAUSSEQUIV_PSFs.fits'

    def _get_COV_mxs(self):
        """ """

        COV_mxs = OrderedDict()

        for tn in self.testnames:
 
            COV_mxs[tn] = OrderedDict()

            PT = self.ParsedTable[tn]

            for jY in range(self.NSLICES_FPA):
                for iX in range(self.NCOLS_FPA):
                    Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                    COV_mxs[tn][Ckey] = OrderedDict()

                    locator = self.fpa.FPA_MAP[Ckey]
                    block = locator[0]
                    CCDk = locator[1]

                    COV_mxs[tn][Ckey]['block'] = block
                    COV_mxs[tn][Ckey]['CCD'] = CCDk

                    ixblock = np.where(PT['BLOCK'] == block)

                    _covkey = PT['BFCOV_KEY'][ixblock][0]

                    _COV = self.products['BFCOV'][_covkey][CCDk]

                    Ncols = len(_COV.keys())

                    for ic in range(Ncols):
                        colk = 'col%03i' % (ic+1,)

                        COV_mxs[tn][Ckey][colk] = OrderedDict()

                        COV_mxs[tn][Ckey][colk]['Npairs'] = _COV[colk]['Npairs']

                        COV_mxs[tn][Ckey][colk]['mu'] = OrderedDict()
                        COV_mxs[tn][Ckey][colk]['var'] = OrderedDict()
                        COV_mxs[tn][Ckey][colk]['corrmap'] = OrderedDict()

                        for Q in self.Quads:

                            COV_mxs[tn][Ckey][colk]['mu'][Q] =\
                                 float(_COV[colk]['av_mu'][Q])
                            COV_mxs[tn][Ckey][colk]['var'][Q] =\
                                 float(_COV[colk]['av_var'][Q])
                            COV_mxs[tn][Ckey][colk]['corrmap'][Q] = \
                                _COV[colk]['av_corrmap'][Q].tolist()

        return COV_mxs                       

    def _get_G15P_mxs(self):
        """ """

        G15P_mxs = OrderedDict()

        for tn in self.testnames:
 
            G15P_mxs[tn] = OrderedDict()

            PT = self.ParsedTable[tn]

            metrics = ['fluence','psmooth','e_psmooth','Asol']
            filldefault = OrderedDict()
            for metric in metrics:
                filldefault[metric] = OrderedDict()
            filldefault['fluence']['E'] = np.nan
            filldefault['psmooth']['E'] = [np.nan,np.nan]
            filldefault['e_psmooth']['E'] = [np.nan, np.nan]
            filldefault['Asol']['E'] = [[[np.nan]]]
            for metric in metrics:
                for Q in ['F','G','H']:
                    filldefault[metric][Q] = copy.deepcopy(filldefault[metric]['E'])

            for jY in range(self.NSLICES_FPA):
                for iX in range(self.NCOLS_FPA):
                    Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                    G15P_mxs[tn][Ckey] = OrderedDict()

                    locator = self.fpa.FPA_MAP[Ckey]
                    block = locator[0]
                    CCDk = locator[1]

                    G15P_mxs[tn][Ckey]['block'] = block
                    G15P_mxs[tn][Ckey]['CCD'] = CCDk

                    ixblock = np.where(PT['BLOCK'] == block)

                    _bfkey = PT['BFCDP_KEY'][ixblock][0]
                    _, session, srep = st.split(st.replace(_bfkey,'%s_' % tn,''),'_')

                    rep = int(srep)

                    bf = self.inventory[block][tn][rep-1]['dd'].products['BF']

                    cols = bf['ulabels']


                    for ic, colk in enumerate(cols):
                        colk = 'col%03i' % (ic+1,)

                        G15P_mxs[tn][Ckey][colk] = filldefault.copy()

                        for Q in self.Quads:

                            if len(bf[CCDk][Q][colk])>0:

                                G15P_mxs[tn][Ckey][colk]['fluence'][Q] =\
                                    float(bf[CCDk][Q][colk]['fluence'])
                                
                                G15P_mxs[tn][Ckey][colk]['psmooth'][Q] =\
                                     bf[CCDk][Q][colk]['psmooth'][0].tolist()
                                G15P_mxs[tn][Ckey][colk]['e_psmooth'][Q] =\
                                     bf[CCDk][Q][colk]['psmooth'][1].tolist()
                                G15P_mxs[tn][Ckey][colk]['Asol'][Q] =\
                                     bf[CCDk][Q][colk]['Asol'].tolist()


        return G15P_mxs                       



    def _get_g15_gpsf_cdp(self, inCDP_header=None):
        """ """
        
        CDP_header = OrderedDict()
        if CDP_header is not None:
            CDP_header.update(inCDP_header)

        cdpname = self.outcdps['G15_GAUSS_PSF']
        path = self.cdpspath

        gpsf_cdp = cdp.FitsTables_CDP()
        gpsf_cdp.rootname = os.path.splitext(cdpname)[0]
        gpsf_cdp.path = path

        gpsf_meta = OrderedDict()

        gpsf_data = OrderedDict()

        cols = ['col', 'CCD', 'Q', 'fluence', 'FWHMx', 'FWHMy', 'e']

        for tn in self.testnames:

            PT = self.ParsedTable[tn]

            for block in self.flight_blocks:
                ixblock = np.where(PT['BLOCK'] == block)
                _bfkey = PT['BFCDP_KEY'][ixblock][0]

                BFdf = self.products['BF'][_bfkey]['data']['BF'].copy()

                rawd=BFdf.to_dict(orient='list')
                dd = OrderedDict()
                for col in cols: dd[col] = np.array(rawd[col]).copy()
                

                gpsf_data['%s_%s' % (tn, block)] = dd.copy()



        CDP_header = self.FITSify_CDP_header(CDP_header)

        gpsf_cdp.ingest_inputs(data=gpsf_data.copy(),
            meta=gpsf_meta.copy(),
            header=CDP_header.copy())
        gpsf_cdp.init_HL_and_fillAll()

        gpsf_cdp.hdulist[0].header.insert(CDP_header.keys()[0],
            ('title','BF0X: GAUSS PSF EQUIVALENTS TO BFE-G+15'))

        return gpsf_cdp



    def dump_aggregated_results(self):
        """ """

        
        if self.report is not None:
            self.report.add_Section(keyword='dump', 
                Title='Aggregated Results', level=0)
            
            self.add_DataAlbaran2Report()

        doAll = True

        doFWHMvWAVE = doAll
        doFWHMvFLU = doAll
        doFWHMvFLUalt = doAll
        doTestMAPs = doAll
        doELLMAP = doAll
        doCDPs = doAll
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()



        # FWHM? vs. WAVE

        BLOCKcolors = cm.rainbow(np.linspace(0, 1, len(self.flight_blocks)))

        comcorekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            for iCCD in self.CCDs:
                for kQ in self.Quads:
                    comcorekwargs['%s_CCD%i_%s' % (block, iCCD, kQ)] = dict(linestyle='-',
                        marker='.', color=jcolor)
        comcorekwargs['Niemi'] = dict(linestyle='--', marker='o', color='k')


        if doFWHMvWAVE:


            for dim in ['x','y']:

                figkey1 = 'FWHM%s_vs_WAVE' % dim.upper()
                figname1 = self.figs[figkey1]


                FWHMzWAVEdict = self._get_FWHMZ_vs_WAVE(orientation=dim.upper())

                FWHMzWAVEkwargs = dict(
                    title='FWHM-%s vs. Wavelength' % dim.upper(),
                    doLegend=False,
                    xlabel='Wavelength [nm]',
                    ylabel='FWHM%s [um]' % dim,
                    ylim=[7.5, 12.5],
                    figname=figname1)


                FWHMzWAVEkwargs['corekwargs'] = comcorekwargs

                self.plot_XY(FWHMzWAVEdict, **FWHMzWAVEkwargs)

                if self.report is not None:
                    self.addFigure2Report(figname1, 
                        figkey=figkey1, 
                        caption='FWHM%s vs. Wavelength' % dim, 
                        texfraction=0.7)

        if doFWHMvFLU:

            # FWHM? vs. FLUENCE

            for dim in ['x','y']:

                figkey2 = 'FWHM%s_vs_FLU' % dim.upper()
                figname2 = self.figs[figkey2]


                FWHMzFLUdict = self._get_FWHMZ_vs_FLU('BF01', orientation=dim.upper())

                FWHMzFLUkwargs = dict(
                    title='FWHM-%s (800nm) vs. Fluence' % dim.upper(),
                    doLegend=False,
                    xlabel='Fluence [ke-]',
                    ylabel='FWHM%s [um]' % dim,
                    ylim=[6., 12.5],
                    figname=figname2)

                FWHMzFLUkwargs['corekwargs'] = comcorekwargs

                self.plot_XY(FWHMzFLUdict, **FWHMzFLUkwargs)

                if self.report is not None:
                    self.addFigure2Report(figname2, 
                        figkey=figkey2, 
                        caption='FWHM%s vs. Wavelength' % dim, 
                        texfraction=0.7)


        if doFWHMvFLUalt:

            # FWHM? vs. FLUENCE
            # Using a fixed-fluence kernel this time

            for dim in ['x','y']:

                figkey21 = 'FWHM%s_vs_FLU_fixkernel' % dim.upper()
                figname21 = self.figs[figkey21]


                FWHMzFLUFKdict = self._get_FWHMZ_vs_FLU('BF01', orientation=dim.upper(),
                    fixFluKernel=True)

                FWHMzFLUFKkwargs = dict(
                    title='FWHM-%s (800nm) vs. Fluence, Fixed Kernel' % dim.upper(),
                    doLegend=False,
                    xlabel='Fluence [ke-]',
                    ylabel='FWHM%s [um]' % dim,
                    ylim=[6., 12.5],
                    figname=figname21)

                FWHMzFLUFKkwargs['corekwargs'] = comcorekwargs

                self.plot_XY(FWHMzFLUFKdict, **FWHMzFLUFKkwargs)

                if self.report is not None:
                    self.addFigure2Report(figname21, 
                        figkey=figkey21, 
                        caption='FWHM%s vs. Wavelength, using the kernel at mid-fluence.' % dim, 
                        texfraction=0.7)


        # MAPS

        if doTestMAPs:

            for testname in self.testnames:

                stestname = st.replace(testname, '_', '\_')

                for dim in ['x','y']:

                    # FWHMz CURVES MAP

                    figkey3 = 'FWHM%s_curves_MAP_%s' \
                                % (dim.upper(),testname)

                    figname3 = self.figs[figkey3]

                    corekwargs_fwhmzmap = dict(E=dict(linestyle='', marker='.', color='r'),
                                               Efit=dict(linestyle='--', marker='', color='r'),
                                               F=dict(linestyle='', marker='.', color='g'),
                                               Ffit=dict(linestyle='--', marker='', color='g'),
                                               G=dict(linestyle='', marker='.', color='b'),
                                               Gfit=dict(linestyle='--', marker='', color='b'),
                                               H=dict(linestyle='', marker='.', color='m'),
                                               Hfit=dict(linestyle='--', marker='', color='m'))

                    FWHMZ_fits_MAP = self._get_FWHMzfitsMAP_from_PT(testname, orientation=dim)

                    self.plot_XYMAP(FWHMZ_fits_MAP, **dict(
                        suptitle='%s: FWHM%s Fits' % (stestname,dim),
                        doLegend=True,
                        xlabel='kADU',
                        ylim=[6., 12.5],
                        ylabel='FWH%s [um]' % dim,
                        corekwargs=corekwargs_fwhmzmap,
                        figname=figname3
                        #figname = ''
                    ))

                    if self.report is not None:
                        self.addFigure2Report(figname3, 
                            figkey=figkey3, 
                            caption='%s' % stestname, 
                            texfraction=0.7)


                    # SLOPEZ MAP

                    figkey4 = 'SLOPE%s_MAP_%s' % (dim.upper(),testname)
                    figname4 = self.figs[figkey4]

                    SLOPEz_MAP = self.get_FPAMAP_from_PT(
                        self.ParsedTable[testname],
                        extractor=self._get_GEN_extractor_fromPT('FWHM%s_SLOPE' % \
                            dim.upper()))

                    self.plot_SimpleMAP(SLOPEz_MAP, **dict(
                        suptitle='%s: Slope-%s, um/10kADU' % (stestname,dim),
                        figname=figname4))

                    if self.report is not None:
                        self.addFigure2Report(figname4, 
                            figkey=figkey4, 
                            caption='%s' % stestname, 
                            texfraction=0.7)

                    # FWHMZ MAP

                    figkey5 = 'FWHM%s_MAP_%s' % (dim.upper(),testname)
                    figname5 = self.figs[figkey5]

                    FWHMZ_MAP = self.get_FPAMAP_from_PT(
                        self.ParsedTable[testname],
                        extractor=self._get_GEN_extractor_fromPT('FWHM%s_HWC' % \
                            dim.upper()))

                    stestname = st.replace(testname, '_', '\_')
                    self.plot_SimpleMAP(FWHMZ_MAP, **dict(
                        suptitle='%s: FWHM-%s at HWC' % (stestname,dim.upper()),
                        figname=figname5))

                    if self.report is not None:
                        self.addFigure2Report(figname5, 
                            figkey=figkey5, 
                            caption='%s' % stestname, 
                            texfraction=0.7)
            
        if doELLMAP:

            # ELL MAP

            figkey6 = 'ELL_MAP_%s' % testname
            figname6 = self.figs[figkey6]

            ELL_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname],
                            extractor=self._get_GEN_extractor_fromPT('ELL_HWC'))

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(ELL_MAP, **dict(
                suptitle='%s: ellipticity at HWC' % stestname,
                figname=figname6))

            if self.report is not None:
                self.addFigure2Report(figname6,
                    figkey=figkey6,
                    caption='%s: Ellipticity of the Gaussian-equivalent kernel.' % stestname, 
                    texfraction=0.7)


        # CDPs

        if doCDPs:

            # Covariance Matrices
            # json

            
            COVs_mx = self._get_COV_mxs()

            cov_header = OrderedDict()
            cov_header['title'] = 'Covariance Matrices'
            cov_header['test'] = 'BF0X'
            cov_header.update(CDP_header)

            cov_meta = OrderedDict()
            cov_meta['structure'] = dict(
                    TEST=dict(
                        CCDID=dict(
                            block="BLOCK ID",
                            CCD="ROE CCD",
                            colXXX=dict(
                                Npairs="Nr. of pairs of images",
                                mu="dict(Quad), median fluence, [ADU]",
                                var="dict(Quad), variances [ADU^2]",
                                corrmap="dict(Quad), correlations, [ADIM.]"
                                )
                            )
                        )
                    )


            cov_cdp = cdp.Json_CDP(rootname=self.outcdps['COV'],
                path=self.cdpspath)

            cov_cdp.ingest_inputs(data=COVs_mx,
                header=cov_header,
                meta=cov_meta)
            cov_cdp.savehardcopy()

            

            # G15 Parameters
            # json

            G15P_mx = self._get_G15P_mxs()

            g15p_header = OrderedDict()
            g15p_header['title'] = 'Guyonnet+15 Parameters'
            g15p_header['test'] = 'BF0X'
            g15p_header.update(CDP_header)

            g15p_meta = OrderedDict()
            g15p_meta['structure'] = dict(
                    TEST=dict(
                        CCDID=dict(
                            block="BLOCK ID",
                            CCD="ROE CCD",
                            colXXX=dict(
                                fluence="dict(Q), Fluence, [ADU]",
                                psmooth="dict(Q), (p0, p1), [p0,p1]=[L, 1/L]",
                                e_psmooth="dict(Q), (unc(p0), unc(p1))",
                                Asol="dict(Q), matrix of a_{i,j}^X coefficients, [ADU^{-1}]."
                                )
                            )
                        )
                    )


            g15p_cdp = cdp.Json_CDP(rootname=self.outcdps['G15_PARS'],
                path=self.cdpspath)

            g15p_cdp.ingest_inputs(data=G15P_mx,
                header=g15p_header,
                meta=g15p_meta)
            g15p_cdp.savehardcopy()


            # G15 G-PSF Equivalent Parameters
            # FITS table

            g15_gpsf_cdp = self._get_g15_gpsf_cdp(inCDP_header=CDP_header)

            for i in range(2,len(g15_gpsf_cdp.hdulist)):
                extname = g15_gpsf_cdp.hdulist[i].header['EXTNAME']
                block = st.split(extname,'_')[-1]
                tn = st.join(st.split(extname,'_')[:-1],'_')

                g15_gpsf_cdp.hdulist[i].header['TEST'] = tn
                g15_gpsf_cdp.hdulist[i].header['BLOCK'] = block

                for CCD in self.CCDs:
                    CCDk = 'CCD%i' % CCD
                    CCDID = self.fpa.get_Ckey_from_BlockCCD(block, CCD)
                    g15_gpsf_cdp.hdulist[i].header[CCDk] = CCDID

                g15_gpsf_cdp.hdulist[0].header['EXT_%i' % (i+1,)] = extname

            g15_gpsf_cdp.savehardcopy()

