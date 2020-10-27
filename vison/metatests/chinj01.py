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
import pandas as pd

from vison.datamodel import cdp
from vison.support import files
from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal, utils
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.inject import lib as ilib


import matplotlib.cm as cm
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
    'HK_CCD1_TEMP_T',
    'HK_CCD2_TEMP_T',
    'HK_CCD3_TEMP_T',
    'HK_CCD1_TEMP_B',
    'HK_CCD2_TEMP_B',
    'HK_CCD3_TEMP_B',
    'HK_CCD1_OD_T',
    'HK_CCD2_OD_T',
    'HK_CCD3_OD_T',
    'HK_CCD1_OD_B',
    'HK_CCD2_OD_B',
    'HK_CCD3_OD_B',
    'HK_COMM_RD_T',
    'HK_COMM_RD_B',
    'HK_CCD1_IG1_T',
    'HK_CCD2_IG1_T',
    'HK_CCD3_IG1_T',
    'HK_CCD1_IG1_B',
    'HK_CCD2_IG1_B',
    'HK_CCD3_IG1_B',
    'HK_COMM_IG2_T',
    'HK_COMM_IG2_B',
    'HK_FPGA_BIAS_ID2',
    'HK_VID_PCB_TEMP_T',
    'HK_VID_PCB_TEMP_B',
    'HK_RPSU_TEMP1',
    'HK_FPGA_PCB_TEMP_T',
    'HK_FPGA_PCB_TEMP_B',
    'HK_RPSU_TEMP_2',
    'HK_RPSU_28V_PRI_I',
    'chk_NPIXOFF',
    'chk_NPIXSAT',
    'offset_pre',
    'offset_ove',
    'std_pre',
    'std_ove']


class MetaChinj01(MetaCal):
    """ """

    def __init__(self, **kwargs):
        """ """

        super(MetaChinj01, self).__init__(**kwargs)

        self.testnames = ['CHINJ01']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['METAFIT'] = OrderedDict()
        self.products['VERPROFILES'] = OrderedDict()
        self.products['HORPROFILES'] = OrderedDict()

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
        #   OFFSETS: pre, img, ove
        #   RON: pre, img, ove
        #   REFERENCES TO PROFILES

        CHAMBER = sidd.meta['inputs']['CHAMBER']
        CHAMBER_key = CHAMBER[0]
        chamber_v = np.array([CHAMBER_key])
        sidd.addColumn(chamber_v, 'CHAMBERKEY', IndexS, ix=0)

        block_v = np.array([block])
        sidd.addColumn(block_v, 'BLOCK', IndexS, ix=0)

        test_v = np.array([jrep + 1])
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)

        test_v = np.array([session])
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)

        test_v = np.array([testname])
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)

        productspath = os.path.join(inventoryitem['resroot'], 'products')

        metafitcdp_pick = os.path.join(productspath,
                                       os.path.split(sidd.products['METAFIT_CDP'])[-1])
        metafitcdp = files.cPickleRead(metafitcdp_pick)
        metafit = copy.deepcopy(metafitcdp['data']['ANALYSIS'])

        metafitkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['METAFIT'][metafitkey] = copy.deepcopy(metafit)
        metafitkey_v = np.array([metafitkey])
        sidd.addColumn(metafitkey_v, 'METAFIT', IndexS, ix=0)

        metacdp_pick = os.path.join(productspath, os.path.split(
            sidd.products['META_CDP'])[-1])  # change to META_CDP
        metacdp = files.cPickleRead(metacdp_pick)
        meta = metacdp['data']['ANALYSIS']  # this is a pandas DataFrame


        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        bgd_adu_v = tmp_v_CQ.copy()
        ig1_thresh_v = tmp_v_CQ.copy()
        ig1_notch_v = tmp_v_CQ.copy()
        slope_v = tmp_v_CQ.copy()
        n_adu_v = tmp_v_CQ.copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                ixloc = np.where((meta['CCD'] == iCCD + 1) & (meta['Q'] == kQ + 1))

                bgd_adu_v[0, iCCD, kQ] = meta['BGD_ADU'][ixloc[0][0]]
                ig1_thresh_v[0, iCCD, kQ] = meta['IG1_THRESH'][ixloc[0][0]]
                ig1_notch_v[0, iCCD, kQ] = meta['IG1_NOTCH'][ixloc[0][0]]
                slope_v[0, iCCD, kQ] = meta['S'][ixloc[0][0]]
                n_adu_v[0, iCCD, kQ] = meta['N_ADU'][ixloc[0][0]]

        sidd.addColumn(bgd_adu_v, 'FIT_BGD_ADU', IndexCQ)
        sidd.addColumn(ig1_thresh_v, 'FIT_IG1_THRESH', IndexCQ)
        sidd.addColumn(ig1_notch_v, 'FIT_IG1_NOTCH', IndexCQ)
        sidd.addColumn(slope_v, 'FIT_SLOPE', IndexCQ)
        sidd.addColumn(n_adu_v, 'FIT_N_ADU', IndexCQ)

        # charge injection profiles

        verprofspick = os.path.join(productspath,
            os.path.split(sidd.products['PROFS_ALCOL'])[-1])
        verprofs = files.cPickleRead(verprofspick)

        vprofkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)

        self.products['VERPROFILES'][vprofkey] = verprofs.copy()

        vprofskeys_v = np.zeros((1),dtype='U50')

        vprofskeys_v[0] = vprofkey

        sidd.addColumn(vprofskeys_v, 'VERPROFS_KEY', IndexS)

        horprofspick = os.path.join(productspath,
            os.path.split(sidd.products['PROFS_ALROW'])[-1])
        horprofs = files.cPickleRead(horprofspick)        

        hprofkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)

        self.products['HORPROFILES'][hprofkey] = horprofs.copy()

        hprofskeys_v = np.zeros((1),dtype='U50')

        hprofskeys_v[0] = hprofkey

        sidd.addColumn(hprofskeys_v, 'HORPROFS_KEY', IndexS)


        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _get_extractor_NOTCH_fromPT(self, units):
        """ """

        def _extract_NOTCH_fromPT(PT, block, CCDk, Q):

            ixblock = self.get_ixblock(PT, block)
            column = 'FIT_N_ADU_%s_Quad%s' % (CCDk, Q)

            if units == 'ADU':
                unitsConvFactor = 1
            elif units == 'E':
                unitsConvFactor = self.cdps['GAIN'][block][CCDk][Q][0]

            Notch = np.nanmedian(PT[column][ixblock]) * unitsConvFactor
            return Notch

        return _extract_NOTCH_fromPT

    def _get_injcurve(self, _chfitdf, ixCCD, ixQ, IG1raw, gain):
        """ """
        ixsel = np.where((_chfitdf['CCD'] == ixCCD) & (_chfitdf['Q'] == ixQ))

        pars = ['BGD', 'K', 'XT', 'XN', 'A', 'N']
        trans = dict(BGD='b', K='k', XT='xt', XN='xN', A='a', N='N')

        parsdict = dict()
        for par in pars:
            parsdict[trans[par]] = _chfitdf[par].values[ixsel][0]

        parsdict['IG1'] = IG1raw.copy()

        inj = ilib.f_Inj_vs_IG1_ReLU(**parsdict) * 2.**16  # ADU

        inj_kel = inj * gain / 1.E3

        return inj_kel

    def _get_CHIG1_MAP_from_PT(self, kind='CAL'):
        """ """

        CHIG1MAP = OrderedDict()
        CHIG1MAP['labelkeys'] = self.Quads

        PT = self.ParsedTable['CHINJ01']
        column = 'METAFIT'

        IG1s = [2.5, 6.75]
        dIG1 = 0.05

        NIG1 = (IG1s[1] - IG1s[0]) / dIG1 + 1
        IG1raw = np.arange(NIG1) * dIG1 + IG1s[0]

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                CHIG1MAP[Ckey] = OrderedDict()

                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]

                jCCD = int(CCDk[-1])

                ixblock = np.where(PT['BLOCK'] == block)

                if len(ixblock[0]) == 0:
                    CHIG1MAP[Ckey] = OrderedDict(x=OrderedDict(),
                                                 y=OrderedDict())
                    for Q in self.Quads:
                        CHIG1MAP[Ckey]['x'][Q] = []
                        CHIG1MAP[Ckey]['y'][Q] = []
                    continue

                _chkey = PT[column][ixblock][0]

                _chfitdf = self.products['METAFIT'][_chkey]

                _ccd_chfitdict = OrderedDict(x=OrderedDict(),
                                             y=OrderedDict())

                for kQ, Q in enumerate(self.Quads):

                    roeVCal = self.roeVCals[block]

                    IG1cal = roeVCal.fcal_HK(IG1raw, 'IG1', jCCD, Q)

                    gain = self.cdps['GAIN'][block][CCDk][Q][0]

                    inj_kel = self._get_injcurve(_chfitdf, jCCD, kQ + 1, IG1raw, gain)

                    if kind == 'CAL':
                        _IG1 = IG1cal.copy()
                    elif kind == 'RAW':
                        _IG1 = IG1raw.copy()

                    _ccd_chfitdict['x'][Q] = _IG1.copy()
                    _ccd_chfitdict['y'][Q] = inj_kel.copy()

                CHIG1MAP[Ckey] = _ccd_chfitdict.copy()

        return CHIG1MAP

    def _get_XYdict_INJ(self, kind='CAL'):

        x = dict()
        y = dict()

        PT = self.ParsedTable['CHINJ01']
        column = 'METAFIT'

        IG1s = [2.5, 6.75]
        dIG1 = 0.05

        NIG1 = (IG1s[1] - IG1s[0]) / dIG1 + 1
        IG1raw = np.arange(NIG1) * dIG1 + IG1s[0]

        labelkeys = []

        for block in self.flight_blocks:
            ixblock = np.where(PT['BLOCK'] == block)
            ch_key = PT[column][ixblock][0]
            chfitdf = self.products['METAFIT'][ch_key]

            for iCCD, CCD in enumerate(self.CCDs):

                CCDk = 'CCD%i' % CCD

                for kQ, Q in enumerate(self.Quads):

                    roeVCal = self.roeVCals[block]

                    IG1cal = roeVCal.fcal_HK(IG1raw, 'IG1', iCCD + 1, Q)

                    gain = self.cdps['GAIN'][block][CCDk][Q][0]

                    if kind == 'CAL':
                        _IG1 = IG1cal.copy()
                    elif kind == 'RAW':
                        _IG1 = IG1raw.copy()

                    pkey = '%s_%s_%s' % (block, CCDk, Q)

                    inj_kel = self._get_injcurve(chfitdf, iCCD + 1, kQ + 1, IG1raw, gain)

                    x[pkey] = _IG1.copy()
                    y[pkey] = inj_kel.copy()
                    labelkeys.append(pkey)

        CHdict = dict(x=x, y=y, labelkeys=labelkeys)

        return CHdict

    def _extract_INJCURVES_PAR_fromPT(self,PT,block,CCDk,Q):
        """ """
        ixblock = self.get_ixblock(PT,block)
        column = 'METAFIT'
        ch_key = PT[column][ixblock][0]
        chfitdf = self.products['METAFIT'][ch_key]

        ixCCD = ['CCD1','CCD2','CCD3'].index(CCDk)+1
        ixQ = ['E','F','G','H'].index(Q)+1

        ixsel = np.where((chfitdf['CCD'] == ixCCD) & (chfitdf['Q'] == ixQ))

        pars = ['BGD', 'K', 'XT', 'XN', 'A', 'N']
        trans = dict(BGD='b', K='k', XT='xt', XN='xN', A='a', N='N')

        parsdict = dict()
        for par in pars:
            parsdict[trans[par]] = '%.3e' % chfitdf[par].values[ixsel][0]

        return parsdict

    def _get_XYdict_PROFS(self,proftype, IG1=4.5, Quads=None, doNorm=False, xrangeNorm=None):
        """ """

        if Quads is None:
            Quads = self.Quads

        x = dict()
        y = dict()

        labelkeys = []

        PT = self.ParsedTable['CHINJ01']

        profcol = '%sPROFS_KEY' % proftype.upper()
        prodkey = '%sPROFILES' % proftype.upper()

        for block in self.flight_blocks:

            ixsel = np.where(PT['BLOCK'] == block)
            prof_key = PT[profcol][ixsel][0]

            i_Prof = self.products[prodkey][prof_key].copy()

            IG1key = 'IG1_%.2fV' % IG1

            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' % CCD

                for kQ, Q in enumerate(Quads):

                    pkey = '%s_%s_%s' % (block, CCDk, Q)

                    _pcq = i_Prof['data'][CCDk][Q].copy()

                    _x = _pcq['x'][IG1key].copy()
                    _y = _pcq['y'][IG1key].copy()

                    x[pkey] = _x
                    if doNorm:

                        if xrangeNorm is not None:
                            norm = np.nanmedian(_y[xrangeNorm[0]:xrangeNorm[1]])
                        else:
                            norm = np.nanmedian(_y)

                        y[pkey] = _y / norm

                    labelkeys.append(pkey)

        Pdict = dict(x=x,y=y,labelkeys=labelkeys)

        return Pdict

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['NOTCH_ADU_MAP'] = os.path.join(self.figspath,
                                                  'NOTCH_ADU_MAP.png')

        self.figs['NOTCH_ELE_MAP'] = os.path.join(self.figspath,
                                                  'NOTCH_ELE_MAP.png')

        self.figs['CHINJ01_curves_IG1_RAW'] = os.path.join(self.figspath,
                                                           'CHINJ01_CURVES_IG1_RAW.png')

        self.figs['CHINJ01_curves_IG1_CAL'] = os.path.join(self.figspath,
                                                           'CHINJ01_CURVES_IG1_CAL.png')

        self.figs['CHINJ01_curves_MAP_IG1_CAL'] = os.path.join(self.figspath,
                                                               'CHINJ01_CURVES_MAP_IG1_CAL.png')

        for proftype in ['ver','hor']:
            for ccdhalf in ['top','bot']:
                figkey = 'PROFS_%s_%s' % (proftype.upper(),ccdhalf.upper())
                self.figs[figkey] = os.path.join(self.figspath,
                        'CHINJ01_%s_%s_PROFILES.png' % \
                            (proftype.upper(),ccdhalf.upper()))

        for ccdhalf in ['top','bot']:
            figkey = 'PROFS_ver_%s_ZOOM' % (ccdhalf.upper(),)
            self.figs[figkey] = os.path.join(self.figspath,
                    'CHINJ01_ver_%s_ZOOM_PROFILES.png' % \
                        (ccdhalf.upper()),)

    def init_outcdpnames(self):

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        self.outcdps['INJCURVES'] = 'CHINJ01_INJCURVES_PAR.json'
        self.outcdps['INJPROF_XLSX_HOR'] = 'CHINJ01_INJPROFILES_HOR.xlsx'
        self.outcdps['INJPROF_XLSX_VER'] = 'CHINJ01_INJPROFILES_VER.xlsx'
        self.outcdps['INJPROF_FITS_HOR'] = 'CHINJ01_INJPROFILES_HOR.fits'
        self.outcdps['INJPROF_FITS_VER'] = 'CHINJ01_INJPROFILES_VER.fits'


    def _extract_NUNHOR_fromPT(self, PT, block, CCDk, Q):
        """ """

        IG1 = 4.5

        ixblock = self.get_ixblock(PT, block)

        profcol = 'HORPROFS_KEY'
        prodkey = 'HORPROFILES'

        prof_key = PT[profcol][ixblock][0]

        i_Prof = self.products[prodkey][prof_key].copy()

        IG1key = 'IG1_%.2fV' % IG1

        _pcq = i_Prof['data'][CCDk][Q].copy()
        _y = _pcq['y'][IG1key].copy()

        return np.nanstd(_y)/np.nanmean(_y)*100.

    def _get_injprof_dfdict(self, direction, pandice=False):
        """ """

        injprofs = OrderedDict()

        Quads = self.Quads
        PT = self.ParsedTable['CHINJ01']
        profcol = '{}PROFS_KEY'.format(direction.upper())
        prodkey = '{}PROFILES'.format(direction.upper())


        for ib, block in enumerate(self.flight_blocks):
            injprofs[block] = OrderedDict()

            ixsel = np.where(PT['BLOCK'] == block)
            prof_key = PT[profcol][ixsel][0]

            i_Prof = self.products[prodkey][prof_key].copy()

            if ib==0:
                rawIG1keys = list(i_Prof['data']['CCD1']['E']['x'].keys())
                IG1values = [float(item.replace('IG1_','').replace('V','')) for item in rawIG1keys]
                _order = np.argsort(IG1values)
                IG1keys = np.array(rawIG1keys)[_order].tolist()
                IG1values = np.array(IG1values)[_order].tolist()
            
            for IG1key in IG1keys:
            
                for iCCD, CCD in enumerate(self.CCDs):
                    CCDk = 'CCD%i' % CCD

                    Ckey = self.fpa.get_Ckey_from_BlockCCD(block, CCD)

                    for kQ, Q in enumerate(Quads):

                        _pcq = i_Prof['data'][CCDk][Q].copy()

                        _x = _pcq['x'][IG1key].copy()
                        _y = _pcq['y'][IG1key].copy()

                        #_y /= np.nanmedian(_y)

                        if iCCD==0 and kQ==0:
                            injprofs[block]['pixel'] = _x.copy()
                        injprofs[block]['%s_%s_%s' % (Ckey,Q,IG1key)] = _y.copy() 

        if pandice:
            for block in self.flight_blocks:
                injprofs[block] = pd.DataFrame.from_dict(injprofs[block])

        return injprofs, IG1values


    def get_injprof_xlsx_cdp(self, direction, inCDP_header=None):
        """ """

        CDP_header = OrderedDict()
        if CDP_header is not None:
            CDP_header.update(inCDP_header)

        cdpname = self.outcdps['INJPROF_XLSX_%s' % direction.upper()]
        path = self.cdpspath

        injprof_cdp = cdp.Tables_CDP()
        injprof_cdp.rootname = os.path.splitext(cdpname)[0]
        injprof_cdp.path = path

        injprofs_meta = OrderedDict()

        injprofs, IG1values = self._get_injprof_dfdict(direction, pandice=True)        

        injprofs_meta['IG1'] = IG1values.__repr__()
        #injprofs_meta['norm'] = 'median'

        injprof_cdp.ingest_inputs(data=injprofs.copy(),
            meta=injprofs_meta.copy(),
            header=CDP_header.copy())
        injprof_cdp.init_wb_and_fillAll(
            header_title='CHINJ01: INJPROFS-%s' % direction.upper())

        return injprof_cdp

    def get_injprof_fits_cdp(self, direction, inCDP_header=None):
        """ """

        CDP_header = OrderedDict()
        if inCDP_header is not None:
            CDP_header.update(inCDP_header)

        cdpname = self.outcdps['INJPROF_FITS_%s' % direction.upper()]
        path = self.cdpspath

        injprof_cdp = cdp.FitsTables_CDP()
        injprof_cdp.rootname = os.path.splitext(cdpname)[0]
        injprof_cdp.path = path

        injprofs_meta = OrderedDict()

        injprofs, IG1values = self._get_injprof_dfdict(direction, pandice=False)

        injprofs_meta['IG1'] = IG1values.__repr__()
        #injprofs_meta['norm'] = 'median'

        CDP_header = self.FITSify_CDP_header(CDP_header)

        injprof_cdp.ingest_inputs(data=injprofs.copy(),
            meta=injprofs_meta.copy(),
            header=CDP_header.copy())

        injprof_cdp.init_HL_and_fillAll()

        injprof_cdp.hdulist[0].header.insert(list(CDP_header.keys())[0],
            ('title', 'CHINJ01: INJPROFS-%s' % direction.upper()))

        return injprof_cdp


    def dump_aggregated_results(self):
        """ """

        if self.report is not None:
            self.report.add_Section(keyword='dump', 
                Title='Aggregated Results', level=0)

            self.add_DataAlbaran2Report()

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        # Histogram of Slopes [ADU/electrons]

        # Histogram of Notch [ADU/electrons]

        # Histogram of IG1_THRESH

        # Injection level vs. Calibrated IG1, MAP

        CURVES_IG1CAL_MAP = self._get_CHIG1_MAP_from_PT(kind='CAL')

        figkey1 = 'CHINJ01_curves_MAP_IG1_CAL'
        figname1 = self.figs[figkey1]

        self.plot_XYMAP(CURVES_IG1CAL_MAP, **dict(
                        suptitle='Charge Injection Curves - Calibrated IG1',
                        doLegend=True,
                        ylabel='Inj [kel]',
                        xlabel='IG1 [V]',
                        corekwargs=dict(E=dict(linestyle='-', marker='', color='r'),
                                        F=dict(linestyle='-', marker='', color='g'),
                                        G=dict(linestyle='-', marker='', color='b'),
                                        H=dict(linestyle='-', marker='', color='m')),
                        figname=figname1
                        ))


        if self.report is not None:
            self.addFigure2Report(figname1, 
                        figkey=figkey1, 
                        caption='CHINJ01: Charge injection level [ke-] as a function of '+\
                            'calibrated IG1 voltage.', 
                        texfraction=0.7)

        # saving charge injection parameters to a json CDP

        ICURVES_PAR_MAP = self.get_FPAMAP_from_PT(
                    self.ParsedTable['CHINJ01'],
                    extractor=self._extract_INJCURVES_PAR_fromPT)

        ic_header = OrderedDict()
        ic_header['title'] = 'Injection Curves Parameters'
        ic_header['test'] = 'CHINJ01'
        ic_header.update(CDP_header)

        ic_meta = OrderedDict()
        ic_meta['units'] ='/2^16 ADU',
        ic_meta['model'] = 'I=b+1/(1+exp(-K(IG1-XT))) * (-A*(IG1-XN)[IG1<XN] + N)'
        ic_meta['structure'] = ''

        ic_cdp = cdp.Json_CDP(rootname=self.outcdps['INJCURVES'],
                         path=self.cdpspath)

        ic_cdp.ingest_inputs(data=ICURVES_PAR_MAP,
                header = ic_header,
                meta=ic_meta)
        ic_cdp.savehardcopy()


        # Injection level vs. Calibrated IG1, single plot

        IG1CAL_Singledict = self._get_XYdict_INJ(kind='CAL')

        figkey2 = 'CHINJ01_curves_IG1_CAL'
        figname2 = self.figs[figkey2]

        IG1CAL_kwargs = dict(
            title='Charge Injection Curves - Calibrated IG1',
            doLegend=False,
            xlabel='IG1 (Calibrated) [V]',
            ylabel='Injection [kel]',
            figname=figname2)

        corekwargs = dict()
        for block in self.flight_blocks:
            for iCCD in self.CCDs:
                corekwargs['%s_CCD%i_E' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#FF4600')  # red
                corekwargs['%s_CCD%i_F' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#61FF00')  # green
                corekwargs['%s_CCD%i_G' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#00FFE0')  # cyan
                corekwargs['%s_CCD%i_H' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#1700FF')  # blue

        IG1CAL_kwargs['corekwargs'] = corekwargs.copy()

        self.plot_XY(IG1CAL_Singledict, **IG1CAL_kwargs)

        if self.report is not None:
            self.addFigure2Report(figname2, 
                        figkey=figkey2, 
                        caption='CHINJ01: Charge injection level [ke-] as a function of '+\
                            'calibrated IG1 voltage.', 
                        texfraction=0.7)        

        # Injection level vs. Non-Calibrated IG1, single plot


        IG1RAW_Singledict = self._get_XYdict_INJ(kind='RAW')

        figkey3 = 'CHINJ01_curves_IG1_RAW'
        figname3 = self.figs[figkey3]

        IG1RAW_kwargs = dict(
            title='Charge Injection Curves - RAW IG1',
            doLegend=False,
            xlabel='IG1 (RAW) [V]',
            ylabel='Injection [kel]',
            figname=figname3)

        corekwargs = dict()
        for block in self.flight_blocks:
            for iCCD in self.CCDs:
                corekwargs['%s_CCD%i_E' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#FF4600')  # red
                corekwargs['%s_CCD%i_F' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#61FF00')  # green
                corekwargs['%s_CCD%i_G' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#00FFE0')  # cyan
                corekwargs['%s_CCD%i_H' % (block, iCCD)] = dict(linestyle='-',
                                                                marker='', color='#1700FF')  # blue

        IG1RAW_kwargs['corekwargs'] = corekwargs.copy()

        self.plot_XY(IG1RAW_Singledict, **IG1RAW_kwargs)

        if self.report is not None:
            self.addFigure2Report(figname3, 
                        figkey=figkey3, 
                        caption='CHINJ01: Charge injection level [ke-] as a function of '+\
                            'Non-calibrated IG1 voltage.', 
                        texfraction=0.7)        

        # Notch level vs. calibrated IG2

        # Notch level vs. calibrated IDL

        # Notch level vs. calibrated OD

        # Notch injection map, ADUs

        NOTCHADUMAP = self.get_FPAMAP_from_PT(
            self.ParsedTable['CHINJ01'],
            extractor=self._get_extractor_NOTCH_fromPT(
                units='ADU'))

        figkey4 = 'NOTCH_ADU_MAP'
        figname4 = self.figs[figkey4]

        self.plot_SimpleMAP(NOTCHADUMAP, **dict(
            suptitle='CHINJ01: NOTCH INJECTION [ADU]',
            ColorbarText='ADU',
            figname=figname4))

        if self.report is not None:
            self.addFigure2Report(figname4, 
                        figkey=figkey4, 
                        caption='CHINJ01: notch injection level, in ADU.', 
                        texfraction=0.7)

        # Notch injection map, ELECTRONs

        NOTCHEMAP = self.get_FPAMAP_from_PT(self.ParsedTable['CHINJ01'],
                                            extractor=self._get_extractor_NOTCH_fromPT(units='E'))

        figkey5 = 'NOTCH_ELE_MAP'
        figname5 = self.figs[figkey5]

        self.plot_SimpleMAP(NOTCHEMAP, **dict(
            suptitle='CHINJ01: NOTCH INJECTION [ELECTRONS]',
            ColorbarText='electrons',
            figname=figname5))

        if self.report is not None:
            self.addFigure2Report(figname5, 
                        figkey=figkey5, 
                        caption='CHINJ01: notch injection level, in electrons.', 
                        texfraction=0.7)

        # Average injection profiles

        IG1profs = 4.5

        xlabels_profs = dict(hor='column [pix]',
                            ver='row [pix]')

        ylabels_profs = dict(hor='Injection level [Normalized]',
                            ver='Injection level [ADU]',)

        proftypes = ['hor','ver']
        ccdhalves = ['top','bot']

        BLOCKcolors = cm.rainbow(np.linspace(0, 1, len(self.flight_blocks)))


        pointcorekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            for iCCD in self.CCDs:
                for kQ in self.Quads:
                    pointcorekwargs['%s_CCD%i_%s' % (block, iCCD, kQ)] = dict(
                        linestyle='', marker='.', color=jcolor, ms=2.0)


        for ccdhalf in ccdhalves:

            if ccdhalf == 'top':
                _Quads = ['G','H']
            elif ccdhalf == 'bot':
                _Quads = ['E','F']

            for proftype in proftypes:

                if proftype == 'hor':
                    xrangeNorm = None
                elif proftype == 'ver':
                    xrangeNorm = [10,20]


                XY_profs = self._get_XYdict_PROFS(proftype=proftype,
                    IG1=IG1profs,Quads=_Quads, doNorm=True, 
                    xrangeNorm=xrangeNorm)

                figkey6 = 'PROFS_%s_%s' % (proftype.upper(),ccdhalf.upper())
                figname6 = self.figs[figkey6]

                title = 'CHINJ01: Direction: %s, CCDHalf: %s' % \
                    (proftype.upper(),ccdhalf.upper()),

                if proftype == 'ver':
                    xlim=[0,50]
                    ylim=None
                elif proftype == 'hor':
                    xlim=None
                    ylim=[0.5,1.5]

                profkwargs = dict(
                    title=title,
                    doLegend=False,
                    xlabel=xlabels_profs[proftype],
                    xlim=xlim,
                    ylim=ylim,
                    ylabel=ylabels_profs[proftype],
                    figname=figname6,
                    corekwargs=pointcorekwargs)

                self.plot_XY(XY_profs, **profkwargs)

                if proftype == 'ver':
                    captemp = 'CHINJ01: Average (normalized) injection profiles in vertical direction (along CCD columns) '+\
                        'for IG1=%.2fV. Only the 2 channels in the CCD %s-half are shown '+\
                        '(%s, %s). Each colour corresponds to a '+\
                        'different block (2x3 quadrant-channels in each colour).'
                elif proftype == 'hor':
                    captemp = 'CHINJ01: Average injection profiles in horizontal direction (along CCD rows) '+\
                        'for IG1=%.2fV. The profiles have been normalized by the median injection level. '+\
                        'Only the 2 channels in the CCD %s-half are shown (%s, %s). Each colour corresponds to a '+\
                        'different block (2x3 quadrant-channels in each colour).'

                if self.report is not None:
                    self.addFigure2Report(figname6, 
                        figkey=figkey6, 
                        caption= captemp % (IG1profs, ccdhalf, _Quads[0],_Quads[1]),
                        texfraction=0.7)


        # Average injection vertical profiles, zoomed in to highlight
        # non-perfect charge injection shut-down.

        
        pointcorekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            for iCCD in self.CCDs:
                for kQ in self.Quads:
                    pointcorekwargs['%s_CCD%i_%s' % (block, iCCD, kQ)] = dict(
                        linestyle='', marker='.', color=jcolor, ms=2.0)


        for ccdhalf in ccdhalves:

            if ccdhalf == 'top':
                _Quads = ['G','H']
            elif ccdhalf == 'bot':
                _Quads = ['E','F']


            XY_profs = self._get_XYdict_PROFS(proftype='ver',
                IG1=IG1profs,Quads=_Quads, doNorm=True,
                xrangeNorm=[10,20])

            figkey7 = 'PROFS_ver_%s_ZOOM' % (ccdhalf.upper(),)
            figname7 = self.figs[figkey7]

            title = 'CHINJ01: Direction: ver, CCDHalf: %s, ZOOM-in' % \
                (ccdhalf.upper(),),
            
            xlim=[25,50]
            ylim=[0,4.e-3]
            
            profkwargs = dict(
                title=title,
                doLegend=False,
                xlabel=xlabels_profs[proftype],
                xlim=xlim,
                ylim=ylim,
                ylabel=ylabels_profs[proftype],
                figname=figname7,
                corekwargs=pointcorekwargs)

            self.plot_XY(XY_profs, **profkwargs)

            captemp = 'CHINJ01: Average injection profiles in vertical direction (along CCD columns) '+\
                    'for IG1=%.2fV. Only the 2 channels in the CCD %s-half are shown '+\
                    '(%s, %s). Each colour corresponds to a '+\
                    'different block (2x3 quadrant-channels in each colour). Zoomed in '+\
                    'to highlight injection shutdown profile.'
            
            if self.report is not None:
                self.addFigure2Report(figname7, 
                    figkey=figkey7, 
                    caption= captemp % (IG1profs, ccdhalf, _Quads[0],_Quads[1]),
                    texfraction=0.7)

        # creating and saving INJ PROFILES CDPs.

        for direction in ['hor','ver']:

            _injprof_xlsx_cdp = self.get_injprof_xlsx_cdp(direction=direction,
                inCDP_header=CDP_header)
            _injprof_xlsx_cdp.savehardcopy()

            _injprof_fits_cdp = self.get_injprof_fits_cdp(direction=direction,
                inCDP_header=CDP_header)
            _injprof_fits_cdp.savehardcopy()


        # reporting non-uniformity of injection lines to report

        if self.report is not None:

            NUN_HOR = self.get_FPAMAP_from_PT(self.ParsedTable['CHINJ01'],
                        extractor=self._extract_NUNHOR_fromPT)

            nun_cdpdict = dict(
                caption='CHINJ01: Non-Uniformity of the injection lines, rms, as percentage.',
                valformat='%.2f')

            ignore = self.add_StdQuadsTable2Report( 
                            Matrix = NUN_HOR,
                            cdpdict = nun_cdpdict)

