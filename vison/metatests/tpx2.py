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

from vison.support import files
from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse


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
    'std_ove',
    'chk_mea_inject',
    'chk_med_inject',
    'chk_std_inject']


class MetaTPX2(MetaCal):
    """ """

    def __init__(self, **kwargs):
        """ """

        super(MetaTPX2, self).__init__(**kwargs)

        self.testnames = ['TPX2']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['CHINJ'] = OrderedDict()
        self.products['CHINJNOISE'] = OrderedDict()
        self.products['MERGEDL1'] = OrderedDict()
        self.products['MERGEDL2'] = OrderedDict()

        self.init_fignames()

    def load_block_results(self, inventoryfile=None):
        """ """

        self.testnames = ['TP02', 'TP21']

        super(MetaTPX2, self).load_block_results(inventoryfile)

        # Renaming tests for convenience in the inventory

        for block in self.blocks:

            oldtestkey = self.inventory[block].keys()[0]
            self.inventory[block]['TPX2'] = copy.deepcopy(self.inventory[block][oldtestkey])
            self.inventory[block].pop(oldtestkey)

        self.testnames = ['TPX2']

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

        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        inject_top_v = tmp_v_CQ.copy()
        inject_top_std_v = tmp_v_CQ.copy()

        inject_bot_v = tmp_v_CQ.copy()
        inject_bot_std_v = tmp_v_CQ.copy()

        id_dly = inventoryitem['dd'].mx['id_dly'][:].copy()
        s_tpump = inventoryitem['dd'].mx['s_tpump'][:].copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                if Q in ['G', 'H']:

                    ixsel = np.where((id_dly[:, iCCD] == id_dly.min()) &
                                     (s_tpump[:, iCCD] == 1))

                    inject_top_v[0, iCCD, kQ] = np.nanmedian(
                        inventoryitem['dd'].mx['chk_med_inject'][ixsel, iCCD, kQ])
                    inject_top_std_v[0, iCCD, kQ] = np.nanmedian(
                        inventoryitem['dd'].mx['chk_std_inject'][ixsel, iCCD, kQ])

                elif Q in ['E', 'F']:

                    ixsel = np.where((id_dly[:, iCCD] == id_dly.max()) &
                                     (s_tpump[:, iCCD] == 1))

                    inject_bot_v[0, iCCD, kQ] = np.nanmedian(
                        inventoryitem['dd'].mx['chk_med_inject'][ixsel, iCCD, kQ])
                    inject_bot_std_v[0, iCCD, kQ] = np.nanmedian(
                        inventoryitem['dd'].mx['chk_std_inject'][ixsel, iCCD, kQ])

        sidd.addColumn(inject_top_v, 'INJ_TOP', IndexCQ)
        sidd.addColumn(inject_top_std_v, 'INJ_STD_TOP', IndexCQ)

        sidd.addColumn(inject_bot_v, 'INJ_BOT', IndexCQ)
        sidd.addColumn(inject_bot_std_v, 'INJ_STD_BOT', IndexCQ)

        #productspath = os.path.join(inventoryitem['resroot'],'products')

        chinjkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['CHINJ'][chinjkey] = sidd.products['chinj_mx'].copy()
        sidd.addColumn(np.array([chinjkey]), 'CHINJ', IndexS)
        #print('added CHINJ')

        chinjnoisekey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['CHINJNOISE'][chinjnoisekey] = sidd.products['chinjnoise_mx'].copy()
        sidd.addColumn(np.array([chinjnoisekey]), 'CHINJNOISE', IndexS)
        #print('added CHINJNOISE')

        mergedL1key = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['MERGEDL1'][mergedL1key] = copy.deepcopy(sidd.products['mergedL1_df'])
        sidd.addColumn(np.array([mergedL1key]), 'MERGEDL1', IndexS)
        #print('added MERGEDL1')

        mergedL2key = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['MERGEDL2'][mergedL2key] = copy.deepcopy(sidd.products['mergedL2_df'])
        sidd.addColumn(np.array([mergedL2key]), 'MERGEDL2', IndexS)
        #print('added MERGEDL2')

        dip_count_v = tmp_v_CQ.copy()
        dip_tau_v = tmp_v_CQ.copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                dip_count_v[0, iCCD, kQ] = sidd.products['mergedL2_df'].loc[CCDk, Q]['N']
                dip_tau_v[0, iCCD, kQ] = sidd.products['mergedL2_df'].loc[CCDk, Q]['<tau>']

        sidd.addColumn(dip_count_v, 'DIP_COUNT', IndexCQ)
        sidd.addColumn(dip_tau_v, 'DIP_TAU', IndexCQ)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _extract_NDIP_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(self, PT, block)
        column = 'DIP_COUNT_%s_Quad%s' % (CCDk, Q)
        NDIP = np.nanmedian(PT[column][ixblock])

        return NDIP

    def _extract_TAU_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(self, PT, block)
        column = 'DIP_TAU_%s_Quad%s' % (CCDk, Q)
        AVTAU = np.nanmedian(PT[column][ixblock])
        if np.isnan(AVTAU):
            AVTAU = 0.
        return AVTAU

    def _extract_INJ_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(self, PT, block)
        if Q in ['E', 'F']:
            column = 'INJ_BOT_%s_Quad%s' % (CCDk, Q)
        elif Q in ['G', 'H']:
            column = 'INJ_TOP_%s_Quad%s' % (CCDk, Q)

        INJ = np.nanmedian(PT[column][ixblock])
        return INJ

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['NDIP_MAP'] = os.path.join(self.figspath,
                                             'NDIP_MAP_TPX2.png')

        self.figs['TAU_MAP'] = os.path.join(self.figspath,
                                            'TAU_MAP_TPX2.png')

        self.figs['INJ_MAP'] = os.path.join(self.figspath,
                                            'INJ_MAP_TPX2.png')

    def dump_aggregated_results(self):
        """ """

        # MAP with total numbers of dipoles measured

        NDIPMAP = self.get_FPAMAP_from_PT(self.ParsedTable['TPX2'],
                                          extractor=self._extract_NDIP_fromPT)

        self.plot_SimpleMAP(NDIPMAP, **dict(
            suptitle='TPX2: Nr. OF DIPOLES',
            figname=self.figs['NDIP_MAP']))

        # MAP with tau's

        TAUMAP = self.get_FPAMAP_from_PT(self.ParsedTable['TPX2'],
                                         extractor=self._extract_TAU_fromPT)

        self.plot_SimpleMAP(TAUMAP, **dict(
            suptitle=r'$TPX2:\ <TAU>\ [us]$',
            figname=self.figs['TAU_MAP']))

        # INJECTION level

        INJMAP = self.get_FPAMAP_from_PT(self.ParsedTable['TPX2'],
                                         extractor=self._extract_INJ_fromPT)

        self.plot_SimpleMAP(INJMAP, **dict(
            suptitle=r'$TPX2: CHARGE INJECTION [ADU]]$',
            figname=self.figs['INJ_MAP']))
