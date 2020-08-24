#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 16:28:00 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
from collections import OrderedDict
import string as st
import os

#from vison.support import vjson
from vison.datamodel import cdp as cdpmod
from vison.support import files
from vison.fpa import fpa as fpamod
from vison.support import utils

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
    'offset_img',
    'offset_ove',
    'std_pre',
    'std_img',
    'std_ove']


class MetaMOT(MetaCal):
    """ """

    def __init__(self, **kwargs):
        """ """

        super(MetaMOT, self).__init__(**kwargs)

        self.testnames = ['MOT_WARM']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['RAMPS'] = OrderedDict()

        self.init_fignames()
        self.init_outcdpnames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """

        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']

        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        IndexC = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0]),
                                    vcore.vIndex('CCD', vals=self.CCDs)])

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


        for reg in ['pre','img','ove']:

            _fwd_off_v = tmp_v_CQ.copy()
            _rwdv_off_v = tmp_v_CQ.copy()
            _rwdvs_off_v = tmp_v_CQ.copy()        
            _rwdvs_ron_v = tmp_v_CQ.copy()

            for iCCD, CCDk in enumerate(CCDkeys):
                for kQ, Q in enumerate(self.Quads):

                    _fwd_off_v[0, iCCD, kQ] = sidd.products['fwd_off_%s_mx' % reg][iCCD, kQ]
                    _rwdv_off_v[0, iCCD, kQ] = sidd.products['rwdv_off_%s_mx' % reg][iCCD, kQ]
                    _rwdvs_off_v[0, iCCD, kQ] = sidd.products['rwdvs_off_%s_mx' % reg][iCCD, kQ]
                    _rwdvs_ron_v[0, iCCD, kQ] = sidd.products['rwdvs_ron_%s_mx' % reg][iCCD, kQ]

            sidd.addColumn(_fwd_off_v, 'FWD_%s_OFF' % reg.upper(), IndexCQ)
            sidd.addColumn(_rwdv_off_v, 'RWDV_%s_OFF' % reg.upper(), IndexCQ)
            sidd.addColumn(_rwdvs_off_v, 'RWDVS_%s_OFF' % reg.upper(), IndexCQ)
            sidd.addColumn(_rwdvs_ron_v, 'RWDVS_%s_RON' % reg.upper(), IndexCQ)

        # Extracting Injection level from profiles
        profilespath = os.path.join(inventoryitem['resroot'], 'profiles')
        prof_pick = os.path.join(profilespath, os.path.split(sidd.products['MW_PROFILES'])[-1])

        profiles = files.cPickleRead(prof_pick)

        chinj_on = sidd.meta['inputs']['structure']['col004']['chinj_on']
        chinj_of = sidd.meta['inputs']['structure']['col004']['chinj_of']

        injection_v = tmp_v_CQ.copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                inj_arr = profiles['data']['CHINJ'][CCDk][Q]['y'].copy()

                inj_on = np.median(inj_arr[0:chinj_on])
                inj_of = np.median(inj_arr[chinj_of:chinj_on + chinj_of])

                inj_net = inj_on - inj_of

                injection_v[0, iCCD, kQ] = inj_net

        sidd.addColumn(injection_v, 'INJECTION', IndexCQ)

        rampkeys_v = np.zeros((1, NCCDs), dtype='U50')

        for iCCD, CCDk in enumerate(CCDkeys):

            rampccdkey = '%s_%s_%s_%i_%s' % (testname, block, session, jrep + 1, CCDk)

            RAMP_dict = dict()

            for kQ, Q in enumerate(self.Quads):

                RAMP_dict[Q] = OrderedDict()

                rampprof = profiles['data']['RAMP'][CCDk][Q].copy()

                offset = sidd.mx['offset_pre'][0, iCCD, kQ]

                rampslope = np.median(rampprof['y'][1:100] - rampprof['y'][0:99])
                rampoff = rampprof['y'][0] - offset

                RAMP_dict[Q]['slope'] = rampslope
                RAMP_dict[Q]['offset'] = rampoff

            self.products['RAMPS'][rampccdkey] = RAMP_dict.copy()

            rampkeys_v[0, iCCD] = rampccdkey

        sidd.addColumn(rampkeys_v, 'RAMPKEYS', IndexC)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _extract_INJ_fromPT(self, PT, block, CCDk, Q):
        ixblock = self.get_ixblock(PT, block)
        column = 'INJECTION_%s_Quad%s' % (CCDk, Q)

        injection = PT[column][ixblock][0]
        return injection

    def _extract_RON_fromPT(self, PT, block, CCDk, Q):
        ixblock = self.get_ixblock(PT, block)
        column = 'RWDVS_OVE_RON_%s_Quad%s' % (CCDk, Q)

        ron = PT[column][ixblock][0]
        return ron

    def _get_extractor_OFF_fromPT(self, readout, reg):

        def _extractor(PT, block, CCDk, Q):

            ixblock = self.get_ixblock(PT, block)

            if reg != 'all':
                column = '%s_%s_OFF_%s_Quad%s' % (readout,reg.upper(), CCDk, Q)                            
                off = int(PT[column][ixblock][0])
            elif reg == 'all':
                off = OrderedDict()
                for _reg in ['pre','img','ove']:
                    column = '%s_%s_OFF_%s_Quad%s' % (readout,_reg.upper(), CCDk, Q)
                    off[_reg] = int(PT[column][ixblock][0])

            return off

        return _extractor

    def _extract_RAMP_fromPT(self, PT, block, CCDk, Q):
        ixblock = self.get_ixblock(PT, block)
        column = 'RAMPKEYS_%s' % (CCDk,)
        rampkey = PT[column][ixblock][0]

        rampdict = self.products['RAMPS'][rampkey][Q].copy()
        for key in rampdict:
            rampdict[key] = float(rampdict[key])

        return rampdict

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)


        self.figs['INJ_MAP'] = os.path.join(self.figspath,
                                            'MOT_WARM_INJECTION_MAP.png')


        self.figs['RON_MAP'] = os.path.join(self.figspath,
                                            'MOT_WARM_RON_MAP.png')


        self.figs['OFF_RWDVS_MAP'] = os.path.join(self.figspath,
                                                  'MOT_WARM_OFF_RWDVS_MAP.png')


        self.figs['OFF_PREFWD_MAP'] = os.path.join(self.figspath,
                                                   'MOT_WARM_OFF_PREFWD_MAP.png')


    def init_outcdpnames(self):
        """ """

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        self.outcdps['RAMP_MAP_json'] = 'MOT_WARM_RAMP_MAP.json'


        self.outcdps['INJ_MAP_json'] = 'MOT_WARM_INJECTION_MAP.json'


        self.outcdps['RON_MAP_json'] ='MOT_WARM_RON_MAP.json'

        self.outcdps['OFF_RWDVS_MAP_json'] = 'MOT_WARM_OFF_RWDVS_MAP.json'
        self.outcdps['OFF_RWDV_MAP_json'] = 'MOT_WARM_OFF_RWDV_MAP.json'
        self.outcdps['OFF_FWD_MAP_json'] ='MOT_WARM_OFF_FWD_MAP.json'



    def dump_aggregated_results(self):
        """ """

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))


        # SAVING FWD RAMPs' slope & intercept

        RAMPMAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                          extractor=self._extract_RAMP_fromPT)

        rm_header = OrderedDict()
        rm_header['title'] = 'RAMP MAP'
        rm_header.update(CDP_header)

        rm_cdp = cdpmod.Json_CDP(rootname=self.outcdps['RAMP_MAP_json'],
                              path=self.cdpspath)
        rm_cdp.ingest_inputs(data=RAMPMAP,
                             header = rm_header,
                             meta=dict(units='ADU/row, ADU'))
        rm_cdp.savehardcopy()


        # Injection map., ADUs

        INJMAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                         extractor=self._extract_INJ_fromPT)

        ci_header = OrderedDict()
        ci_header['title'] = 'CHARGE INJECTION MAP'
        ci_header.update(CDP_header)

        ci_cdp = cdpmod.Json_CDP(rootname=self.outcdps['INJ_MAP_json'],
                              path=self.cdpspath)
        ci_cdp.ingest_inputs(data=INJMAP,
                             header = ci_header,
                             meta=dict(units='ADU'))
        ci_cdp.savehardcopy()


        self.plot_SimpleMAP(INJMAP, **dict(
            suptitle='MOT\_WARM: INJECTION LEVEL [ADU]',
            figname=self.figs['INJ_MAP']
        ))

        # RON map, ADUs

        RONMAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                         extractor=self._extract_RON_fromPT)


        rn_header = OrderedDict()
        rn_header['title'] = 'RON MAP [RWDVS]'
        rn_header.update(CDP_header)

        rn_cdp = cdpmod.Json_CDP(rootname=self.outcdps['RON_MAP_json'],
                              path=self.cdpspath)
        rn_cdp.ingest_inputs(data=RONMAP,
                             header = rn_header,
                             meta=dict(units='ADU'))
        rn_cdp.savehardcopy()

        self.plot_SimpleMAP(RONMAP, **dict(
            suptitle='MOT\_WARM: RON [ADU]',
            figname=self.figs['RON_MAP']
        ))

        # RON map, ELECTRONs


        # OFFSET map (RWDV)

        OFF_RWDV_MAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                                extractor=self._get_extractor_OFF_fromPT('RWDV','all'))

        offrwdv_header = OrderedDict()
        offrwdv_header['title'] = 'OFFSET MAP [RWDV]'
        offrwdv_header.update(CDP_header)

        offrwdv_cdp = cdpmod.Json_CDP(rootname=self.outcdps['OFF_RWDV_MAP_json'],
                              path=self.cdpspath)
        offrwdv_cdp.ingest_inputs(data=OFF_RWDV_MAP,
                             header = offrwdv_header,
                             meta=dict(units='ADU'))
        offrwdv_cdp.savehardcopy()



        # OFFSET map (RWDVS)

        OFF_RWDVS_MAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                                extractor=self._get_extractor_OFF_fromPT('RWDVS','all'))

        offrwdvs_header = OrderedDict()
        offrwdvs_header['title'] = 'OFFSET MAP [RWDVS]'
        offrwdvs_header.update(CDP_header)

        offrwdvs_cdp = cdpmod.Json_CDP(rootname=self.outcdps['OFF_RWDVS_MAP_json'],
                              path=self.cdpspath)
        offrwdvs_cdp.ingest_inputs(data=OFF_RWDVS_MAP,
                             header = offrwdvs_header,
                             meta=dict(units='ADU'))
        offrwdvs_cdp.savehardcopy()

        OFF_RWDVS_IMG_MAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                                extractor=self._get_extractor_OFF_fromPT('RWDVS','img'))

        self.plot_SimpleMAP(OFF_RWDVS_IMG_MAP, **dict(
            suptitle='MOT\_WARM: OFFSET-IMG, RWDVS [ADU]',
            figname=self.figs['OFF_RWDVS_MAP']
        ))

        # OFFSET map (FWD, 'all')

        OFF_FWD_MAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                                 extractor=self._get_extractor_OFF_fromPT('FWD','all'))

        offfwd_header = OrderedDict()
        offfwd_header['title'] = 'OFFSET MAP [FWD]'
        offfwd_header.update(CDP_header)

        offfwd_cdp = cdpmod.Json_CDP(rootname=self.outcdps['OFF_FWD_MAP_json'],
                              path=self.cdpspath)
        offfwd_cdp.ingest_inputs(data=OFF_FWD_MAP,
                             header = offfwd_header,
                             meta=dict(units='ADU'))
        offfwd_cdp.savehardcopy()

         # OFFSET map (FWD, PRE)

        OFF_PREFWD_MAP  = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                                                 extractor=self._get_extractor_OFF_fromPT('FWD','pre'))

        self.plot_SimpleMAP(OFF_PREFWD_MAP, **dict(
            suptitle='MOT\_WARM: OFFSET, PRE-SCAN FWD [ADU]',
            figname=self.figs['OFF_PREFWD_MAP']
        ))
