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

from vison.support import vjson
from vison.support import files
from vison.fpa import fpa as fpamod

from vison.fpatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse



from matplotlib import pyplot as plt
plt.switch_backend('TkAgg')
from matplotlib.colors import Normalize
# END IMPORT

cols2keep = ['test', 'sn_ccd1', 'sn_ccd2', 'sn_ccd3', 'sn_roe', 'sn_rpsu', 'exptime', 'vstart', 'vend', 'rdmode', 'flushes', 'siflsh', 'siflsh_p', 'swellw', 'swelldly', 'inisweep', 
         'cdpu_clk', 'chinj', 'chinj_on', 'chinj_of', 'id_wid', 'id_dly', 'chin_dly', 'v_tpump', 's_tpump', 'v_tp_mod', 's_tp_mod', 'v_tp_cnt', 's_tp_cnt', 'dwell_v', 'dwell_s', 'toi_fl', 'toi_tp', 
         'toi_ro', 'toi_ch', 'motr', 'motr_cnt', 'motr_siz', 'source', 'wave', 'mirr_on', 'mirr_pos', 'R1C1_TT', 'R1C1_TB', 'R1C2_TT', 'R1C2_TB', 'R1C3_TT', 'R1C3_TB', 'IDL', 'IDH', 'IG1_1_T', 'IG1_2_T', 
         'IG1_3_T', 'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T', 'IG2_B', 'OD_1_T', 'OD_2_T', 'OD_3_T', 'OD_1_B', 'OD_2_B', 'OD_3_B', 'RD_T', 'RD_B', 'time', 'HK_CCD1_TEMP_T', 'HK_CCD2_TEMP_T', 
         'HK_CCD3_TEMP_T', 'HK_CCD1_TEMP_B', 'HK_CCD2_TEMP_B', 'HK_CCD3_TEMP_B', 'HK_CCD1_OD_T', 'HK_CCD2_OD_T', 'HK_CCD3_OD_T', 'HK_CCD1_OD_B', 'HK_CCD2_OD_B', 'HK_CCD3_OD_B', 'HK_COMM_RD_T', 
         'HK_COMM_RD_B', 'HK_CCD1_IG1_T', 'HK_CCD2_IG1_T', 'HK_CCD3_IG1_T', 'HK_CCD1_IG1_B', 'HK_CCD2_IG1_B', 'HK_CCD3_IG1_B', 'HK_COMM_IG2_T', 'HK_COMM_IG2_B', 'HK_FPGA_BIAS_ID2', 
         'HK_VID_PCB_TEMP_T', 'HK_VID_PCB_TEMP_B', 'HK_RPSU_TEMP1', 'HK_FPGA_PCB_TEMP_T', 'HK_FPGA_PCB_TEMP_B', 'HK_RPSU_TEMP_2', 'HK_RPSU_28V_PRI_I', 'chk_NPIXOFF', 'chk_NPIXSAT',
         'offset_pre','offset_img','offset_ove','std_pre','std_img','std_ove']


class MetaMOT(MetaCal):
    """ """
    
    def __init__(self, **kwargs):
        """ """
        
        super(MetaMOT,self).__init__(**kwargs)
        
        self.testnames = ['MOT_WARM']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
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
        
        IndexS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0])])
        
        
        IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
                                     vcore.vIndex('CCD',vals=self.CCDs),
                                     vcore.vIndex('Quad',vals=self.Quads)])
        
    
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
                
        test_v = np.array([jrep+1])            
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)
        
        test_v = np.array([session])            
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)
        
        test_v = np.array([testname])            
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)
 
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        
        rwdvs_off_v = tmp_v_CQ.copy()
        rwdvs_ron_v = tmp_v_CQ.copy()
        
                
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                rwdvs_off_v[0,iCCD,kQ] = sidd.products['rwdvs_off_mx'][iCCD,kQ]
                rwdvs_ron_v[0,iCCD,kQ] = sidd.products['rwdvs_ron_mx'][iCCD,kQ]
                
                
        sidd.addColumn(rwdvs_off_v, 'RWDVS_OFF', IndexCQ)        
        sidd.addColumn(rwdvs_ron_v, 'RWDVS_RON', IndexCQ)
        
        # Extracting Injection level from profiles
        profilespath = os.path.join(inventoryitem['resroot'],'profiles')
        prof_pick = os.path.join(profilespath,os.path.split(sidd.products['MW_PROFILES'])[-1])
        
        profiles = files.cPickleRead(prof_pick)
        
        chinj_on = sidd.meta['inputs']['structure']['col004']['chinj_on']
        chinj_of = sidd.meta['inputs']['structure']['col004']['chinj_of']
        
        
        injection_v = tmp_v_CQ.copy()
        
                
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                inj_arr = profiles['data']['CHINJ'][CCDk][Q]['y'].copy()
        
                inj_on = np.median(inj_arr[0:chinj_on])
                inj_of = np.median(inj_arr[chinj_of:chinj_on+chinj_of])
                
                inj_net = inj_on - inj_of
                
                injection_v[0,iCCD,kQ] = inj_net
        
        
        sidd.addColumn(injection_v, 'INJECTION', IndexCQ)   
        
        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit

    def _extract_INJ_fromPT(self, PT, block, CCDk, Q):
        ixblock = np.where(PT['BLOCK'].data == block)
        column = 'INJECTION_%s_Quad%s' % (CCDk,Q)
            
        injection = PT[column][ixblock][0]
        return injection
    
    def _extract_RON_fromPT(self, PT, block, CCDk, Q):
        ixblock = np.where(PT['BLOCK'].data == block)
        column = 'RWDVS_RON_%s_Quad%s' % (CCDk,Q)
            
        ron = PT[column][ixblock][0]
        return ron
    
    def _extract_OFF_RWDVS_fromPT(self, PT, block, CCDk, Q):
        ixblock = np.where(PT['BLOCK'].data == block)
        column = 'RWDVS_OFF_%s_Quad%s' % (CCDk,Q)
            
        off = int(PT[column][ixblock][0])
        return off
    
    def _extract_OFF_PREFWD_fromPT(self, PT, block, CCDk, Q):
        ixblock = np.where(PT['BLOCK'].data == block)
        column = 'offset_pre_%s_Quad%s' % (CCDk,Q)
            
        off = int(PT[column][ixblock][0])
        return off
    
    
    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        self.figs['INJ_MAP'] = os.path.join(self.figspath,
                         'INJECTION_MAP.png')
        
        self.figs['INJ_MAP_json'] = os.path.join(self.figspath,
                         'INJECTION_MAP.json')

        self.figs['RON_MAP'] = os.path.join(self.figspath,
                         'RON_MAP.png')
        
        self.figs['RON_MAP_json'] = os.path.join(self.figspath,
                         'RON_MAP.json')
        
        self.figs['OFF_RWDVS_MAP'] = os.path.join(self.figspath,
                         'OFF_RWDVS_MAP.png')
        
        self.figs['OFF_RWDVS_MAP_json'] = os.path.join(self.figspath,
                         'OFF_RWDVS_MAP.json')
        
        self.figs['OFF_PREFWD_MAP'] = os.path.join(self.figspath,
                         'OFF_PREFWD_MAP.png')
        
        self.figs['OFF_PREFWD_MAP_json'] = os.path.join(self.figspath,
                         'OFF_PREFWD_MAP.json')

    def dump_aggregated_results(self):
        """ """
        
        # Injection map., ADUs
        
        
        INJMAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                    extractor=self._extract_INJ_fromPT)
        
        vjson.save_jsonfile(INJMAP,self.figs['INJ_MAP_json'])
        
        
        self.plot_SimpleMAP(INJMAP,kwargs=dict(
                        suptitle='MOT\_WARM: INJECTION LEVEL [ADU]',
                        figname = self.figs['INJ_MAP']
                        ))
        
        # RON map, ADUs 
        
        RONMAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                    extractor=self._extract_RON_fromPT)
        
        vjson.save_jsonfile(RONMAP,self.figs['RON_MAP_json'])
        
        self.plot_SimpleMAP(RONMAP,kwargs=dict(
                        suptitle='MOT\_WARM: RON [ADU]',
                        figname = self.figs['RON_MAP']
                        ))
        
        # RON map, ELECTRONs
        
        # OFFSET map (RWDVS)
        
        OFF_RWDVS_MAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                    extractor=self._extract_OFF_RWDVS_fromPT)
        
        vjson.save_jsonfile(OFF_RWDVS_MAP,self.figs['OFF_RWDVS_MAP_json'])
        
        self.plot_SimpleMAP(OFF_RWDVS_MAP,kwargs=dict(
                        suptitle='MOT\_WARM: OFFSET, RWDVS [ADU]',
                        figname = self.figs['OFF_RWDVS_MAP']
                        ))
        
        # OFFSET map (FWD, PRE)
        
        OFF_PREFWD_MAP = self.get_FPAMAP_from_PT(self.ParsedTable['MOT_WARM'],
                    extractor=self._extract_OFF_PREFWD_fromPT)
        
        vjson.save_jsonfile(OFF_PREFWD_MAP,self.figs['OFF_PREFWD_MAP_json'])
        
        self.plot_SimpleMAP(OFF_PREFWD_MAP,kwargs=dict(
                        suptitle='MOT\_WARM: OFFSET, PRESCAN FWD [ADU]',
                        figname = self.figs['OFF_PREFWD_MAP']
                        ))

        
        
        
        