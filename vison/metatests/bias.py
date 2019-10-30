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
from vison.support import vjson


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
         'offset_pre', 'offset_img', 'offset_ove', 'std_pre', 'std_img', 'std_ove', 'RON']

class MetaBias(MetaCal):
    """ """
    
    def __init__(self, **kwargs):
        """ """
        
        super(MetaBias,self).__init__(**kwargs)
        
        self.testnames = ['BIAS01','BIAS02']
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
        
        off_pre_v = tmp_v_CQ.copy()
        off_img_v = tmp_v_CQ.copy()
        off_ove_v = tmp_v_CQ.copy()
        
        ron_pre_v = tmp_v_CQ.copy()
        ron_img_v = tmp_v_CQ.copy()
        ron_ove_v = tmp_v_CQ.copy()
        
        productspath = os.path.join(inventoryitem['resroot'],'products')
        
        roncdp_pick = os.path.join(productspath,os.path.split(sidd.products['RON_CDP'])[-1])
        roncdp = files.cPickleRead(roncdp_pick)
        
        offcdp_pick = os.path.join(productspath,os.path.split(sidd.products['OFF_CDP'])[-1])
        offcdp = files.cPickleRead(offcdp_pick)
        
        
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                off_pre_v[0,iCCD,kQ] = offcdp['data']['OFF_PRE'][CCDk][Q]
                off_img_v[0,iCCD,kQ] = offcdp['data']['OFF_IMG'][CCDk][Q]
                off_ove_v[0,iCCD,kQ] = offcdp['data']['OFF_OVE'][CCDk][Q]
                
                ron_pre_v[0,iCCD,kQ] = roncdp['data']['RON_PRE'][CCDk][Q]
                ron_img_v[0,iCCD,kQ] = roncdp['data']['RON_IMG'][CCDk][Q]
                ron_ove_v[0,iCCD,kQ] = roncdp['data']['RON_OVE'][CCDk][Q]
                
                
        sidd.addColumn(off_pre_v, 'OFF_PRE', IndexCQ)
        sidd.addColumn(off_img_v, 'OFF_IMG', IndexCQ)
        sidd.addColumn(off_ove_v, 'OFF_OVE', IndexCQ)
        
        sidd.addColumn(ron_pre_v, 'RON_PRE', IndexCQ)
        sidd.addColumn(ron_img_v, 'RON_IMG', IndexCQ)
        sidd.addColumn(ron_ove_v, 'RON_OVE', IndexCQ)

        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit

    def _get_extractor_RON_fromPT(self,units):
        """ """
        
        def _extract_RON_fromPT(PT, block, CCDk, Q):
            ixblock = np.where(PT['BLOCK'].data == block)
            column = 'RON_OVE_%s_Quad%s' % (CCDk,Q)
            
            if units =='ADU':
                unitsConvFactor=1
            elif units == 'E':
                unitsConvFactor = self.cdps['GAIN'][block][CCDk][Q][0]
            
            RON = np.nanmedian(PT[column][ixblock]) * unitsConvFactor
            return RON
        
        return _extract_RON_fromPT
    
    def _extract_OFFSET_fromPT(self,PT, block, CCDk, Q):
        """ """
        ixblock = np.where(PT['BLOCK'].data == block)
        column = 'OFF_OVE_%s_Quad%s' % (CCDk,Q)            
        RON = np.nanmedian(PT[column][ixblock])
        return RON

    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        for testname in self.testnames:
            self.figs['RON_ADU_%s' % testname] = os.path.join(self.figspath,
                         'RON_ADU_MAP_%s.png' % testname)
            self.figs['RONADU_MAP_%s_json' % testname] = os.path.join(self.figspath,
                         'RON_ADU_MAP_%s.json' % testname)            
            self.figs['RON_ELE_%s' % testname] = os.path.join(self.figspath,
                         'RON_ELE_MAP_%s.png' % testname)
            self.figs['OFFSETS_%s' % testname] = os.path.join(self.figspath,
                         'OFFSETS_MAP_%s.png' % testname)
            self.figs['OFFSETS_%s_json' % testname] = os.path.join(self.figspath,
                         'OFFSETS_MAP_%s.json' % testname)
            
            
    def dump_aggregated_results(self):
        """ """
        
        
        # RON maps (all tests/waves)
        
        
        # RON maps, ADUs
        for testname in self.testnames:
            
            RONADUMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], 
                                                extractor=self._get_extractor_RON_fromPT(units='ADU'))
            
            vjson.save_jsonfile(RONADUMAP,self.figs['RONADU_MAP_%s_json' % testname])
            
            stestname = st.replace(testname,'_','\_')
            self.plot_SimpleMAP(RONADUMAP,**dict(
                    suptitle='%s: RON [ADU]' % stestname,
                    figname=self.figs['RON_ADU_%s' % testname]))
        
        # RON maps, ELECTRONs
        for testname in self.testnames:
            
            RONEMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], 
                                              extractor=self._get_extractor_RON_fromPT(units='E'))
        
            stestname = st.replace(testname,'_','\_')
            self.plot_SimpleMAP(RONEMAP,**dict(
                    suptitle='%s: RON [ELECTRONS]' % stestname,
                    figname=self.figs['RON_ELE_%s' % testname]))
        
        # OFFSET maps
        
        for testname in self.testnames:
            
            OFFMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_OFFSET_fromPT)
            
            vjson.save_jsonfile(OFFMAP,self.figs['OFFSETS_%s_json' % testname])
            
            stestname = st.replace(testname,'_','\_')
            self.plot_SimpleMAP(OFFMAP,**dict(
                    suptitle='%s: OFFSET' % stestname,
                    figname=self.figs['OFFSETS_%s' % testname]))
        
        
