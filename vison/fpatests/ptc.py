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

from vison.fpa import fpa as fpamod

from vison.fpatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files


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
         'offset_pre', 'offset_ove', 'deltaoff_pre', 'deltaoff_ove', 'flu_med_img', 'flu_std_img', 'std_pre', 'std_ove']

class MetaPTC(MetaCal):
    """ """
    
    def __init__(self, *args, **kwargs):
        """ """
        
        super(MetaPTC,self).__init__(*args,**kwargs)
        
        self.testnames = ['PTC01','PTC02_590','PTC02_730','PTC02_880']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
    
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
        #   BLOCK, TEST, REPEAT
        #   wavenm, calibrated HK (voltages),
        #   GAIN, EGAIN, ALPHA, BLOOM_ADU, BLOOM_E
        #   REFERENCES TO CURVES

        CHAMBER = sidd.meta['inputs']['CHAMBER']
        
        
        CHAMBER_key = CHAMBER[0]        
        chamber_v = np.array([CHAMBER_key])
        sidd.addColumn(chamber_v, 'CHAMBERKEY', IndexS, ix=0)
        
        ogseobj = ogse.Ogse(CHAMBER=CHAMBER)
        
        wave = sidd.mx['wave'][0,0]
        
        wave_v = np.array([ogseobj.get_wavelength(wave)])
        sidd.addColumn(wave_v, 'WAVENM', IndexS, ix=0)
        
        block_v = np.array([block])            
        sidd.addColumn(block_v, 'BLOCK', IndexS, ix=0)
        
        
        test_v = np.array([jrep+1])            
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)
        
        test_v = np.array([session])            
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)
        
        test_v = np.array([testname])            
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)
 
        gain_mx = sidd.products['gain_mx']
        bloom_mx = sidd.products['bloom_mx']
        
        HER_dict = sidd.compliances['HER']
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        gain_v = tmp_v_CQ.copy()
        egain_v = tmp_v_CQ.copy()
        alpha_v = tmp_v_CQ.copy()
        bloom_ADU_v = tmp_v_CQ.copy()
        bloom_e_v = tmp_v_CQ.copy()
        HER_v = tmp_v_CQ.copy()
        
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                gain_v[0,iCCD,kQ] = gain_mx[CCDk][Q]['gain']
                egain_v[0,iCCD,kQ] = gain_mx[CCDk][Q]['egain']
                alpha_v[0,iCCD,kQ] = gain_mx[CCDk][Q]['alpha']
                
                bloom_ADU_v[0,iCCD,kQ] = bloom_mx[CCDk][Q]['bloom_ADU']
                bloom_e_v[0,iCCD,kQ] = bloom_mx[CCDk][Q]['bloom_e']
                
                HER_v[0,iCCD,kQ] = HER_dict[CCDk][Q][1]
                
        sidd.addColumn(gain_v, 'GAIN', IndexCQ)
        sidd.addColumn(egain_v, 'EGAIN', IndexCQ)
        sidd.addColumn(alpha_v, 'ALPHA', IndexCQ)
        sidd.addColumn(bloom_ADU_v, 'BLOOM_ADU', IndexCQ)
        sidd.addColumn(bloom_e_v, 'BLOOM_E', IndexCQ)
        sidd.addColumn(HER_v, 'HER', IndexCQ)

        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit
        
   
    def _extract_GAIN_fromPT(self,PT,block,CCDk,Q):
        """ """
        ixblock = np.where(PT['BLOCK'].data == block)        
        column = 'GAIN_%s_Quad%s' % (CCDk,Q)        
        G = PT[column][ixblock][0]        
        return G
    
    def _extract_HER_fromPT(self,PT,block,CCDk,Q):
        """ """
        ixblock = np.where(PT['BLOCK'].data == block)        
        column = 'HER_%s_Quad%s' % (CCDk,Q)        
        HER = PT[column][ixblock][0]        
        return HER
        
    def _extract_BADU_fromPT(self,PT,block,CCDk,Q):
        """ """
        ixblock = np.where(PT['BLOCK'].data == block)
        
        column = 'BLOOM_ADU_%s_Quad%s' % (CCDk,Q)
        badu = PT[column][ixblock][0]
        if badu>0:
            return badu
        else:
            return np.nan


    def _extract_BE_fromPT(self,PT,block,CCDk,Q):
        """ """
        ixblock = np.where(PT['BLOCK'].data == block)
        
        column = 'BLOOM_E_%s_Quad%s' % (CCDk,Q)
        be = PT[column][ixblock][0]
        if be>0:
            return be
        else:
            return np.nan
    
    def gen_GAIN_MXdict(self):
        """ """
        
        G_MX = OrderedDict()
        
        for block in self.blocks:
            
            G_MX[block] = OrderedDict()
            for testname in self.testnames:
                    
                PT = self.ParsedTable[testname]
                
                ixblock = np.where(PT['BLOCK'].data == block)
                
                G_MX[block][testname] = OrderedDict()
                
                for iCCD in self.CCDs:
                    CCDk = 'CCD%i' % iCCD
                    G_MX[block][testname][CCDk] = OrderedDict()
                    for Q in self.Quads:
                        G = PT['GAIN_%s_Quad%s' % (CCDk,Q)].data[ixblock][0]
                        EG = PT['EGAIN_%s_Quad%s' % (CCDk,Q)].data[ixblock][0]
                        gpair = (G,EG)
                        G_MX[block][testname][CCDk][Q] = gpair
                            
        return G_MX

    def dump_aggregated_results(self):
        """ """
        
        outpathroot = self.outpath
        
        # GAIN matrix (all blocks and test/waves) to dict() saved as pickle
        
        GAIN_MXdict = self.gen_GAIN_MXdict()
        
        GAIN_MXpick = os.path.join(outpathroot,'GAIN_MX_PTC0X.pick')
        
        files.cPickleDump(GAIN_MXdict,GAIN_MXpick)
        
        
        # GAIN maps (all tests/waves)
        
        for testname in self.testnames:
            
            GMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_GAIN_fromPT)
        
            stestname = st.replace(testname,'_','\_')
                        
            self.plot_SimpleMAP(GMAP,kwargs=dict(
                    suptitle='%s: GAIN e-/ADU' % stestname))
        
        
        
        
        # BLOOM maps (ADU and e-, from PTC01)
        
        for testname in self.testnames:
            
            BADU_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_BADU_fromPT)
        
            stestname = st.replace(testname,'_','\_')
            self.plot_SimpleMAP(BADU_MAP,kwargs=dict(
                    suptitle='%s: BLOOM-ADU [DN]' % stestname))
        
        for testname in self.testnames:
            
            BE_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_BE_fromPT)
        
            stestname = st.replace(testname,'_','\_')
            self.plot_SimpleMAP(BE_MAP,kwargs=dict(
                    suptitle='%s: BLOOM-ELECTRONS' % stestname))
        
        # HER map
        
        
        
        # GAIN vs. detector temperature (PTC01)
        
        # GAIN vs. OD-CAL (PTC01)
        
        # GAIN vs. RD-CAL (PTC01)
        
        
        # Save the ParsedTable(s)
        
        
        
        