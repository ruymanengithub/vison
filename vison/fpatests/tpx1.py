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
         'offset_pre', 'offset_ove', 'std_pre', 'std_ove', 'chk_mea_inject', 'chk_med_inject', 'chk_std_inject']
         
class MetaTPX1(MetaCal):
    """ """
    
    def __init__(self, **kwargs):
        """ """
        
        super(MetaTPX1,self).__init__(**kwargs)
        
        self.testnames = ['TP01','TP11']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy() 
        
        #self.products['METAFIT'] = OrderedDict()
        
    
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
        
        if testname in ['TP01','TP11']:
            nickname = 'TPX1'
        
        test_v = np.array([nickname])            
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)
        
        productspath = os.path.join(inventoryitem['resroot'],'products')
        
        print('RE-WRITTEN UP TO HERE!')
        stop()
        
#        metafitcdp_pick = os.path.join(productspath,os.path.split(sidd.products['METAFIT_CDP'])[-1])
#        metafitcdp = files.cPickleRead(metafitcdp_pick)
#        metafit = copy.deepcopy(metafitcdp['data']['ANALYSIS'])
#        
#        metafitkey = '%s_%s_%s_%i' % (testname, block, session,jrep+1)
#        self.products['METAFIT'][metafitkey] = copy.deepcopy(metafit)
#        metafitkey_v = np.array([metafitkey])
#        sidd.addColumn(metafitkey_v, 'METAFIT', IndexS, ix=0)
#        
#        metacdp_pick = os.path.join(productspath,os.path.split(sidd.products['META_CDP'])[-1]) # change to META_CDP
#        metacdp = files.cPickleRead(metacdp_pick)
#        meta = metacdp['data']['ANALYSIS'] # this is a pandas DataFrame
#        
#        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
#                
#        bgd_adu_v = tmp_v_CQ.copy()
#        ig1_thresh_v = tmp_v_CQ.copy()
#        ig1_notch_v = tmp_v_CQ.copy()
#        slope_v = tmp_v_CQ.copy()
#        n_adu_v = tmp_v_CQ.copy()        
#        
#        
#        for iCCD, CCDk in enumerate(CCDkeys):
#            for kQ, Q in enumerate(self.Quads):
#                
#                ixloc = np.where((meta['CCD'] == iCCD+1) & (meta['Q'] == kQ+1))
#                
#                bgd_adu_v[0,iCCD,kQ] = meta['BGD_ADU'][ixloc[0][0]]
#                ig1_thresh_v[0,iCCD,kQ] = meta['IG1_THRESH'][ixloc[0][0]]
#                ig1_notch_v[0,iCCD,kQ] = meta['IG1_NOTCH'][ixloc[0][0]]
#                slope_v[0,iCCD,kQ] = meta['S'][ixloc[0][0]]
#                n_adu_v[0,iCCD,kQ] = meta['N_ADU'][ixloc[0][0]]
#                       
#        
#        sidd.addColumn(bgd_adu_v, 'FIT_BGD_ADU', IndexCQ)
#        sidd.addColumn(ig1_thresh_v, 'FIT_IG1_THRESH', IndexCQ)
#        sidd.addColumn(ig1_notch_v, 'FIT_IG1_NOTCH', IndexCQ)        
#        sidd.addColumn(slope_v, 'FIT_SLOPE', IndexCQ)
#        sidd.addColumn(n_adu_v, 'FIT_N_ADU', IndexCQ)
#        # flatten sidd to table
#        
#        sit = sidd.flattentoTable()
#        
#        return sit



    def dump_aggregated_results(self):
        """ """
        
        
        outpathroot = self.outpath
        
        
        
        

        
        
        
        