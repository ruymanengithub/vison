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
         'IG1_3_T', 'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T', 'IG2_B', 'OD_1_T', 'OD_2_T', 'OD_3_T', 'OD_1_B', 'OD_2_B', 'OD_3_B', 'RD_T', 'RD_B', 'time',  'chk_NPIXOFF', 'chk_NPIXSAT', 
         'offset_pre', 'offset_ove', 'deltaoff_pre', 'deltaoff_ove', 'flu_med_img', 'flu_std_img', 'std_pre', 'std_ove']

class MetaBF(MetaCal):
    """ """
    
    def __init__(self, *args, **kwargs):
        """ """
        
        super(MetaBF,self).__init__(*args,**kwargs)
        
        self.testnames = ['BF01','BF01_590','BF01_730','BF01_880']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy() 
    
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
                
        productspath = os.path.join(inventoryitem['resroot'],'products')
        BFfitCDP_pick = os.path.join(productspath,os.path.split(sidd.products['BFfitTABLE_CDP'])[-1])
        
        BFfitCDP = files.cPickleRead(BFfitCDP_pick)
        
        BFfit_df = BFfitCDP['data']['BFFIT']
        
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        fwhmx_hwc_v = tmp_v_CQ.copy()
        fwhmy_hwc_v = tmp_v_CQ.copy()
        
        fwhmx_slope_v = tmp_v_CQ.copy()
        fwhmy_slope_v = tmp_v_CQ.copy()
        
        ell_v = tmp_v_CQ.copy()
        
        
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                ixsel = np.where((BFfit_df['CCD'] == iCCD+1) & (BFfit_df['Q'] == kQ+1))
                
                fwhmx_hwc_v[0,iCCD,kQ] = BFfit_df['FWHMx_HWC'].as_matrix()[ixsel][0]
                fwhmy_hwc_v[0,iCCD,kQ] = BFfit_df['FWHMy_HWC'].as_matrix()[ixsel][0]
                
                fwhmx_slope_v[0,iCCD,kQ] = BFfit_df['FWHMx_Slope'].as_matrix()[ixsel][0]
                fwhmy_slope_v[0,iCCD,kQ] = BFfit_df['FWHMy_Slope'].as_matrix()[ixsel][0]
                
                ell_v[0,iCCD,kQ] = BFfit_df['ELL_HWC'].as_matrix()[ixsel][0]
        
        sidd.addColumn(fwhmx_hwc_v, 'FWHMX_HWC', IndexCQ)
        sidd.addColumn(fwhmy_hwc_v, 'FWHMY_HWC', IndexCQ)
        sidd.addColumn(fwhmx_slope_v, 'FWHMX_SLOPE', IndexCQ)
        sidd.addColumn(fwhmy_slope_v, 'FWHMY_SLOPE', IndexCQ)
        sidd.addColumn(ell_v, 'ELL_HWC', IndexCQ)
        
        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit
        
   

    def dump_aggregated_results(self):
        """ """
        
        outpathroot = self.outpath
        
        stop()
        