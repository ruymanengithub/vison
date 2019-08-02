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
         'offset_pre', 'offset_ove', 'deltaoff_pre', 'deltaoff_ove', 'flu_med_img', 'flu_std_img', 'std_pre', 'std_ove']

class MetaPTC(MetaCal):
    """ """
    
    def __init__(self, *args, **kwargs):
        
        super(MetaPTC,self).__init__(*args,**kwargs)
        
        self.testnames = ['PTC01','PTC02_590','PTC02_730','PTC02_880']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
    
    def _parse_single_test(self, jrep, block, testname, inventoryitem):
        
        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        
        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]
        
        IndexS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0])])
        
        IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
                                     vcore.vIndex('CCD',vals=self.CCDs),
                                     vcore.vIndex('Quad',vals=self.Quads)])
        
        idd = copy.deepcopy(inventoryitem['dd'])
        session = inventoryitem['session']
        
        
        CHAMBER = idd.meta['inputs']['CHAMBER']
        
        ogseobj = ogse.Ogse(CHAMBER=CHAMBER)
        
        
        sidd = self.stack_dd(idd,self.incols,
                             indices2keep=['ix','CCD','Quad'],
                             index2stack = 'ix',
                             stacker='median')
        
        sidd.dropColumn('test')
        
        # rename the CCD index values in sidd, for convenience
        
        sidd.indices[sidd.indices.names.index('CCD')].vals = self.CCDs
                    
        for col in sidd.colnames:
            if 'CCD' in sidd.mx[col].indices.get_names():
                _i = sidd.mx[col].indices.names.index('CCD')
                sidd.mx[col].indices[_i].vals = self.CCDs
        
        #resroot = self.inventory[block][testname]['resroot']
        
        
        # TO BE ADDED:            
        #   BLOCK, TEST, REPEAT
        #   wavenm, calibrated HK (voltages),
        #   GAIN, EGAIN, ALPHA, BLOOM_ADU, BLOOM_E
        #   REFERENCES TO CURVES
        
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
        
 
        gain_mx = idd.products['gain_mx']
        bloom_mx = idd.products['bloom_mx']
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        gain_v = tmp_v_CQ.copy()
        egain_v = tmp_v_CQ.copy()
        alpha_v = tmp_v_CQ.copy()
        bloom_ADU_v = tmp_v_CQ.copy()
        bloom_e_v = tmp_v_CQ.copy()
        
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                gain_v[0,iCCD,kQ] = gain_mx[CCDk][Q]['gain']
                egain_v[0,iCCD,kQ] = gain_mx[CCDk][Q]['egain']
                alpha_v[0,iCCD,kQ] = gain_mx[CCDk][Q]['alpha']
                
                bloom_ADU_v[0,iCCD,kQ] = bloom_mx[CCDk][Q]['bloom_ADU']
                bloom_e_v[0,iCCD,kQ] = bloom_mx[CCDk][Q]['bloom_e']
                
        sidd.addColumn(gain_v, 'GAIN', IndexCQ)
        sidd.addColumn(egain_v, 'EGAIN', IndexCQ)
        sidd.addColumn(alpha_v, 'ALPHA', IndexCQ)
        sidd.addColumn(bloom_ADU_v, 'BLOOM_ADU', IndexCQ)
        sidd.addColumn(bloom_e_v, 'BLOOM_E', IndexCQ)
        
        # CALIBRATED HK
        
        try:
            roeVCal = self.roeVCals[block]
        except KeyError:
            print('Voltage calibrations for block %s not found!' % block)
            roeVCal = None
        
        
        
        if roeVCal is not None:
                        
            for cal_key in ['OD','RD','IG1','IG2']:
                
                for CCD in self.CCDs:
                    for Q in self.Quads:
                                                
                        rHKkey= roeVCal.get_HKkey(cal_key, CCD, Q)
                        HKkey = 'HK_%s' % rHKkey.upper()
                        cHKkey = '%s_CAL' % HKkey
                        
                        HKV = sidd.mx[HKkey][0]
                        
                        HKVcal = roeVCal.fcal_HK(HKV, cal_key, CCD, Q)
                        
                        sidd.addColumn(np.zeros(1,dtype=float)+HKVcal, 
                                       cHKkey, IndexS)
            
        else:
            
            
            dummy_roeVCal = vcal.RoeVCal()
        
            cHKkeys = []
        
            for cal_key in ['OD','RD','IG1','IG2']:
            
                for CCD in self.CCDs:
                    for Q in self.Quads:
                                            
                        rHKkey= dummy_roeVCal.get_HKkey(cal_key, CCD, Q)
                        HKkey = 'HK_%s' % rHKkey.upper()
                        cHKkey = '%s_CAL' % HKkey
                    
                        cHKkeys.append(cHKkey)
            
            for cHKkey in cHKkeys:
                sidd.addColumn(np.zeros(1,dtype=float)+np.nan, 
                                       cHKkey, IndexS)
                
        

        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit
        
    
    
    def parse_test_results(self,testname):
        """ """
        
    
        for iblock, block in enumerate(self.blocks):
            
            try:
                Nreps = len(self.inventory[block][testname])
            except KeyError:
                print('block %s not found!' % block)
                continue


            for jrep in range(Nreps):
                
                inventoryitem = self.inventory[block][testname][jrep]
                
                sit = self._parse_single_test(jrep, block, testname, inventoryitem)
            
                
                # MERGING WITH PREVIOUS DDs
                
                if (iblock == 0) and (jrep == 0):
                    pt = copy.deepcopy(sit)
                else:
                    pt = self.stackTables(pt,sit)
        
        self.ParsedTable[testname] = pt
    
    def _extract_GAIN_fromPT(self,PT,block,CCDk,Q):        
        ixblock = np.where(PT['BLOCK'].data == block)        
        column = 'GAIN_%s_Quad%s' % (CCDk,Q)        
        G = PT[column][ixblock][0]        
        return G
    
    def plot_GMAP(self,GMAP,suptitle,filename=''):
        """ """
        
        Gs = []
        for ckey in GMAP.keys():
            for Q in self.Quads:
                Gs.append(GMAP[ckey][Q])
        
        
        normfunction = Normalize(vmin=np.min(Gs),vmax=np.max(Gs),clip=False)
   
        kwargs = dict(doColorbar=True,
              suptitle=suptitle,
              corekwargs=dict(
                      norm = normfunction                      
                      ))
    
        heatmap = plfpa.FpaHeatMap(GMAP, **kwargs)
        heatmap.render()
        
    def _extract_BADU_fromPT(self,PT,block,CCDk,Q):        
        ixblock = np.where(PT['BLOCK'].data == block)
        
        column = 'BLOOM_ADU_%s_Quad%s' % (CCDk,Q)
        badu = PT[column][ixblock][0]
        if badu>0:
            return badu
        else:
            return np.nan
    
    def _extract_BE_fromPT(self,PT,block,CCDk,Q):        
        ixblock = np.where(PT['BLOCK'].data == block)
        
        column = 'BLOOM_E_%s_Quad%s' % (CCDk,Q)
        be = PT[column][ixblock][0]
        if be>0:
            return be
        else:
            return np.nan
    
    def plot_BLOOM_MAP(self,BMAP,suptitle,filename=''):
        """ """
        
        Bs = []
        for ckey in BMAP.keys():
            for Q in self.Quads:
                Bs.append(BMAP[ckey][Q])
        
        
        normfunction = Normalize(vmin=np.nanmin(Bs),vmax=np.nanmax(Bs),clip=False)
   
        kwargs = dict(doColorbar=True,
              suptitle=suptitle,
              corekwargs=dict(
                      norm = normfunction                      
                      ))
    
        heatmap = plfpa.FpaHeatMap(BMAP, **kwargs)
        heatmap.render()    
    
    
    
    def dump_aggregated_results(self):
        """ """
        
        outpathroot = self.outpath
        
        # GAIN maps (all tests/waves)
        
        for testname in self.testnames:
            
            GMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_GAIN_fromPT)
        
            stestname = st.replace(testname,'_','\_')
            self.plot_GMAP(GMAP,suptitle='%s: GAIN e-/ADU' % stestname)
        
        
        # GAIN matrix (all blocks and test/waves)
        
        
        
        # BLOOM maps (ADU and e-, from PTC01)
        
        for testname in self.testnames:
            
            BADU_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_BADU_fromPT)
        
            stestname = st.replace(testname,'_','\_')
            self.plot_BLOOM_MAP(BADU_MAP,suptitle='%s: BLOOM-ADU [DN]' % stestname)
        
        for testname in self.testnames:
            
            BE_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_BE_fromPT)
        
            stestname = st.replace(testname,'_','\_')
            self.plot_BLOOM_MAP(BE_MAP,suptitle='%s: BLOOM-ELECTRONS' % stestname)
        
        
        # GAIN vs. detector temperature (PTC01)
        
        # GAIN vs. OD-CAL (PTC01)
        
        # GAIN vs. RD-CAL (PTC01)
        
        # Save the ParsedTable(s)
        
        
        
        