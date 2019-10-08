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

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files
from vison.datamodel import ccd as ccdmod


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
         'HK_VID_PCB_TEMP_T', 'HK_VID_PCB_TEMP_B', 'HK_RPSU_TEMP1', 'HK_FPGA_PCB_TEMP_T', 'HK_FPGA_PCB_TEMP_B', 'HK_RPSU_TEMP_2', 'HK_RPSU_28V_PRI_I']

class MetaCosmetics(MetaCal):
    """ """
    
    def __init__(self, *args, **kwargs):
        """ """
        
        super(MetaCosmetics,self).__init__(*args,**kwargs)
        
        self.testnames = ['COSMETICS00']        
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        #self.blocks = self.blocks[-2:] # TESTS
        
        self.products['MASKSCOOS'] = OrderedDict()
        self.maskkeys = ['DARK','FLAT','MERGE']
        
        self.init_fignames()
    
    def _extract_badpix_coordinates(self,all_mask_fits, block):
        """ """
        
        coordinates = OrderedDict()
        
        for CCD in self.CCDs:
            CCDkey = 'CCD%i' % CCD
            
            ccdobj = ccdmod.CCD(all_mask_fits[CCDkey])
            
            Ckey = self.fpa.get_Ckey_from_BlockCCD(block,CCD)
            
            if Ckey is None:
                flip = (1,0) # block outside FPA (reserves)
            else:
                flip = flip = self.fpa.FPA_MAP[Ckey][-1]
            
            
            coordinates[CCDkey] = OrderedDict()
            
            img = ccdobj.extensions[-1].data.copy()
            
            flipimg = self.fpa.flip_img(img,flip)
            
            x,y = np.where(flipimg > 0)
            
            coordinates[CCDkey] = OrderedDict(x=x,
                                             y=y)
                        
        
        return coordinates
    
    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """
        
        
        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']
        
        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]
        
        IndexS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0])])
        
        #IndexC = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
        #                             vcore.vIndex('CCD',vals=self.CCDs)])
        
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
        
        
        block_v = np.array([block])            
        sidd.addColumn(block_v, 'BLOCK', IndexS, ix=0)
        
        
        test_v = np.array([jrep+1])            
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)
        
        test_v = np.array([session])            
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)
        
        test_v = np.array([testname])            
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)
        
        # ADDING Nr. pixels masked: DARK, FLAT, MERGEs
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        Ndark_v = tmp_v_CQ.copy()
        Nflat_v = tmp_v_CQ.copy()
        Nmerge_v = tmp_v_CQ.copy()
        Ncols_v = tmp_v_CQ.copy()
        
        productspath = os.path.join(inventoryitem['resroot'],'products')       
        
        NPIX_pick = os.path.join(productspath,os.path.split(sidd.products['DEF_TB_CDP'])[-1])
        
        NPIX_dict = files.cPickleRead(NPIX_pick)['data']['DEFECTS'].copy()
                
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                ixsel = np.where((NPIX_dict['CCD']==(iCCD+1)) & (NPIX_dict['Q']==(kQ+1)))
                
                Ndark_v[0,iCCD,kQ] = NPIX_dict['N_DARK'].as_matrix()[ixsel][0]
                Nflat_v[0,iCCD,kQ] = NPIX_dict['N_FLAT'].as_matrix()[ixsel][0] 
                Nmerge_v[0,iCCD,kQ] = NPIX_dict['N_MERGE'].as_matrix()[ixsel][0] 
                Ncols_v[0,iCCD,kQ] = NPIX_dict['NCOLS_MERGE'].as_matrix()[ixsel][0] 
                
        sidd.addColumn(Ndark_v, 'NDARK', IndexCQ)
        sidd.addColumn(Nflat_v, 'NFLAT', IndexCQ)
        sidd.addColumn(Nmerge_v, 'NMERGE', IndexCQ)
        sidd.addColumn(Ncols_v, 'NCOLS', IndexCQ)
        
        # Maps of bad pixels
        
        tmp_v_str = np.zeros((1),dtype='S50')
        
        for maskkey in self.maskkeys:
            
            mskkey_v = tmp_v_str.copy()
            
            all_mask_fits = dict()
            
            for jCCD, CCDk in enumerate(CCDkeys):
                
                all_mask_fits[CCDk] = os.path.join(productspath,os.path.split(sidd.products[maskkey][CCDk])[-1])
                
            _mskkey = '%s_%s_%s_R%i' % (maskkey, block, session, jrep+1) 
                
            self.products['MASKSCOOS'][_mskkey] = self._extract_badpix_coordinates(all_mask_fits, block)
                
            mskkey_v[0] = _mskkey
        
            sidd.addColumn(mskkey_v,'MASKCOOS_%s' % maskkey,IndexS)
            
        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        
        return sit
        
    def _get_extractor_LOGNDEF_fromPT(self,maskkey):
        
        def _extract_LOGNDEF_fromPT(PT, block, CCDk, Q):
            ixblock = np.where(PT['BLOCK'].data == block)
            column = 'N%s_%s_Quad%s' % (maskkey,CCDk,Q)            
            LOGNDEF = np.log10(max(PT[column][ixblock][0],1))
            return LOGNDEF
        
        return _extract_LOGNDEF_fromPT
    
    def _get_NBADCOLS_fromPT(self, PT, block, CCDk, Q):
        ixblock = np.where(PT['BLOCK'].data == block)
        column = 'NCOLS_%s_Quad%s' % (CCDk,Q)
        return PT[column][ixblock][0]
    
    def _get_DEFECTSMAP_from_PT(self,maskkey):
        """ """
        
        DEFMAP = OrderedDict()
        
        PT = self.ParsedTable['COSMETICS00']
        column = 'MASKCOOS_%s' % maskkey
        
        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                
                Ckey  = 'C_%i%i' % (jY+1,iX+1)
                DEFMAP[Ckey] = OrderedDict()
                
                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                
                ixblock = np.where(PT['BLOCK'] == block)
                
                _maskkey = PT[column][ixblock][0]
                
                _coodict = self.products['MASKSCOOS'][_maskkey]
                
                
                DEFMAP[Ckey] = _coodict[CCDk]
                
        
        return DEFMAP
    
    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        self.figs['DEFECTS_MAP'] = os.path.join(self.figspath,
                         'DEFECTS_MAP.png')
        
        self.figs['DEFECTS_COUNTS_MAP'] = os.path.join(self.figspath,
                         'DEFECTS_COUNTS_MAP.png')
        self.figs['BADCOLS_COUNTS_MAP'] = os.path.join(self.figspath,
                         'BADCOLS_COUNTS_MAP.png')
        
        
    def dump_aggregated_results(self):
        """ """
        
        # Maps of Defects
        
        for maskkey in self.maskkeys:
            
            DEFMAP = self._get_DEFECTSMAP_from_PT(maskkey)
            
            self.plot_XYMAP(DEFMAP,kwargs=dict(
                        suptitle='DEFECTS: %s' % maskkey,
                        figname=self.figs['DEFECTS_MAP']
                        ))
        

        # Defects Counts
        
        for maskkey in self.maskkeys:
            
            NDEFMAP = self.get_FPAMAP_from_PT(self.ParsedTable['COSMETICS00'],
                                              extractor=self._get_extractor_LOGNDEF_fromPT(maskkey))
            
            self.plot_SimpleMAP(NDEFMAP, kwargs=dict(
                        suptitle='NR. DEFECTS: %s' % maskkey,
                        figname=self.figs['DEFECTS_COUNTS_MAP']
                        ))

        # Bad Column Counts
        
        
        NBADCOLSMAP = self.get_FPAMAP_from_PT(self.ParsedTable['COSMETICS00'],
                                              extractor=self._get_NBADCOLS_fromPT)
            
        self.plot_SimpleMAP(NBADCOLSMAP, kwargs=dict(
            suptitle='NR. of BAD COLUMNS [MERGED MASK]',
            figname=self.figs['BADCOLS_COUNTS_MAP']
            ))
        
