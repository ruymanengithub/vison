#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Created on Thu Aug 22 10:33:00 2019

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

class MetaNL(MetaCal):
    """ """
    
    def __init__(self, *args, **kwargs):
        """ """
        
        super(MetaNL,self).__init__(*args,**kwargs)
        
        self.testnames = ['NL02']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy() 
        
        self.products['NL'] = OrderedDict()
    
    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """
                
        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']
        
        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]
        
        IndexS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0])])
        
        IndexC = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
                                     vcore.vIndex('CCD',vals=self.CCDs)])
        
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
        
        
        
        NL_dict = sidd.products['NL'].copy()
 
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        maxNLpc_v = tmp_v_CQ.copy()
        fluADU_maxNLpc_v = tmp_v_CQ.copy()
        fluELE_maxNLpc_v = tmp_v_CQ.copy()        
        stabilitypc_v = tmp_v_CQ.copy()
        
        nlcdpkeys_v = np.zeros((1,NCCDs),dtype='S50')
        
        for iCCD, CCDk in enumerate(CCDkeys):
            
            nlcdpkey = '%s_%s_%s_%i_%s' % (testname, block, session,jrep+1,CCDk)
            
            for kQ, Q in enumerate(self.Quads):
                
                maxNLpc_v[0,iCCD,kQ] = NL_dict[CCDk][Q]['maxNLpc']
                
                gain = self.cdps['GAIN'][block][CCDk][Q][0]
                
                fluADU_maxNLpc_v[0,iCCD,kQ] = NL_dict[CCDk][Q]['flu_maxNLpc']
                fluELE_maxNLpc_v[0,iCCD,kQ] = gain * fluADU_maxNLpc_v[0,iCCD,kQ]
                
                stabilitypc_v[0,iCCD,kQ] = NL_dict[CCDk][Q]['stability_pc']
                
            self.products['NL'][nlcdpkey] = NL_dict.copy()
            
            nlcdpkeys_v[0,iCCD] = nlcdpkey
            
                
        sidd.addColumn(maxNLpc_v, 'MAXNL', IndexCQ)
        sidd.addColumn(fluADU_maxNLpc_v, 'FLUADU_MAXNL', IndexCQ)
        sidd.addColumn(fluELE_maxNLpc_v, 'FLUELE_MAXNL', IndexCQ)
        sidd.addColumn(stabilitypc_v, 'STABILITYPC', IndexCQ)
        
        
        sidd.addColumn(nlcdpkeys_v, 'NLCDP_KEY', IndexC)

        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit

    def dump_aggregated_results(self):
        """ """
        
        outpathroot = self.outpath
        
        stop()
        