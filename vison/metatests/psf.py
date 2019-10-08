#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 14:19:00 2019

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
         'offset_pre', 'offset_ove', 'bgd_img', 'std_pre', 'std_ove', 'chk_x', 'chk_x_ccd', 'chk_y', 'chk_y_ccd', 'chk_peak', 'chk_fluence', 'chk_fwhmx', 'chk_fwhmy']

class MetaPsf(MetaCal):
    """ """
    
    def __init__(self, *args, **kwargs):
        """ """
        
        super(MetaPsf,self).__init__(*args,**kwargs)
        
        self.Spots = ['ALPHA','BRAVO','CHARLIE','DELTA','ECHO']
        self.testnames = ['PSF01_590','PSF01_730','PSF01_800','PSF01_880']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        #self.blocks = self.blocks[1:] # TESTS!
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy() 
        
        self.products['XTALK'] = OrderedDict()
        
        self.init_fignames()
    
    
    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """
                
        
        #NCCDs = len(self.CCDs)
        #NQuads = len(self.Quads)
        session = inventoryitem['session']
        
        #CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]
        
        IndexS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0])])
        
        #IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
        #                             vcore.vIndex('CCD',vals=self.CCDs),
        #                             vcore.vIndex('Quad',vals=self.Quads)])
        
        #IndexCQS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
        #                             vcore.vIndex('CCD',vals=self.CCDs),
        #                             vcore.vIndex('Quad',vals=self.Quads),
        #                             vcore.vIndex('Spot', vals=self.Spots)])
    
        
        #idd = copy.deepcopy(inventoryitem['dd'])
        sidd = self.parse_single_test_gen(jrep, block, testname, inventoryitem)
        
        # TEST SCPECIFIC
        # TO BE ADDED:            
        #   BLOCK, TEST, REPEAT
        #   wavenm, calibrated HK (voltages),
        

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
 
        
        xtalkpath = os.path.join(inventoryitem['resroot'],'xtalk')
        
        
        ctalkcdp_pick = os.path.join(xtalkpath,os.path.split(sidd.products['CTALK'])[-1])        
        ctalkcdp = files.cPickleRead(ctalkcdp_pick)
                
        ctalkkey = 'XTALK_%s_%s_%s_%i' % (testname, block, session,jrep+1)
        self.products['XTALK'][ctalkkey] = ctalkcdp['data'].copy()
        ctalkkey_v = np.array([ctalkkey])
        sidd.addColumn(ctalkkey_v, 'XTALK', IndexS, ix=0)
        
        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit
    
    def get_XTALKDICT_from_PT(self,testname):
        """ """
        
        XTALKdict = OrderedDict()
        XTALKdict['blocks'] = self.flight_blocks
        
        PT = self.ParsedTable[testname]
                 
        for block in self.flight_blocks:
            ixblock = np.where(PT['BLOCK'] == block)[0][0]
            ctalkkey = PT['XTALK'][ixblock]
            XTALKdict[block] = self.products['XTALK'][ctalkkey].copy()
        
        return XTALKdict
    
    def plot_XtalkMAP(self,XTALKs, kwargs):
        """ """
        
        xtalkmap = plfpa.XtalkPlot(XTALKs,**kwargs)
        if 'figname' in kwargs:
            figname = kwargs['figname']
        else:
            figname = ''
        xtalkmap.render(figname=figname)
        
    
    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        for testname in self.testnames:
            self.figs['XTALK_MAP_%s' % testname] = os.path.join(self.figspath,
                         'XTALK_MAP_%s.png' % testname)
   
    def dump_aggregated_results(self):
        """ """
        
        for testname in self.testnames:
            
            XTALKs = self.get_XTALKDICT_from_PT(testname)            
        
            stestname = st.replace(testname,'_','\_')
            self.plot_XtalkMAP(XTALKs,kwargs=dict(
                    scale='ADU',
                    showvalues=False,
                    title='%s: XTALK [ADU]' % stestname,
                    figname=self.figs['XTALK_MAP_%s' % testname]))
        