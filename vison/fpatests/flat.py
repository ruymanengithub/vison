#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 11:49:00 2019

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
             'toi_ro', 'toi_ch', 'motr', 'motr_cnt', 'motr_siz', 'source', 'wave', 'mirr_on', 'mirr_pos', 'R1C1_TT', 'R1C1_TB', 'R1C2_TT', 'R1C2_TB', 'R1C3_TT', 'R1C3_TB', 'IDL', 
             'IDH', 'IG1_1_T', 'IG1_2_T', 'IG1_3_T', 'IG1_1_B', 'IG1_2_B', 'IG1_3_B', 'IG2_T', 'IG2_B', 'OD_1_T', 'OD_2_T', 'OD_3_T', 'OD_1_B', 'OD_2_B', 'OD_3_B', 'RD_T', 'RD_B', 
             'HK_CCD1_TEMP_T', 'HK_CCD2_TEMP_T', 'HK_CCD3_TEMP_T', 'HK_CCD1_TEMP_B', 'HK_CCD2_TEMP_B', 'HK_CCD3_TEMP_B', 'HK_CCD1_OD_T', 'HK_CCD2_OD_T', 'HK_CCD3_OD_T', 'HK_CCD1_OD_B', 'HK_CCD2_OD_B', 
             'HK_CCD3_OD_B', 'HK_COMM_RD_T', 'HK_COMM_RD_B', 'HK_CCD1_IG1_T', 'HK_CCD2_IG1_T', 'HK_CCD3_IG1_T', 'HK_CCD1_IG1_B', 'HK_CCD2_IG1_B', 'HK_CCD3_IG1_B', 'HK_COMM_IG2_T', 'HK_COMM_IG2_B', 
             'HK_FPGA_BIAS_ID2', 'HK_VID_PCB_TEMP_T', 'HK_VID_PCB_TEMP_B', 'HK_RPSU_TEMP1', 'HK_FPGA_PCB_TEMP_T', 'HK_FPGA_PCB_TEMP_B', 'HK_RPSU_TEMP_2', 'HK_RPSU_28V_PRI_I', 'chk_NPIXOFF', 
             'chk_NPIXSAT', 'offset_pre', 'offset_ove', 'flu_med_img', 'flu_std_img', 'std_pre', 'std_ove']


class MetaFlat(MetaCal):
    """ """
    
    def __init__(self, **kwargs):
        """ """
        
        super(MetaFlat,self).__init__(**kwargs)
        
        self.testnames = ['FLAT01','FLAT02_590', 'FLAT02_730','FLAT02_880']
        
        self.colkeys = dict()
        for test in self.testnames:
            if test == 'FLAT01':
                self.colkeys[test] = ['col001','col002','col003']
            else:
                self.colkeys[test] = ['col001','col002']
        
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy() 
        
        self.products['MASTERFLATS'] = OrderedDict()
        
        self.Ncols = dict()
        
        self.init_fignames()
        
    
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
        
        Ncols = sidd.meta['structure']['Ncols']
        if testname not in self.Ncols:
            self.Ncols[testname] = Ncols
        
        # TEST SCPECIFIC
        
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
        
        # ADDING PRNU AND FLUENCES
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
                
        prnu_pc_vs = OrderedDict()
        prnu_fluadu_vs = OrderedDict()
        prnu_fluele_vs = OrderedDict()
        for icol in range(1,Ncols+1):
            colkey = 'col%03i' % icol
            
            prnu_pc_vs[colkey] = tmp_v_CQ.copy()        
            prnu_fluadu_vs[colkey] = tmp_v_CQ.copy()        
            prnu_fluele_vs[colkey] = tmp_v_CQ.copy()
        
        productspath = os.path.join(inventoryitem['resroot'],'products')
        
        prnucdp_pick = os.path.join(productspath,os.path.split(sidd.products['PRNU_TB_CDP'])[-1])
        prnucdp = files.cPickleRead(prnucdp_pick)
        
        for icol in range(1, Ncols+1):
            
            colkey = 'col%03i' % icol
            
            prnutb = prnucdp['data']['PRNU_%s' % colkey]
            
            for jCCD, CCDk in enumerate(CCDkeys):
                for kQ, Q in enumerate(self.Quads):
                    
                    ixsel = np.where((prnutb['CCD']==(jCCD+1)) & (prnutb['Q']==(kQ+1)))                    
                    prnu_pc_vs[colkey][0,jCCD,kQ] = prnutb['PRNU_PC'].as_matrix()[ixsel][0]
                    
                    G = self.cdps['GAIN'][block][CCDk][Q][0]
                    
                    prnu_fluadu_vs[colkey][0,jCCD,kQ] = prnutb['AVFLUENCE'].as_matrix()[ixsel][0]
                    
                    prnu_fluele_vs[colkey][0,jCCD,kQ] = G * prnu_fluadu_vs[colkey][0,jCCD,kQ]
                    
        for icol in range(1, Ncols+1):
            
            colkey = 'col%03i' % icol
                    
            sidd.addColumn(prnu_fluadu_vs[colkey], 'PRNU_FLUADU_%s' % colkey.upper(), IndexCQ)
            sidd.addColumn(prnu_fluele_vs[colkey], 'PRNU_FLUELE_%s' % colkey.upper(), IndexCQ)
            sidd.addColumn(prnu_pc_vs[colkey], 'PRNU_PC_%s' % colkey.upper(), IndexCQ)
        
        # ADDING REFERENCES TO MASTER FLAT-FIELDS
        
        tmp_v_C = np.zeros((1,NCCDs),dtype='S50')
        
        for icol in range(1, Ncols+1):
            colkey = 'col%03i' % icol
            
            mfkey_v = tmp_v_C.copy()
            
            for jCCD, CCDk in enumerate(CCDkeys):
                
                mfkey = '%s_%s_%s_%i_%s_%s' % (testname, block, session,jrep+1,colkey,CCDk)
                
                mfkey_v[0,jCCD] = mfkey
                
                self.products['MASTERFLATS'][mfkey] = os.path.split(sidd.products['MasterFFs'][colkey][CCDk])[-1]
            
            sidd.addColumn(mfkey_v, 'MASTERFLATS_%s' % colkey.upper(), IndexC)
        

        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit

    def _get_extractor_PRNU_fromPT(self,colkey):
        """ """
        
        def _extract_PRNU_fromPT(PT, block, CCDk, Q):
            ixblock = np.where(PT['BLOCK'].data == block)
            column = 'PRNU_PC_%s_%s_Quad%s' % (colkey.upper(),CCDk,Q)
            
            PRNU = PT[column][ixblock][0]
            return PRNU
        
        return _extract_PRNU_fromPT

    def _get_XYdict_PRNUFLU(self,PT,Ncols):
        """ """
                
        x = dict()
        y = dict()
        
        labelkeys = []
        
        for icol in range(1,Ncols+1):
            colkey = 'col%03i' % icol
            labelkeys.append(colkey)
            
            x[colkey] = []
            y[colkey] = []
            
            for iCCD,CCD in enumerate(self.CCDs):
                for iQ, Q in enumerate(self.Quads):
                    Xcolumn = 'PRNU_FLUELE_%s_CCD%i_Quad%s' % (colkey.upper(),CCD,Q)
                    Ycolumn = 'PRNU_PC_%s_CCD%i_Quad%s' % (colkey.upper(),CCD,Q)
                    
                    x[colkey] += (PT[Xcolumn]/1.E3).tolist()
                    y[colkey] += PT[Ycolumn].tolist()

            x[colkey] = np.array(x[colkey])
            y[colkey] = np.array(y[colkey])
        
        XYdict = dict(x=x,y=y,labelkeys=labelkeys)
        
        return XYdict

    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        for testname in self.testnames:
            self.figs['PRNU_vs_FLU_%s' % testname] = os.path.join(self.figspath,
                         'PRNU_vs_FLU_%s.png' % testname)
        
        for testname in self.testnames:
            self.figs['PRNU_MAP_%s' % testname] = os.path.join(self.figspath,
                         'PRNU_MAP_%s.png' % testname)
            
    
    def dump_aggregated_results(self):
        """ """
        
        
        # MASTER FLATS DISPLAY
        
        for testname in self.testnames:
            
            stestname = st.replace(testname,'_', '\_')
            
            for colkey in self.colkeys[testname]:
                
                MFdict = self._get_MFdict(testname, colkey)             
                MFkwargs = dict(suptitle='%s-%s: Master Flat-Field' % (stestname, colkey),
                                figname = self.figs['MF_%s_%s' % (testname, colkey)])
                
                self.plot_ImgFPA(MFdict, kwargs=MFkwargs)
        
        
        # PRNU vs. FLUENCE
        
        
        for testname in self.testnames:
            
            Ncols = self.Ncols[testname]
            stestname = st.replace(testname,'_', '\_')
                        
            XYdict = self._get_XYdict_PRNUFLU(self.ParsedTable[testname],Ncols)      
            
            self.plot_XY(XYdict,kwargs=dict(
                    title='%s: PRNU' % (stestname,),
                    doLegend=True,
                    xlabel='Fluence [ke-]',
                    ylabel='PRNU',
                    figname=self.figs['PRNU_vs_FLU_%s' % testname]))


        # PRNU maps
        
        for testname in self.testnames:
            
            Ncols = self.Ncols[testname]
            stestname = st.replace(testname,'_','\_')
            
            for icol in range(1,Ncols+1):
                
                colkey = 'col%03i' % icol
                
                PRNUMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], 
                                extractor=self._get_extractor_PRNU_fromPT(colkey))
            
                
                self.plot_SimpleMAP(PRNUMAP,kwargs=dict(
                        suptitle='%s [%s]: PRNU' % (stestname,colkey),
                        figname = self.figs['PRNU_MAP_%s' % testname]
                        ))
        

        
        
        
        