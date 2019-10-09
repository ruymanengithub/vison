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
        
        self.batches_highPRNU = ['14313','14471']
        
        self.products['HER_CURVES'] = OrderedDict()
        
        self.censored = []
        
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
        
        productspath = os.path.join(inventoryitem['resroot'],'products')
        her_pick = os.path.join(productspath,os.path.split(sidd.products['HER_PROFILES'])[-1])
        her_profs = files.cPickleRead(her_pick)['data'].copy()
        
        
        herprofkeys_v = np.zeros((1),dtype='S50')
        
        for iCCD, CCDk in enumerate(CCDkeys):
            
            herkey = '%s_%s_%s_%i_%s' % (testname, block, session,jrep+1,CCDk)
                
            self.products['HER_CURVES'][herkey] = her_profs.copy()
            
            herprofkeys_v[0] = herkey
        
        sidd.addColumn(herprofkeys_v, 'HERPROF_KEY', IndexS)

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

    def init_fignames(self):
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        self.figs['GvsT'] = os.path.join(self.figspath,
                         'GAIN_vs_Det_Temp.png')
        
        self.figs['GvsWave'] = os.path.join(self.figspath,
                         'GAIN_vs_Wavelength.png')
        
        self.figs['GvsOD'] = os.path.join(self.figspath,
                         'GAIN_vs_OD_CAL.png')
        
        self.figs['GvsRD'] = os.path.join(self.figspath,
                         'GAIN_vs_RD_CAL.png')
        
        for testname in self.testnames:
            self.figs['GAIN_MAP_%s' % testname] = os.path.join(self.figspath,
                         'GAIN_MAP_%s.png' % testname)
            
            self.figs['BLOOM_ADU_MAP_%s' % testname] = os.path.join(self.figspath,
                         'BLOOM_ADU_MAP_%s.png' % testname)
            
            self.figs['BLOOM_ELE_MAP_%s' % testname] = os.path.join(self.figspath,
                         'BLOOM_ELE_MAP_%s.png' % testname)
            
            self.figs['HER_MAP_%s' % testname] = os.path.join(self.figspath,
                         'HER_MAP_%s.png' % testname)
            
            self.figs['HER_curves_%s' % testname] = os.path.join(self.figspath,
                         'HER_curves_%s.png' % testname)

    
    def _get_XYdict_GvsLAM(self):
        """ """
        
        x = []
        y = []
        ey = []
        
        for testname in self.testnames:
            
            PT = self.ParsedTable[testname]
            
            x.append(PT['WAVENM'][0])
            
            _y = []
            
            for block in self.flight_blocks:
                
                ixblock = np.where(PT['BLOCK'] == block)
                
                for CCD in self.CCDs:
                                        
                    for Q in self.Quads:
                    
                        Gcol = 'GAIN_CCD%i_Quad%s' % (CCD,Q)
                        
                        _y.append(PT[Gcol][ixblock][0])
                        
            y.append(np.nanmean(_y))
            ey.append(np.nanstd(_y))
        
        sort =np.argsort(x)
        x = np.array(x)[sort]
        y = np.array(y)[sort]
        ey = np.array(ey)[sort]
                
        XYdict = dict(x=x,y=y,ey=ey)
        
        return XYdict
    
    def _get_XYdict_GvsT(self):
        """ """
        
        x = dict()
        y = dict()
        
        labelkeys = []
        
        for testname in self.testnames:
            labelkeys.append(testname)
            
            PT = self.ParsedTable[testname]
            
            x[testname] = []
            y[testname] = []
            
            for iCCD,CCD in enumerate(self.CCDs):
                TCCD = (PT['HK_CCD%s_TEMP_T' % CCD]+PT['HK_CCD%s_TEMP_B' % CCD])/2.
                for iQ, Q in enumerate(self.Quads):
                    
                    Ycolumn = 'GAIN_CCD%s_Quad%s' % (CCD,Q)
                    
                    _Ydata = PT[Ycolumn].copy()
                    
                    x[testname] += TCCD.tolist()
                    y[testname] += _Ydata.tolist()
                    
                    
            x[testname] = np.array(x[testname])
            y[testname] = np.array(y[testname])
        
        XYdict = dict(x=x,y=y,labelkeys=labelkeys)
        
        return XYdict

    def _get_XYdict_GvsOD(self):
        """ """
        
        x = dict()
        y = dict()
        
        labelkeys = []
        
        for testname in self.testnames:
            labelkeys.append(testname)
            
            PT = self.ParsedTable[testname]
            
            x[testname] = []
            y[testname] = []
            
            for iCCD,CCD in enumerate(self.CCDs):
                
                for iQ, Q in enumerate(self.Quads):
                    
                    if Q in ['E','F']:
                        Xcolumn = 'HK_CCD%i_OD_T_CAL' % (CCD,)
                    else:
                        Xcolumn = 'HK_CCD%i_OD_B_CAL' % (CCD,) 
                    
                    _Xdata = PT[Xcolumn].copy()
                    
                    Ycolumn = 'GAIN_CCD%s_Quad%s' % (CCD,Q)                    
                    _Ydata = PT[Ycolumn].copy()
                    
                    x[testname] += _Xdata.tolist()
                    y[testname] += _Ydata.tolist()
                    
                    
            x[testname] = np.array(x[testname])
            y[testname] = np.array(y[testname])
        
        XYdict = dict(x=x,y=y,labelkeys=labelkeys)
        
        return XYdict

    def _get_XYdict_GvsRD(self):
        """ """
        
        x = dict()
        y = dict()
        
        labelkeys = []
        
        for testname in self.testnames:
            labelkeys.append(testname)
            
            PT = self.ParsedTable[testname]
            
            x[testname] = []
            y[testname] = []
            
            for iCCD,CCD in enumerate(self.CCDs):
                
                for iQ, Q in enumerate(self.Quads):
                                        
                    if Q in ['E','F']:
                        Xcolumn = 'HK_COMM_RD_T_CAL' 
                    else:
                        Xcolumn = 'HK_COMM_RD_B_CAL' 
                    
                    _Xdata = PT[Xcolumn].copy()
                    
                    Ycolumn = 'GAIN_CCD%s_Quad%s' % (CCD,Q)                    
                    _Ydata = PT[Ycolumn].copy()
                    
                    x[testname] += _Xdata.tolist()
                    y[testname] += _Ydata.tolist()
                    
                    
            x[testname] = np.array(x[testname])
            y[testname] = np.array(y[testname])
        
        XYdict = dict(x=x,y=y,labelkeys=labelkeys)
        
        return XYdict

    def _get_XYdict_HER(self, testname):
        """ """
        
        x = dict()
        y = dict()
        
        PT = self.ParsedTable[testname]
        
        labelkeys = []
        
        for block in self.flight_blocks:
            ixsel = np.where(PT['BLOCK'] == block)
            herprof_key = PT['HERPROF_KEY'][ixsel][0]
            i_her = self.products['HER_CURVES'][herprof_key].copy()
        
            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' %  CCD
                
                for kQ, Q in enumerate(self.Quads):
                    
                    pkey = '%s_%s_%s' % (block, CCDk, Q)
                    
                    _x = i_her[CCDk][Q]['x'].copy()
                    _x -= _x.min()
                    _y = i_her[CCDk][Q]['y'].copy()
                    
                    
                    if pkey not in self.censored:
                        
                        x[pkey] = _x.copy()
                        y[pkey] = _y.copy()
                        labelkeys.append(pkey)
        
        
        HERdict = dict(x=x,y=y,labelkeys=labelkeys)
                
        return HERdict
    
    
    def dump_aggregated_results(self):
        """ """
        
        outpathroot = self.outpathroot
        
        doGainMaps=True
        doBloomMaps=True
        doHERMaps=True
        doHERcurves=True
        doGvsWave=True
        doGvsT=True
        doGvsOD=True
        doGvsRD=True
        
        # GAIN matrix (all blocks and test/waves) to dict() saved as pickle
        
        GAIN_MXdict = self.gen_GAIN_MXdict()
        
        GAIN_MXpick = os.path.join(outpathroot,'GAIN_MX_PTC0X.pick')
        
        files.cPickleDump(GAIN_MXdict,GAIN_MXpick)
        
        
        # GAIN maps (all tests/waves)
        
        if doGainMaps:
            
            for testname in self.testnames:
                
                GMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_GAIN_fromPT)
                
                avgG = self.get_stat_from_FPAMAP(GMAP,np.nanmean)
                print('Average G [%s]: %.2f' % (testname,avgG))
                
                stestname = st.replace(testname,'_','\_')
                
                self.plot_SimpleMAP(GMAP,kwargs=dict(
                        suptitle='%s: GAIN e-/ADU' % stestname,
                        figname=self.figs['GAIN_MAP_%s' % testname]))
        
        if doBloomMaps:

            # BLOOM maps (ADU and e-, from PTC01)
            
            for testname in self.testnames:
                
                BADU_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_BADU_fromPT)
            
                stestname = st.replace(testname,'_','\_')
                self.plot_SimpleMAP(BADU_MAP,kwargs=dict(
                        suptitle='%s: BLOOM-ADU [DN]' % stestname,
                        figname=self.figs['BLOOM_ADU_MAP_%s' % testname]))#,
                        #corekwargs=dict(norm = Normalize(vmin=3e4,vmax=2**16, clip=False))))
            
            for testname in self.testnames:
                
                BE_MAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_BE_fromPT)
            
                stestname = st.replace(testname,'_','\_')
                self.plot_SimpleMAP(BE_MAP,kwargs=dict(
                        suptitle='%s: BLOOM-ELECTRONS' % stestname,
                        figname=self.figs['BLOOM_ELE_MAP_%s' % testname]))#,
                        #corekwargs=dict(norm = Normalize(vmin=1e5,vmax=2.2E5, clip=False))))
        
        # HER map
        
        if doHERMaps:
        
            for testname in self.testnames:
                
                HERMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname], extractor=self._extract_HER_fromPT)
                
                stestname = st.replace(testname,'_','\_')
                            
                self.plot_SimpleMAP(HERMAP,kwargs=dict(
                        suptitle='%s: Hard Edge Response Factor' % stestname,
                        figname=self.figs['HER_MAP_%s' % testname]))
        
        # HER Curves
        
        if doHERcurves:
            
            for testname in self.testnames:
                
                HERSingledict = self._get_XYdict_HER(testname)
                
                stestname = st.replace(testname,'_','\_')
                                        
                self.plot_XY(HERSingledict,kwargs=dict(
                    title='%s: H.E.R. CURVES' % stestname,
                    doLegend=False,
                    xlabel='Pixel',
                    ylabel='HER [frac]',
                    ylim=[-2.E-4,5.e-4],
                    xlim=[9,15],
                    corekwargs=dict(linestyle='-',marker=''),
                    figname=self.figs['HER_curves_%s' % testname]))
                
                
        # GAIN vs. Wavelength
        
        if doGvsWave:
        
            GvsWavedict = self._get_XYdict_GvsLAM()
            
            self.plot_XY(GvsWavedict,kwargs=dict(
                        title='Gain vs. Wavelength',
                        doLegend=False,
                        doYErrbars=True,
                        xlabel='Wavelength',
                        ylabel='Gain [e-/ADU]',
                        ylim=[3.3,3.7],
                        corekwargs=dict(linestyle='',marker='.'),
                        #figname = ''))
                        figname=self.figs['GvsWave']))
        
        
        # GAIN vs. detector temperature (PTC01)
        
        if doGvsT:
        
            GvsTdict = self._get_XYdict_GvsT()
            
            self.plot_XY(GvsTdict,kwargs=dict(
                        title='Gain vs. Detector Temperature',
                        doLegend=True,
                        xlabel='Detector Temperature',
                        ylabel='Gain [e-/ADU]',                    
                        corekwargs=dict(linestyle='',marker='.'),                        
                        figname=self.figs['GvsT']))
        
        # GAIN vs. OD-CAL (PTC01)
        
        if doGvsOD:
        
            GvsODdict = self._get_XYdict_GvsOD()
            
            self.plot_XY(GvsODdict,kwargs=dict(
                        title='Gain vs. Output Drain',
                        doLegend=True,
                        xlabel='Output Drain (Calibrated) [V]',
                        ylabel='Gain [e-/ADU]',   
                        corekwargs=dict(linestyle='',marker='.'),
                        #figname=''))
                        figname=self.figs['GvsOD']))
        
#        # GAIN vs. RD-CAL (PTC01)
        
        if doGvsRD:        

            GvsRDdict = self._get_XYdict_GvsRD()
            
            self.plot_XY(GvsRDdict,kwargs=dict(
                        title='Gain vs. Reset Drain',
                        doLegend=True,
                        xlabel='Reset Drain (Calibrated) [V]',
                        ylabel='Gain [e-/ADU]',
                        corekwargs=dict(linestyle='',marker='.'),
                        #figname=''))
                        figname=self.figs['GvsRD']))
        
        # Save the ParsedTable(s)
        
        
        
        