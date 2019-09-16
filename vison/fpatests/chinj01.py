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
from vison.inject import lib as ilib


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
         'offset_pre', 'offset_ove', 'std_pre', 'std_ove']
         
class MetaChinj01(MetaCal):
    """ """
    
    def __init__(self, **kwargs):
        """ """
        
        super(MetaChinj01,self).__init__(**kwargs)
        
        self.testnames = ['CHINJ01']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        
        allgains = files.cPickleRead(kwargs['cdps']['gain'])
        
        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy() 
        
        self.products['METAFIT'] = OrderedDict()
        
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
        
        productspath = os.path.join(inventoryitem['resroot'],'products')
        
        
        metafitcdp_pick = os.path.join(productspath,os.path.split(sidd.products['METAFIT_CDP'])[-1])
        metafitcdp = files.cPickleRead(metafitcdp_pick)
        metafit = copy.deepcopy(metafitcdp['data']['ANALYSIS'])
        
        metafitkey = '%s_%s_%s_%i' % (testname, block, session,jrep+1)
        self.products['METAFIT'][metafitkey] = copy.deepcopy(metafit)
        metafitkey_v = np.array([metafitkey])
        sidd.addColumn(metafitkey_v, 'METAFIT', IndexS, ix=0)
        
        metacdp_pick = os.path.join(productspath,os.path.split(sidd.products['META_CDP'])[-1]) # change to META_CDP
        metacdp = files.cPickleRead(metacdp_pick)
        meta = metacdp['data']['ANALYSIS'] # this is a pandas DataFrame
        
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
                
        bgd_adu_v = tmp_v_CQ.copy()
        ig1_thresh_v = tmp_v_CQ.copy()
        ig1_notch_v = tmp_v_CQ.copy()
        slope_v = tmp_v_CQ.copy()
        n_adu_v = tmp_v_CQ.copy()        
        
        
        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):
                
                ixloc = np.where((meta['CCD'] == iCCD+1) & (meta['Q'] == kQ+1))
                
                bgd_adu_v[0,iCCD,kQ] = meta['BGD_ADU'][ixloc[0][0]]
                ig1_thresh_v[0,iCCD,kQ] = meta['IG1_THRESH'][ixloc[0][0]]
                ig1_notch_v[0,iCCD,kQ] = meta['IG1_NOTCH'][ixloc[0][0]]
                slope_v[0,iCCD,kQ] = meta['S'][ixloc[0][0]]
                n_adu_v[0,iCCD,kQ] = meta['N_ADU'][ixloc[0][0]]
                       
        
        sidd.addColumn(bgd_adu_v, 'FIT_BGD_ADU', IndexCQ)
        sidd.addColumn(ig1_thresh_v, 'FIT_IG1_THRESH', IndexCQ)
        sidd.addColumn(ig1_notch_v, 'FIT_IG1_NOTCH', IndexCQ)        
        sidd.addColumn(slope_v, 'FIT_SLOPE', IndexCQ)
        sidd.addColumn(n_adu_v, 'FIT_N_ADU', IndexCQ)
        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit

    def _get_extractor_NOTCH_fromPT(self,units):
        """ """
        
        def _extract_NOTCH_fromPT(PT, block, CCDk, Q):
            
            ixblock = np.where(PT['BLOCK'].data == block)
            column = 'FIT_N_ADU_%s_Quad%s' % (CCDk,Q)
            
            if units =='ADU':
                unitsConvFactor=1
            elif units == 'E':
                unitsConvFactor = self.cdps['GAIN'][block][CCDk][Q][0]
            
            Notch = np.nanmedian(PT[column][ixblock]) * unitsConvFactor
            return Notch
        
        return _extract_NOTCH_fromPT
    
    def _get_injcurve(self, _chfitdf, ixCCD, ixQ, IG1raw, gain):
        """ """
        ixsel = np.where((_chfitdf['CCD'] == ixCCD) & (_chfitdf['Q'] == ixQ))
        
        pars = ['BGD','K','XT','XN','A','N']
        trans = dict(BGD='b',K='k',XT='xt', XN='xN', A='a', N='N')
        
        parsdict = dict()
        for par in pars:
            parsdict[trans[par]] = _chfitdf[par].as_matrix()[ixsel][0]
        
        parsdict['IG1'] = IG1raw.copy()
        
        inj = ilib.f_Inj_vs_IG1_ReLU(**parsdict) * 2.**16 # ADU
        
        inj_kel = inj * gain / 1.E3
        
        return inj_kel

    def _get_CHIG1_MAP_from_PT(self,kind='CAL'):
        """ """
        
        CHIG1MAP = OrderedDict()
        CHIG1MAP['labelkeys'] = self.Quads
        
        PT = self.ParsedTable['CHINJ01']
        column = 'METAFIT'
        
        IG1s = [2.5, 6.75]
        dIG1 = 0.05
        
        NIG1 = (IG1s[1]-IG1s[0])/dIG1+1
        IG1raw = np.arange(NIG1)*dIG1+IG1s[0]
        
        
        
        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                
                Ckey  = 'C_%i%i' % (jY+1,iX+1)
                CHIG1MAP[Ckey] = OrderedDict()
                
                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                
                jCCD = int(CCDk[-1])
                
                ixblock = np.where(PT['BLOCK'] == block)
                
                
                if len(ixblock[0])==0:
                    CHIG1MAP[Ckey] = OrderedDict(x=OrderedDict(),
                                          y=OrderedDict())
                    for Q in self.Quads:
                        CHIG1MAP[Ckey]['x'][Q] = []
                        CHIG1MAP[Ckey]['y'][Q] = []
                    continue
                
                _chkey = PT[column][ixblock][0]
                
                _chfitdf = self.products['METAFIT'][_chkey]
                
                _ccd_chfitdict = OrderedDict(x=OrderedDict(),
                                          y=OrderedDict())
                
                for kQ, Q in enumerate(self.Quads):
                    
                    roeVCal = self.roeVCals[block]
                
                    IG1cal = roeVCal.fcal_HK(IG1raw, 'IG1', jCCD, Q)
                                        
                    gain = self.cdps['GAIN'][block][CCDk][Q][0]
                    
                    inj_kel = self._get_injcurve(_chfitdf, jCCD, kQ+1, IG1raw, gain)
                    
                    if kind == 'CAL':
                        _IG1 = IG1cal.copy()
                    elif kind == 'RAW':
                        _IG1 = IG1raw.copy()
                    
                    _ccd_chfitdict['x'][Q] = _IG1.copy()
                    _ccd_chfitdict['y'][Q] = inj_kel.copy()
                
                CHIG1MAP[Ckey] = _ccd_chfitdict.copy()
        
        return CHIG1MAP

    def _get_XYdict_INJ(self, kind='CAL'):
        
        x = dict()
        y = dict()
        
        PT = self.ParsedTable['CHINJ01']
        column = 'METAFIT'
        
        IG1s = [2.5, 6.75]
        dIG1 = 0.05
        
        NIG1 = (IG1s[1]-IG1s[0])/dIG1+1
        IG1raw = np.arange(NIG1)*dIG1+IG1s[0]
        
        labelkeys = []
        
        for block in self.flight_blocks:
            ixblock = np.where(PT['BLOCK'] == block)
            ch_key = PT[column][ixblock][0]
            chfitdf = self.products['METAFIT'][ch_key]
            
        
            for iCCD, CCD in enumerate(self.CCDs):
                
                CCDk = 'CCD%i' %  CCD
                
                for kQ, Q in enumerate(self.Quads):
                    
                    roeVCal = self.roeVCals[block]
                
                    IG1cal = roeVCal.fcal_HK(IG1raw, 'IG1', iCCD+1, Q)
                    
                    gain = self.cdps['GAIN'][block][CCDk][Q][0]
                    
                    if kind == 'CAL':
                        _IG1 = IG1cal.copy()
                    elif kind == 'RAW':
                        _IG1 = IG1raw.copy()
                    
                    pkey = '%s_%s_%s' % (block, CCDk, Q)
                    
                    inj_kel = self._get_injcurve(chfitdf, iCCD+1, kQ+1, IG1raw, gain)
                    
                    x[pkey] = _IG1.copy()
                    y[pkey] = inj_kel.copy()
                    labelkeys.append(pkey)
        
        
        CHdict = dict(x=x,y=y,labelkeys=labelkeys)
                
        return CHdict
    


    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
            
        self.figs['NOTCH_ADU_MAP'] = os.path.join(self.figspath,
                     'NOTCH_ADU_MAP.png')
        
        self.figs['NOTCH_ELE_MAP'] = os.path.join(self.figspath,
                     'NOTCH_ELE_MAP.png')
        
        self.figs['CHINJ01_curves_IG1_RAW'] = os.path.join(self.figspath,
                     'CHINJ01_CURVES_IG1_RAW.png')
        
        self.figs['CHINJ01_curves_IG1_CAL'] = os.path.join(self.figspath,
                     'CHINJ01_CURVES_IG1_CAL.png')
        
        
        self.figs['CHINJ01_curves_MAP_IG1_CAL'] = os.path.join(self.figspath,
                     'CHINJ01_CURVES_MAP_IG1_CAL.png')
        

    def dump_aggregated_results(self):
        """ """
        
        
        # Histogram of Slopes [ADU/electrons]
        
        # Histogram of Notch [ADU/electrons]

        # Histogram of IG1_THRESH
        
        # Injection level vs. Calibrated IG1, all channels
        
        CURVES_IG1CAL_MAP = self._get_CHIG1_MAP_from_PT(kind='CAL')
        
        self.plot_XYMAP(CURVES_IG1CAL_MAP,kwargs=dict(
                        suptitle='Charge Injection Curves - Calibrated IG1',
                        doLegend=True,
                        ylabel='Inj [kel]',
                        xlabel = 'IG1 [V]',
                        corekwargs = dict(E=dict(linestyle='-',marker='',color='r'),
                                          F=dict(linestyle='-',marker='',color='g'),
                                          G=dict(linestyle='-',marker='',color='b'),
                                          H=dict(linestyle='-',marker='',color='m')),
                        figname = self.figs['CHINJ01_curves_MAP_IG1_CAL']
                        ))
        
        IG1CAL_Singledict = self._get_XYdict_INJ(kind='CAL')
        
        IG1CAL_kwargs = dict(
                    title='Charge Injection Curves - Calibrated IG1',
                    doLegend=False,
                    xlabel='IG1 (Calibrated) [V]',
                    ylabel='Injection [kel]',                    
                    figname=self.figs['CHINJ01_curves_IG1_CAL'])
        
        corekwargs = dict()
        for block in self.flight_blocks:
            for iCCD in self.CCDs:
                corekwargs['%s_CCD%i_E' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#FF4600') # red
                corekwargs['%s_CCD%i_F' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#61FF00') # green
                corekwargs['%s_CCD%i_G' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#00FFE0') # cyan
                corekwargs['%s_CCD%i_H' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#1700FF') # blue
            
        
        IG1CAL_kwargs['corekwargs'] = corekwargs.copy()
        
        self.plot_XY(IG1CAL_Singledict,kwargs=IG1CAL_kwargs)

        IG1RAW_Singledict = self._get_XYdict_INJ(kind='RAW')
        
        IG1RAW_kwargs = dict(
                    title='Charge Injection Curves - RAW IG1',
                    doLegend=False,
                    xlabel='IG1 (RAW) [V]',
                    ylabel='Injection [kel]',                    
                    figname=self.figs['CHINJ01_curves_IG1_RAW'])
        
        corekwargs = dict()
        for block in self.flight_blocks:
            for iCCD in self.CCDs:
                corekwargs['%s_CCD%i_E' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#FF4600') # red
                corekwargs['%s_CCD%i_F' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#61FF00') # green
                corekwargs['%s_CCD%i_G' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#00FFE0') # cyan
                corekwargs['%s_CCD%i_H' % (block,iCCD)] = dict(linestyle='-',
                           marker='',color='#1700FF') # blue
        
        
        IG1RAW_kwargs['corekwargs'] = corekwargs.copy()
        
        self.plot_XY(IG1RAW_Singledict,kwargs=IG1RAW_kwargs)
        
        # Notch level vs. calibrated IG2
        
        # Notch level vs. calibrated IDL
        
        # Notch level vs. calibrated OD
        
        
        
        # Notch injection map, ADUs
            
        NOTCHADUMAP = self.get_FPAMAP_from_PT(self.ParsedTable['CHINJ01'], 
                                            extractor=self._get_extractor_NOTCH_fromPT(units='ADU'))
        
        
        self.plot_SimpleMAP(NOTCHADUMAP,kwargs=dict(
                suptitle='CHINJ01: NOTCH INJECTION [ADU]',
                figname=self.figs['NOTCH_ADU_MAP']))
        
        # Notch injection map, ELECTRONs
        
        NOTCHEMAP = self.get_FPAMAP_from_PT(self.ParsedTable['CHINJ01'], 
                                            extractor=self._get_extractor_NOTCH_fromPT(units='E'))
        
        self.plot_SimpleMAP(NOTCHEMAP,kwargs=dict(
                suptitle='CHINJ01: NOTCH INJECTION [ELECTRONS]',
                figname=self.figs['NOTCH_ELE_MAP']))
        
        
        
        
        
        
        