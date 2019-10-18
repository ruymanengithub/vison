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
import matplotlib.cm as cm

from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files
from vison.xtalk import xtalk as xtalkmod

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
        self.products['XTALK_RT'] = files.cPickleRead(kwargs['cdps']['xtalk_roetab'])
        
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
        
        # XTALK-from ROE-TALK
        
        
        
        
        
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
    
    
    def _get_xtalk(self,xtalk_dict,mode='sign'):
            
            coefs = xtalk_dict['coefs']
            xsource = np.linspace(0, 2.**16, 200)
            yvictim = xtalkmod.f_fitXT(xsource, *coefs)

            ixmax = np.argmax(np.abs(yvictim))
            ixtalk = yvictim[ixmax]
            if mode == 'sign':
                return ixtalk
            elif mode == 'abs':
                return np.abs(ixtalk)
    
    
    def _get_XYdict_XT(self,TEST1, TEST2, mode='sign'):
        """ """
        
        x = dict()
        y = dict()
        
                
        PT1 = self.ParsedTable[TEST1]
        
        if TEST2 != 'RT':
            PT2 = self.ParsedTable[TEST2]
        
        
        for block in self.flight_blocks:
            
            _x = []
            _y = []
            
            ixblock1 = np.where(PT1['BLOCK'] == block)[0][0]
            ctalkkey1 = PT1['XTALK'][ixblock1]
            ctalk1 = self.products['XTALK'][ctalkkey1].copy()
            
            if TEST2 == 'RT':            
                ctalk2 = self.products['XTALK_RT'][block].copy()
            else:
                ixblock2 = np.where(PT2['BLOCK'] == block)[0][0]
                ctalkkey2 = PT2['XTALK'][ixblock2]
                ctalk2 = self.products['XTALK'][ctalkkey2].copy()
            
            for iCr, CCDref in enumerate(self.CCDs):
                for iQr, Qref in enumerate(self.Quads):
            
                    for iC, CCD in enumerate(self.CCDs):
                        for iQ, Q in enumerate(self.Quads):
                            CCDreftag = 'CCD%i' % CCDref
                            CCDtag = 'CCD%i' % CCD                            
                            try:
                                _x.append(self._get_xtalk(ctalk1[CCDreftag][Qref][CCDtag][Q],mode))
                                _y.append(self._get_xtalk(ctalk2[CCDreftag][Qref][CCDtag][Q],mode))
                            except KeyError:
                                pass
            x[block] = np.array(_x)
            y[block] = np.array(_y)
        
        x['oneone'] = [-100,100]
        y['oneone'] = [-100,100]
        
        labelkeys = self.flight_blocks + ['oneone']
        
        XTdict = dict(x=x,y=y,labelkeys=labelkeys)
            
        return XTdict
            
    
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
        
        self.figs['XTALK_MAP_RT'] = os.path.join(self.figspath,
                         'XTALK_MAP_RT.png')
        
        self.figs['XTALK_RTvs800'] = os.path.join(self.figspath,
                         'XTALK_RT_vs_800nm.png')
        
        self.figs['XTALK_RTvs800_ABS'] = os.path.join(self.figspath,
                         'XTALK_RT_vs_800nm_abs.png')
        
        for wave in [590,730,880]:
        
            self.figs['XTALK_%ivs800' % wave] = os.path.join(self.figspath,
                 'XTALK_%inm_vs_800nm.png' % wave)
        
        for testname in self.testnames:
            self.figs['XTALK_MAP_%s' % testname] = os.path.join(self.figspath,
                         'XTALK_MAP_%s.png' % testname)
    
    def dump_aggregated_results(self):
        """ """
        
        # XTALK MAP (ROE-TAB)
        
        XTALKs_RT = self.products['XTALK_RT'].copy()
        
        self.plot_XtalkMAP(XTALKs_RT,kwargs=dict(
                scale='ADU',
                showvalues=False,
                title='XTALK [ADU] - ROE-TAB',
                figname=self.figs['XTALK_MAP_RT']))
        
        # XTALK MAPS (OPTICAL)
        
        for testname in self.testnames:
            
            XTALKs = self.get_XTALKDICT_from_PT(testname)
            
            stestname = st.replace(testname,'_','\_')
            self.plot_XtalkMAP(XTALKs,kwargs=dict(
                    scale='ADU',
                    showvalues=False,
                    title='%s: XTALK [ADU]' % stestname,
                    figname=self.figs['XTALK_MAP_%s' % testname]))
        
        # XTALK: 800-optical vs. RT  (with SIGN)
        
        XT_RTvs800 = self._get_XYdict_XT('PSF01_800','RT',mode='sign')
        
        XTkwargs = dict(
                    title='Cross-Talk Comparison',
                    doLegend=False,
                    xlabel='Xtalk - Opt. 800nm',
                    ylabel='Xtalk - ROE-TAB',
                    xlim=[-20,50],
                    ylim=[-20,50],                    
                    figname=self.figs['XTALK_RTvs800'])
        
        BLOCKcolors = cm.rainbow(np.linspace(0,1,len(self.flight_blocks)))
        
        xtcorekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            xtcorekwargs['%s' % (block,)] = dict(linestyle='',
                       marker='.',color=jcolor)
        
        xtcorekwargs['oneone'] = dict(linestyle='--',marker='',color='k')
        
        XTkwargs['corekwargs'] = xtcorekwargs
        
        self.plot_XY(XT_RTvs800,kwargs=XTkwargs)
        
        
        # XTALK: 800-optical vs. RT (ABS-VALUE)
        
        XT_RTvs800_abs = self._get_XYdict_XT('PSF01_800','RT',mode='abs')
        
        XTABSkwargs = dict(
                    title='Cross-Talk Comparison - ABS. Value',
                    doLegend=False,
                    xlabel='Abs(Xtalk - Opt. 800 nm)',
                    ylabel='Abs(Xtalk - ROE-TAB)', 
                    xlim=[-20,50],
                    ylim=[-20,50],                               
                    figname=self.figs['XTALK_RTvs800_ABS'])
        
        
        XTABSkwargs['corekwargs'] = xtcorekwargs
        
        self.plot_XY(XT_RTvs800_abs,kwargs=XTABSkwargs)
        
        
        # XTALK: 800-optical vs. OTHER-opt (with SIGN)
        
        for wave in [590,730,880]:
        
            XT_NMvs800 = self._get_XYdict_XT('PSF01_800','PSF01_%i' % wave,mode='sign')
        
            XTNMvs800kwargs = dict(
                        title='Cross-Talk Comparison - With Sign',
                        doLegend=False,
                        xlabel='Xtalk - Opt. 800 nm',
                        ylabel='Xtalk - Opt. %i nm' % wave,
                        xlim=[-20,50],
                        ylim=[-20,50],
                        figname=self.figs['XTALK_%ivs800' % wave])
        
        
            XTNMvs800kwargs['corekwargs'] = xtcorekwargs
                            
            self.plot_XY(XT_NMvs800,kwargs=XTNMvs800kwargs)
        
        