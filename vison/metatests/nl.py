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
import matplotlib.cm as cm

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
        
        self.censored = ['NIELS_CCD2_G', 'BORN_CCD2_E', 
                         'OWEN_CCD1_E', 'OWEN_CCD1_F','OWEN_CCD1_G','OWEN_CCD1_H',
                         'OWEN_CCD2_E', 'OWEN_CCD2_F','OWEN_CCD2_G','OWEN_CCD2_H',
                         'OWEN_CCD3_E', 'OWEN_CCD3_F','OWEN_CCD3_G','OWEN_CCD3_H']
        
        self.init_fignames()
    
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
        
                
        NL_dict = sidd.products['NL'].copy()
 
        tmp_v_CQ = np.zeros((1,NCCDs,NQuads))
        
        maxNLpc_v = tmp_v_CQ.copy()
        fluADU_maxNLpc_v = tmp_v_CQ.copy()
        fluELE_maxNLpc_v = tmp_v_CQ.copy()        
        stabilitypc_v = tmp_v_CQ.copy()
        
        nlcdpkeys_v = np.zeros((1),dtype='S50')
        
        for iCCD, CCDk in enumerate(CCDkeys):
            
            nlcdpkey = '%s_%s_%s_%i_%s' % (testname, block, session,jrep+1,CCDk)
            
            for kQ, Q in enumerate(self.Quads):
                
                maxNLpc_v[0,iCCD,kQ] = NL_dict[CCDk][Q]['maxNLpc']
                
                gain = self.cdps['GAIN'][block][CCDk][Q][0]
                
                fluADU_maxNLpc_v[0,iCCD,kQ] = NL_dict[CCDk][Q]['flu_maxNLpc']
                fluELE_maxNLpc_v[0,iCCD,kQ] = gain * fluADU_maxNLpc_v[0,iCCD,kQ]
                
                stabilitypc_v[0,iCCD,kQ] = NL_dict[CCDk][Q]['stability_pc']
                
            self.products['NL'][nlcdpkey] = NL_dict.copy()
            
            nlcdpkeys_v[0] = nlcdpkey
            
                
        sidd.addColumn(maxNLpc_v, 'MAXNL', IndexCQ)
        sidd.addColumn(fluADU_maxNLpc_v, 'FLUADU_MAXNL', IndexCQ)
        sidd.addColumn(fluELE_maxNLpc_v, 'FLUELE_MAXNL', IndexCQ)
        sidd.addColumn(stabilitypc_v, 'STABILITYPC', IndexCQ)
        
        
        sidd.addColumn(nlcdpkeys_v, 'NLCDP_KEY', IndexS)

        # flatten sidd to table
        
        sit = sidd.flattentoTable()
        
        return sit
    
    def _get_NLMAP_from_PT(self,mode='fit'):
        """ """
                
        NLMAP = OrderedDict()
        NLMAP['labelkeys'] = self.Quads
        
        PT = self.ParsedTable['NL02']
        column = 'NLCDP_KEY'
        
        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                
                Ckey  = 'C_%i%i' % (jY+1,iX+1)
                NLMAP[Ckey] = OrderedDict()
                
                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                
                ixblock = np.where(PT['BLOCK'] == block)
                
                
                if len(ixblock[0])==0:
                    NLMAP[Ckey] = OrderedDict(x=OrderedDict(),
                                          y=OrderedDict())
                    for Q in self.Quads:
                        NLMAP[Ckey]['x'][Q] = []
                        NLMAP[Ckey]['y'][Q] = []
                    continue
                
                _nlkey = PT[column][ixblock][0]

                
                _nldict = self.products['NL'][_nlkey]
                
                _ccd_nldict = OrderedDict(x=OrderedDict(),
                                          y=OrderedDict())
                
                for Q in self.Quads:
                    
                    gain = self.cdps['GAIN'][block][CCDk][Q][0]
                    
                    if mode == 'fit':
                        _y = _nldict[CCDk][Q]['outputcurve']['Y'].copy()
                        _x = _nldict[CCDk][Q]['outputcurve']['X'] * gain / 1.E3
                    elif mode == 'data':
                        _y = _nldict[CCDk][Q]['inputcurve']['Y'].copy()
                        _x = _nldict[CCDk][Q]['inputcurve']['X'] * gain / 1.E3
                    
                    _ccd_nldict['x'][Q] = _x.copy()
                    _ccd_nldict['y'][Q] = _y.copy()
                                                
                NLMAP[Ckey] = _ccd_nldict.copy()
                
        
        return NLMAP
    
    def _get_XYdict_NL(self, mode='fit', scale='rel'):
        
        x = dict()
        y = dict()
        
        PT = self.ParsedTable['NL02']
        
        labelkeys = []
        
        
        if mode == 'data':
            curvekey = 'inputcurve'
        elif mode == 'fit':
            curvekey = 'outputcurve'
        
        for block in self.flight_blocks:
            ixsel = np.where(PT['BLOCK'] == block)
            nlcdp_key = PT['NLCDP_KEY'][ixsel][0]
            i_NL = self.products['NL'][nlcdp_key].copy()
        
            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' %  CCD
                
                for kQ, Q in enumerate(self.Quads):
                    
                    pkey = '%s_%s_%s' % (block, CCDk, Q)
                    
                    gain = self.cdps['GAIN'][block][CCDk][Q][0]
                    
                    xfluadu = i_NL[CCDk][Q][curvekey]['X'].copy()                                                
                    xfluKele = xfluadu * gain / 1.E3
                    
                    if pkey not in self.censored:
                        
                        x[pkey] = xfluKele.copy()
                        
                        if scale == 'rel':                       
                            y[pkey] = i_NL[CCDk][Q][curvekey]['Y'].copy()
                        elif scale == 'abs':                            
                            yrel = i_NL[CCDk][Q][curvekey]['Y'].copy()
                            y[pkey] = yrel / 100. * xfluKele*1.E3
                            
                        labelkeys.append(pkey)
        
        
        NLdict = dict(x=x,y=y,labelkeys=labelkeys)
                
        return NLdict
    
    def _get_compare_NLs_CCD(self, BLOCK1, BLOCK2, CCD):
        """ """
        
        Qdict = dict()
        
        PT = self.ParsedTable['NL02']
                
        
        #if mode == 'data':
        #    curvekey = 'inputcurve'
        #elif mode == 'fit':
        #    curvekey = 'outputcurve'
        
        
        for Q in self.Quads:
            Qdict[Q] = dict(x=dict(),y=dict())
        
        
        CCDk = 'CCD%i' % CCD
        
        labelkeys = []
        
        for block in [BLOCK1, BLOCK2]:
            
            ixblock = np.where(PT['BLOCK'] == block)
            
            timetag = '%s/%s/%s' % PT['time_CCD1'][ixblock][0].isocalendar()
            
            root_labelkey = '%s, %s' % (block[0:6], timetag)
            
            
            nlcdp_key = PT['NLCDP_KEY'][ixblock][0]
            i_NL = self.products['NL'][nlcdp_key].copy()
        
            for Q in self.Quads:
                
                gain = self.cdps['GAIN'][block][CCDk][Q][0]
                
                xfluadu_inp = i_NL[CCDk][Q]['inputcurve']['X'].copy()                                                
                xfluKele_inp = xfluadu_inp * gain / 1.E3
                y_inp = i_NL[CCDk][Q]['inputcurve']['Y'].copy()
                
                xfluadu_out = i_NL[CCDk][Q]['outputcurve']['X'].copy()                                                
                xfluKele_out = xfluadu_out * gain / 1.E3
                y_out = i_NL[CCDk][Q]['outputcurve']['Y'].copy()
                
                labelkeyin = '%s (data)' % root_labelkey
                
                Qdict[Q]['x'][labelkeyin] = xfluKele_inp.copy()
                Qdict[Q]['y'][labelkeyin] = y_inp.copy()
                
                labelkeyout = '%s (fit)' % root_labelkey
                
                Qdict[Q]['x'][labelkeyout] = xfluKele_out.copy()
                Qdict[Q]['y'][labelkeyout] = y_out.copy()
                
                if Q == self.Quads[0]:
                    labelkeys += [labelkeyin, labelkeyout]
                                    
        Qdict['labelkeys'] = labelkeys  
        
        return Qdict
    
    
    def init_fignames(self):
        """ """
        
        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)
        
        
        self.figs['NL_curves_MAP'] = os.path.join(self.figspath,
                         'NL_curves_FPAMAP.png')
        
        self.figs['NL_curves'] = os.path.join(self.figspath,
                         'NL_curves_rel_fit.png')
        
        self.figs['NL_curves_data'] = os.path.join(self.figspath,
                         'NL_curves_rel_data.png')
        
        self.figs['NL_curves_abs'] = os.path.join(self.figspath,
                         'NL_curves_abs_fit.png')
        
        for CCD in [1, 2]:
            self.figs['HvsK_CCD%s' % CCD] = os.path.join(self.figspath,
                         'NL_HEISEN_SKLODO_CCD%s.png' % CCD)
            
        for CCD in [1,3]:
            self.figs['JvsJ2_CCD%s' % CCD] = os.path.join(self.figspath,
                         'NL_JULES_JULES2_CCD%s.png' % CCD)

    def dump_aggregated_results(self):
        """ """
        
        # HeatMap of maximum non-linearities
        
        # PLOT All NL curves in Map-of-CCDs
        
        NLMAP = self._get_NLMAP_from_PT(mode='fit')
        
        self.plot_XYMAP(NLMAP,kwargs=dict(
                        suptitle='Non-Linearity Curves',
                        doLegend=True,
                        ylim = [-3.,7.],
                        corekwargs = dict(E=dict(linestyle='-',marker='',color='r'),
                                          F=dict(linestyle='-',marker='',color='g'),
                                          G=dict(linestyle='-',marker='',color='b'),
                                          H=dict(linestyle='-',marker='',color='m')),
                        figname = self.figs['NL_curves_MAP']
                        ))
        
        # PLOT All NL curves in single Plot - Relative, best fit curve
        
        NLSingledict = self._get_XYdict_NL(mode='fit',scale='rel')
        
        NLkwargs = dict(
                    title='NON-LINEARITY CURVES',
                    doLegend=False,
                    xlabel='Fluence [ke-]',
                    ylabel='Non-Linearity [pc]',
                    ylim=[-3.,7.],                    
                    figname=self.figs['NL_curves'])
        
        BLOCKcolors = cm.rainbow(np.linspace(0,1,len(self.flight_blocks)))
        
        linecorekwargs = dict()
        pointcorekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            for iCCD in self.CCDs:
                for kQ in self.Quads:
                    linecorekwargs['%s_CCD%i_%s' % (block,iCCD, kQ)] = dict(linestyle='-',
                               marker='',color=jcolor) 
                    pointcorekwargs['%s_CCD%i_%s' % (block,iCCD, kQ)] = dict(linestyle='',
                               marker=',',color=jcolor)
        
        NLkwargs['corekwargs'] = linecorekwargs
        
        self.plot_XY(NLSingledict,kwargs=NLkwargs)
        
        
        # PLOT All NL curves in single Plot - Relative, data points
        
        NLSingleDatadict = self._get_XYdict_NL(mode='data',scale='rel')
        
        NLDatakwargs = dict(
                    title='NON-LINEARITY DATA POINTS',
                    doLegend=False,
                    xlabel='Fluence [ke-]',
                    ylabel='Non-Linearity [pc]',
                    ylim=[-3.,7.],             
                    figname=self.figs['NL_curves_data'])
        
        
        NLDatakwargs['corekwargs'] = pointcorekwargs
        
        self.plot_XY(NLSingleDatadict,kwargs=NLDatakwargs)
        
        
        # PLOT All NL (best fit) curves in single Plot - Abs
        
        NLSingledict_abs = self._get_XYdict_NL(mode='fit',scale='abs')
        
        NLabskwargs = dict(
                    title='NON-LINEARITY CURVES - ABS. DEVIATION',
                    doLegend=False,
                    xlabel='Fluence [ke-]',
                    ylabel='Non-Linearity [$\Delta$ electrons]',
                    #ylim=[-3.,7.],                    
                    figname=self.figs['NL_curves_abs'])
        
        NLabskwargs['corekwargs'] = linecorekwargs
        
        self.plot_XY(NLSingledict_abs,kwargs=NLabskwargs)
        
        
        # Compare NL curves of same CCD in different Epochs (re-calibrations)
                
        # HEISENBERG CCD1, CCD2 vs. SKLODOWSKA CCD1, CCD2
        # JULES CCD1, CCD3 vs. JULES2 CCD1, CCD3
        
        
        for CCD in [1, 2]:
            
            XYCCD_HvsSK_CCDX = self._get_compare_NLs_CCD('HEISENBERG','SKLODOWSKA',CCD)
        
            HvsSK_kwargs = dict(
                        title='HEISENBERG/SKLODOWSKA, %s' % CCD,
                        doLegend=True,
                        xlabel='Fluence [ke-]',
                        ylabel='Non-Linearity [pc]',
                        ylim=[-2.,7.],                    
                        figname=self.figs['HvsK_CCD%s' % CCD])
            
            
            _labelkeys = XYCCD_HvsSK_CCDX['labelkeys']
            
            
            HvsSK_corekwargs = dict()
            for _labelkey in _labelkeys:
                if 'data' in _labelkey:
                    linestyle=''
                    marker='.'
                elif 'fit' in _labelkey:
                    linestyle='-'
                    marker=''
                if 'HEISEN' in _labelkey:
                    color = 'r'
                elif 'SLKODO' in _labelkey:
                    color='b'
                HvsSK_corekwargs[_labelkey] = \
                                dict(linestyle=linestyle,
                                     marker=marker,color=color)
                
            
            HvsSK_kwargs['corekwargs'] = HvsSK_corekwargs
            
            self.plot_XYCCD(XYCCD_HvsSK_CCDX, kwargs=HvsSK_kwargs)
        
        
        for CCD in [1, 3]:
            
            XYCCD_JvsJ2_CCDX = self._get_compare_NLs_CCD('JULES','JULES2',CCD)
            
            JvsJ2_kwargs = dict(
                        title='JULES/JULES2, CCD%s' % CCD,
                        doLegend=True,
                        xlabel='Fluence [ke-]',
                        ylabel='Non-Linearity [pc]',  
                        ylim=[-2.,7.],     
                        figname=self.figs['JvsJ2_CCD%s' % CCD])
            
            
            _labelkeys = XYCCD_JvsJ2_CCDX['labelkeys']
            
            JvsJ2_corekwargs = dict()
            for _labelkey in _labelkeys:
                if 'data' in _labelkey:
                    linestyle=''
                    marker='.'
                elif 'fit' in _labelkey:
                    linestyle='-'
                    marker=''
                if 'JULES,' in _labelkey:
                    color = 'r'
                elif 'JULES2,' in _labelkey:
                    color='b'
                JvsJ2_corekwargs[_labelkey] = \
                                dict(linestyle=linestyle,
                                     marker=marker,color=color)
            JvsJ2_kwargs['corekwargs'] = JvsJ2_corekwargs
            
            self.plot_XYCCD(XYCCD_JvsJ2_CCDX, kwargs=JvsJ2_kwargs)
        
        