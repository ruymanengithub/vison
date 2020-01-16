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

from vison.datamodel import cdp
from vison.support import utils
from vison.support import files
from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse
#from vison.support import vjson

import matplotlib.cm as cm
from matplotlib import pyplot as plt
plt.switch_backend('TkAgg')
from matplotlib.colors import Normalize
# END IMPORT

cols2keep = [
    'test',
    'sn_ccd1',
    'sn_ccd2',
    'sn_ccd3',
    'sn_roe',
    'sn_rpsu',
    'exptime',
    'vstart',
    'vend',
    'rdmode',
    'flushes',
    'siflsh',
    'siflsh_p',
    'swellw',
    'swelldly',
    'inisweep',
    'cdpu_clk',
    'chinj',
    'chinj_on',
    'chinj_of',
    'id_wid',
    'id_dly',
    'chin_dly',
    'v_tpump',
    's_tpump',
    'v_tp_mod',
    's_tp_mod',
    'v_tp_cnt',
    's_tp_cnt',
    'dwell_v',
    'dwell_s',
    'toi_fl',
    'toi_tp',
    'toi_ro',
    'toi_ch',
    'motr',
    'motr_cnt',
    'motr_siz',
    'source',
    'wave',
    'mirr_on',
    'mirr_pos',
    'R1C1_TT',
    'R1C1_TB',
    'R1C2_TT',
    'R1C2_TB',
    'R1C3_TT',
    'R1C3_TB',
    'IDL',
    'IDH',
    'IG1_1_T',
    'IG1_2_T',
    'IG1_3_T',
    'IG1_1_B',
    'IG1_2_B',
    'IG1_3_B',
    'IG2_T',
    'IG2_B',
    'OD_1_T',
    'OD_2_T',
    'OD_3_T',
    'OD_1_B',
    'OD_2_B',
    'OD_3_B',
    'RD_T',
    'RD_B',
    'time',
    'HK_CCD1_TEMP_T',
    'HK_CCD2_TEMP_T',
    'HK_CCD3_TEMP_T',
    'HK_CCD1_TEMP_B',
    'HK_CCD2_TEMP_B',
    'HK_CCD3_TEMP_B',
    'HK_CCD1_OD_T',
    'HK_CCD2_OD_T',
    'HK_CCD3_OD_T',
    'HK_CCD1_OD_B',
    'HK_CCD2_OD_B',
    'HK_CCD3_OD_B',
    'HK_COMM_RD_T',
    'HK_COMM_RD_B',
    'HK_CCD1_IG1_T',
    'HK_CCD2_IG1_T',
    'HK_CCD3_IG1_T',
    'HK_CCD1_IG1_B',
    'HK_CCD2_IG1_B',
    'HK_CCD3_IG1_B',
    'HK_COMM_IG2_T',
    'HK_COMM_IG2_B',
    'HK_FPGA_BIAS_ID2',
    'HK_VID_PCB_TEMP_T',
    'HK_VID_PCB_TEMP_B',
    'HK_RPSU_TEMP1',
    'HK_FPGA_PCB_TEMP_T',
    'HK_FPGA_PCB_TEMP_B',
    'HK_RPSU_TEMP_2',
    'HK_RPSU_28V_PRI_I',
    'chk_NPIXOFF',
    'chk_NPIXSAT',
    'offset_pre',
    'offset_img',
    'offset_ove',
    'std_pre',
    'std_img',
    'std_ove',
    'RON']


class MetaBias(MetaCal):
    """ """

    def __init__(self, **kwargs):
        """ """

        super(MetaBias, self).__init__(**kwargs)

        self.testnames = ['BIAS01', 'BIAS02']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['MB_PROFILES'] = OrderedDict()

        self.init_fignames()
        self.init_outcdpnames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """

        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']

        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0]),
                                     vcore.vIndex('CCD', vals=self.CCDs),
                                     vcore.vIndex('Quad', vals=self.Quads)])

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

        test_v = np.array([jrep + 1])
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)

        test_v = np.array([session])
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)

        test_v = np.array([testname])
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)

        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        off_pre_v = tmp_v_CQ.copy()
        off_img_v = tmp_v_CQ.copy()
        off_ove_v = tmp_v_CQ.copy()

        ron_pre_v = tmp_v_CQ.copy()
        ron_img_v = tmp_v_CQ.copy()
        ron_ove_v = tmp_v_CQ.copy()

        productspath = os.path.join(inventoryitem['resroot'], 'products')
        profspath = os.path.join(inventoryitem['resroot'], 'profiles')

        roncdp_pick = os.path.join(productspath, os.path.split(sidd.products['RON_CDP'])[-1])
        roncdp = files.cPickleRead(roncdp_pick)

        offcdp_pick = os.path.join(productspath, os.path.split(sidd.products['OFF_CDP'])[-1])
        offcdp = files.cPickleRead(offcdp_pick)

        profs_pick = os.path.join(profspath, os.path.split(sidd.products['MB_PROFILES'])[-1])
        profs_dict = files.cPickleRead(profs_pick)

        profkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)

        self.products['MB_PROFILES'][profkey] = profs_dict.copy()

        profskeys_v = np.zeros((1),dtype='S50')
        
        profskeys_v[0] = profkey

        for iCCD, CCDk in enumerate(CCDkeys):

            for kQ, Q in enumerate(self.Quads):

                off_pre_v[0, iCCD, kQ] = offcdp['data']['OFF_PRE'][CCDk][Q]
                off_img_v[0, iCCD, kQ] = offcdp['data']['OFF_IMG'][CCDk][Q]
                off_ove_v[0, iCCD, kQ] = offcdp['data']['OFF_OVE'][CCDk][Q]

                ron_pre_v[0, iCCD, kQ] = roncdp['data']['RON_PRE'][CCDk][Q]
                ron_img_v[0, iCCD, kQ] = roncdp['data']['RON_IMG'][CCDk][Q]
                ron_ove_v[0, iCCD, kQ] = roncdp['data']['RON_OVE'][CCDk][Q]
                

        sidd.addColumn(off_pre_v, 'OFF_PRE', IndexCQ)
        sidd.addColumn(off_img_v, 'OFF_IMG', IndexCQ)
        sidd.addColumn(off_ove_v, 'OFF_OVE', IndexCQ)

        sidd.addColumn(ron_pre_v, 'RON_PRE', IndexCQ)
        sidd.addColumn(ron_img_v, 'RON_IMG', IndexCQ)
        sidd.addColumn(ron_ove_v, 'RON_OVE', IndexCQ)

        sidd.addColumn(profskeys_v, 'MBPROFS_KEY', IndexS)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _get_extractor_RON_fromPT(self, units, reg):
        """ """

        def _extract_RON_fromPT(PT, block, CCDk, Q):
            
            ixblock = self.get_ixblock(PT, block)
            
            if units == 'ADU':
                unitsConvFactor = 1
            elif units == 'E':
                unitsConvFactor = self.cdps['GAIN'][block][CCDk][Q][0]
            
            if reg != 'all':            
                column = 'RON_%s_%s_Quad%s' % (reg.upper(),CCDk, Q)
                RON = np.nanmedian(PT[column][ixblock]) * unitsConvFactor
            else:
                RON=OrderedDict()
                for _reg in ['pre','img','ove']:
                    column = 'RON_%s_%s_Quad%s' % (_reg.upper(),CCDk, Q)
                    RON[_reg] = np.nanmedian(PT[column][ixblock]) * unitsConvFactor
            
            return RON

        return _extract_RON_fromPT
    
    def _get_extractor_OFFSET_fromPT(self, reg):
        """ """
    
        def _extractor(PT, block, CCDk, Q):
            """ """
            ixblock = self.get_ixblock(PT, block)
            
            if reg != 'all':            
                column = 'OFF_%s_%s_Quad%s' % (reg.upper(),CCDk, Q)
                OFF = np.nanmedian(PT[column][ixblock])
            else:
                OFF=OrderedDict()
                for _reg in ['pre','img','ove']:
                    column = 'OFF_%s_%s_Quad%s' % (_reg.upper(),CCDk, Q)
                    OFF[_reg] = np.nanmedian(PT[column][ixblock])
            return OFF
        
        return _extractor

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        for testname in self.testnames:
            self.figs['RON_ADU_%s' % testname] = os.path.join(self.figspath,
                                                              '%s_RON_ADU_MAP.png' % testname)
            self.figs['RON_ELE_%s' % testname] = os.path.join(self.figspath,
                                                              '%s_RON_ELE_MAP.png' % testname)
            self.figs['OFFSETS_%s' % testname] = os.path.join(self.figspath,
                                                              '%s_OFFSETS_MAP.png' % testname)
        profdirs = ['hor','ver']

        for testname in self.testnames:
            for profdir in profdirs:
                self.figs['PROFS_%s_%s' % (testname, profdir)] =  os.path.join(self.figspath,
                                    '%s_PROFS_%s.png' % (testname,profdir))


    def init_outcdpnames(self):
        """ """

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        for testname in self.testnames:
            self.outcdps['RON_ADU_%s' % testname] = '%s_RON_ADU_MAP.json' % testname          
            self.outcdps['OFFSETS_%s' % testname] = '%s_OFFSETS_MAP.json' % testname

    def _get_XYdict_PROFS(self,test,proftype):
        """ """

        x = dict()
        y = dict()

        labelkeys = []

        PT = self.ParsedTable[test]

        for block in self.flight_blocks:

            ixsel = np.where(PT['BLOCK'] == block)
            prof_key = PT['MBPROFS_KEY'][ixsel][0]

            i_Prof = self.products['MB_PROFILES'][prof_key].copy()

            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' % CCD

                for kQ, Q in enumerate(self.Quads):

                    pkey = '%s_%s_%s' % (block, CCDk, Q)

                    _pcq = i_Prof['data'][proftype][CCDk][Q].copy()

                    x[pkey] = _pcq['x'].copy()
                    y[pkey] = _pcq['y'] - np.nanmedian(_pcq['y'])

                    labelkeys.append(pkey)

        Pdict = dict(x=x,y=y,labelkeys=labelkeys)

        return Pdict

    

    def dump_aggregated_results(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='dump', Title='Aggregated Results', level=0)

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        # RON maps (all tests/waves)

        # RON maps, ADUs
        for testname in self.testnames:

            RONADUMAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_extractor_RON_fromPT(
                    units='ADU',
                    reg='all'))
            
            rn_header = OrderedDict()
            rn_header['title'] = 'RON MAP'
            rn_header.update(CDP_header)
            
            rn_cdp = cdp.Json_CDP(rootname=self.outcdps['RON_ADU_%s' % testname],
                              path=self.cdpspath)
            rn_cdp.ingest_inputs(data=RONADUMAP,
                             header = rn_header,
                             meta=dict(units='ADU'))
            rn_cdp.savehardcopy()
            
            RON_OVE_ADUMAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_extractor_RON_fromPT(
                    units='ADU',
                    reg='ove'))

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(RON_OVE_ADUMAP, **dict(
                suptitle='%s [OVERSCAN]: RON [ADU]' % stestname,
                figname=self.figs['RON_ADU_%s' % testname]))

        # RON maps, ELECTRONs
        for testname in self.testnames:

            RONEMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname],
                                              extractor=self._get_extractor_RON_fromPT(units='E',
                                                                                       reg='ove'))

            stestname = st.replace(testname, '_', '\_')
            figkey = 'RON_ELE_%s' % testname
            figname = self.figs[figkey]

            self.plot_SimpleMAP(RONEMAP, **dict(
                suptitle='%s: RON [ELECTRONS]' % stestname,
                figname=figname))

            self.addFigure2Report(figname, 
                    figkey=figkey, 
                    caption='%s: RON in e- (rms).' % stestname, 
                    texfraction=0.7)
            
            # reporting tables of RON-e to report

            if self.report is not None:

                eroncdpdict = dict(
                    caption='%s: RON (e-, rms) in the over-scan region.' % testname,
                    valformat='%.2f')

                ignore = self.add_StdQuadsTable2Report( 
                                Matrix = RONEMAP,
                                cdpdict = eroncdpdict)



        # OFFSET maps

        for testname in self.testnames:

            OFFMAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_extractor_OFFSET_fromPT(reg='all'))
            
            off_header = OrderedDict()
            off_header['title'] = 'OFFSETS MAP'
            off_header.update(CDP_header)
            
            off_cdp = cdp.Json_CDP(rootname=self.outcdps['OFFSETS_%s' % testname],
                              path=self.cdpspath)
            off_cdp.ingest_inputs(data=OFFMAP,
                             header = off_header,
                             meta=dict(units='ADU'))
            off_cdp.savehardcopy()
            
            OFF_OVE_MAP = self.get_FPAMAP_from_PT(
                self.ParsedTable[testname],
                extractor=self._get_extractor_OFFSET_fromPT(reg='ove'))
            
            figkey = 'OFFSETS_%s' % testname

            stestname = st.replace(testname, '_', '\_')
            self.plot_SimpleMAP(OFF_OVE_MAP, **dict(
                suptitle='%s [OVERSCAN]: OFFSET' % stestname,
                figname=self.figs[figkey]))

            self.addFigure2Report(self.figs[figkey], 
                    figkey=figkey, 
                    caption='%s: Avg. offsets in ADU. Measured in serial over-scan.' % stestname, 
                    texfraction=0.7)

            if self.report is not None:

                offcdpdict = dict(
                    caption='%s: Offset levels [ADU].' % testname,
                    valformat='%.2f')
                
                def _getOFF_val(Ckey, Q):
                    return OFFMAP[Ckey][Q]['ove']

                ignore = self.add_StdQuadsTable2Report( 
                                extractor = _getOFF_val,
                                cdpdict = offcdpdict)

        # Vertical and Horizonal average profiles

        xlabels_profs = dict(hor='column [pix]',
                            ver='row [pix]')

        proftypes = ['hor','ver']

        BLOCKcolors = cm.rainbow(np.linspace(0, 1, len(self.flight_blocks)))

        pointcorekwargs = dict()
        for jblock, block in enumerate(self.flight_blocks):
            jcolor = BLOCKcolors[jblock]
            for iCCD in self.CCDs:
                for kQ in self.Quads:
                    pointcorekwargs['%s_CCD%i_%s' % (block, iCCD, kQ)] = dict(
                        linestyle='', marker='.', color=jcolor, ms=0.8)

        for testname in self.testnames:

            for proftype in proftypes:

                XY_profs = self._get_XYdict_PROFS(test=testname,proftype=proftype)

                if proftype=='hor':
                    xlim=[0,100]
                else:
                    xlim=None

                figkey = 'PROFS_%s_%s' % (testname,proftype)

                profkwargs = dict(
                    title='%s: Avg. profiles, direction: %s' % (testname, proftype),
                    doLegend=False,
                    xlabel=xlabels_profs[proftype],
                    ylabel=r'$\delta ADU$',
                    ylim=[-20,20],
                    xlim=xlim,
                    figname=self.figs[figkey],
                    corekwargs=pointcorekwargs)

                self.plot_XY(XY_profs, **profkwargs)

                if proftype == 'ver':
                    captemp = '%s: Stacked profiles of Master Bias quadrant'+\
                    ' images in "parallel" direction. Median value of profile has been subtracted for '+\
                    'clarity. Each colour corresponds to a different block (each with 3x4 quadrants).'
                elif proftype == 'hor':
                    captemp = '%s: Stacked profiles of Master Bias quadrant'+\
                    ' images in "serial" direction. Median value of profile has been subtracted for '+\
                    'clarity. Only the first 100 columns are shown, because that is where most of the '+\
                    'intersting structure of this type of profile is. Each colour corresponds to a'+\
                    ' different block (each with 3x4 quadrants).'

                self.addFigure2Report(self.figs[figkey], 
                    figkey=figkey, 
                    caption= captemp % (stestname,),
                    texfraction=0.7)