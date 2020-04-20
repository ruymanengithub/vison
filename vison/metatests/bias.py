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
from scipy import ndimage
from skimage import exposure
import pandas as pd

from vison.datamodel import cdp
from vison.datamodel import ccd as ccdmod
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
        self.products['MASTERBIAS'] = OrderedDict()

        self.init_fignames()
        self.init_outcdpnames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """

        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']

        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        IndexC = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0]),
                                    vcore.vIndex('CCD', vals=self.CCDs)])

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


        # ADDING REFERENCES TO MASTER BIAS

        tmp_v_C = np.zeros((1, NCCDs), dtype='S50')
        
        mbkey_v = tmp_v_C.copy()

        for jCCD, CCDk in enumerate(CCDkeys):

            mbkey = '%s_%s_%s_%i_%s' % (testname, block, session, jrep + 1, CCDk)

            mbkey_v[0, jCCD] = mbkey

            self.products['MASTERBIAS'][mbkey] = os.path.split(
                sidd.products['MASTERBIAS_%s' % CCDk])[-1]

        sidd.addColumn(mbkey_v, 'MASTERBIAS', IndexC)


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
        for testname in self.testnames:
            self.figs['MB_%s' % testname] = os.path.join(self.figspath, \
                'MASTER_BIAS_%s.png' % testname)

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

        for testname in self.testnames:
            self.outcdps['MB_%s' % testname] = os.path.join(
                self.cdpspath, 'MASTER_BIAS_%s.fits' % testname)

    def _get_MBdict(self, testname):
        """ """

        PT = self.ParsedTable[testname]

        MBdict = dict()

        for jY in range(self.NCOLS_FPA):
            for iX in range(self.NSLICES_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                locator = self.fpa.FPA_MAP[Ckey]

                block = locator[0]
                CCDk = locator[1]
                flip = locator[2]

                inventoryitem = self.inventory[block][testname][-1] # last test

                productspath = os.path.join(inventoryitem['resroot'], 'products')

                ixblock = np.where(PT['BLOCK'] == block)
                

                cdpkey = PT['MASTERBIAS_%s' % CCDk][ixblock][0]

                masterpick = os.path.join(productspath, self.products['MASTERBIAS'][cdpkey])
                masterfits = st.replace(masterpick,'.pick','.fits')

                ccdobj = ccdmod.CCD(infits=masterfits, getallextensions=True, withpover=True,
                                    overscan=20)

                img = ccdobj.extensions[1].data.transpose().copy()


                # coarse hot pixel masking (just for display)
                medval = np.median(img)
                img[np.where(img>100)] = medval

                # image smoothing

                simg = ndimage.filters.gaussian_filter(img, sigma=5.,
                                                       mode='constant',
                                                       cval=medval)

                esimg = exposure.equalize_hist(simg, nbins=256)

                MBdict[Ckey] = dict(img=self.fpa.flip_img(esimg, flip))

                #ccdobj = None
                #img = None
                #simg = None
                #esimg = None

        return MBdict


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

    def _do_OffRonvsTime(self, testname):
        """ """

        datadict = OrderedDict()
        cols = ['CCDQID', 'BLOCK', 'CCD', 'time', 'REP', 'OFF', 'RON']
        for col in cols:
            datadict[col] = []

        rawptdf = self.ParsedTable[testname].to_pandas()
            
        for block in self.blocks:
            
            ixsel = np.where(rawptdf.BLOCK == block)

            for iix in ixsel[0]:

                for CCD in self.CCDs:

                    for Q in self.Quads:

                        ccdqid = '%s_%s' % (rawptdf['sn_ccd%i_CCD%i' % (CCD,CCD)][iix],Q)
                        datadict['CCDQID'].append(ccdqid)
                        datadict['BLOCK'].append(rawptdf.BLOCK[iix])
                        datadict['CCD'].append(CCD)
                        datadict['time'].append(rawptdf['time_CCD%i' % CCD][iix])
                        datadict['REP'].append(rawptdf.REP[iix])                        
                        datadict['OFF'].append(rawptdf['OFF_OVE_CCD%i_Quad%s' % (CCD,Q)][iix])
                        datadict['RON'].append(rawptdf['RON_OVE_CCD%i_Quad%s' % (CCD,Q)][iix])
        for col in cols:
            datadict[col] = np.array(datadict[col])

        ptdf = pd.DataFrame(datadict)
        
        
        uCCDQIDs = np.unique(ptdf.CCDQID)

        # discriminating CCDs that were calibrated more than once and 
        # those that were not.

        Lreps = []
        Lnoreps = []
        for uCCDQID in uCCDQIDs:
            utime = ptdf.loc[ptdf.CCDQID==uCCDQID].time
            dtime = utime.max()-utime.min()
            if dtime.days > 10:
                Lreps.append(uCCDQID)
            else:
                Lnoreps.append(uCCDQID)
    
        noreps = pd.DataFrame(dict(CCDQID=Lnoreps))
        reps = pd.DataFrame(dict(CCDQID=Lreps))

        ptnorep_df = ptdf[ptdf.CCDQID.isin(noreps.CCDQID)]
        ptrep_df = ptdf[ptdf.CCDQID.isin(reps.CCDQID)]

        def get_avtempmagstat(df, col, ccdqids, stat, debug=False):
            fstats = dict(mean=np.nanmean,
                            std=np.nanstd)
            statsvec = np.array([fstats[stat](df.loc[df.CCDQID==iid][col]) for iid in ccdqids])
            if debug: stop()
            return np.mean(statsvec)

        stdOFFnorep = get_avtempmagstat(ptnorep_df, 'OFF', Lnoreps, 'std')
        stdOFFrep = get_avtempmagstat(ptrep_df, 'OFF', Lreps, 'std')

        meanRONnorep = get_avtempmagstat(ptnorep_df, 'RON', Lnoreps, 'mean') 
        stdRONnorep = get_avtempmagstat(ptnorep_df, 'RON', Lnoreps, 'std')
        meanRONrep = get_avtempmagstat(ptrep_df, 'RON', Lreps, 'mean') 
        stdRONrep = get_avtempmagstat(ptrep_df, 'RON', Lreps, 'std')

        res = dict(NOREPCAL=dict(
                        stdOFF=stdOFFnorep,
                        meanRON=meanRONnorep,
                        stdRON=stdRONnorep,
                        N=len(ptnorep_df)),
                  REPCAL=dict(
                    stdOFF=stdOFFrep,
                    meanRON=meanRONrep,
                    stdRON=stdRONrep,
                    N=len(ptrep_df)))
        
        return res

    def dump_aggregated_results(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='dump', 
                Title='Aggregated Results', level=0)
            
            self.add_DataAlbaran2Report()

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))


        doAll = True
        doRONMaps = doAll
        doOffMaps = doAll
        doProfs = doAll
        doMaster = doAll
        doTemp = doAll

        if doRONMaps:

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
                    ColorbarText='ADU',
                    figname=self.figs['RON_ADU_%s' % testname]))



            # RON maps, ELECTRONs
            for testname in self.testnames:

                RONEMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname],
                                                  extractor=self._get_extractor_RON_fromPT(units='E',
                                                                                           reg='ove'))

                stestname = st.replace(testname, '_', '\_')
                figkey1 = 'RON_ELE_%s' % testname
                figname1 = self.figs[figkey1]

                self.plot_SimpleMAP(RONEMAP, **dict(
                    suptitle='%s: RON [ELECTRONS]' % stestname,
                    ColorbarText='electrons',
                    figname=figname1))


                if self.report is not None:
                    self.addFigure2Report(figname1, 
                            figkey=figkey1, 
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


        if doOffMaps:

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
                
                figkey2 = 'OFFSETS_%s' % testname
                figname2 = self.figs[figkey2]

                stestname = st.replace(testname, '_', '\_')
                self.plot_SimpleMAP(OFF_OVE_MAP, **dict(
                    suptitle='%s [OVERSCAN]: OFFSET' % stestname,
                    ColorbarText='ADU',
                    figname=figname2))

                if self.report is not None:
                    self.addFigure2Report(figname2, 
                            figkey=figkey2, 
                            caption='%s: Avg. offsets in ADU. Measured in serial over-scan.' % stestname, 
                            texfraction=0.7)

                if self.report is not None:

                    offcdpdict = dict(
                        caption='%s: Offset levels [ADU].' % testname,
                        valformat='%.2f')
                    
                    def _getOFF_val(Ckey, Q):
                        return OFFMAP[Ckey][Q]['ove']

                    _ = self.add_StdQuadsTable2Report( 
                                    extractor = _getOFF_val,
                                    cdpdict = offcdpdict)

        if doProfs:

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

                stestname = st.replace(testname, '_', '\_')

                for proftype in proftypes:

                    XY_profs = self._get_XYdict_PROFS(test=testname,proftype=proftype)

                    if proftype=='hor':
                        xlim=[0,100]
                        ylim=[-5,25]
                    else:
                        xlim=None
                        ylim=[-5,5]

                    figkey3 = 'PROFS_%s_%s' % (testname,proftype)
                    figname3 = self.figs[figkey3]

                    profkwargs = dict(
                        title='%s: Avg. profiles, direction: %s' % (testname, proftype),
                        doLegend=False,
                        xlabel=xlabels_profs[proftype],
                        ylabel=r'$\delta ADU$',
                        ylim=ylim,
                        xlim=xlim,
                        figname=figname3,
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
                        'interesting structure of this type of profile is. Each colour corresponds to a'+\
                        ' different block (each with 3x4 quadrants).'
                    if self.report is not None:

                        self.addFigure2Report(figname3, 
                            figkey=figkey3, 
                            caption= captemp % (stestname,),
                            texfraction=0.7)

        # MASTER BIAS

        if doMaster:

            for testname in self.testnames:

                figkey4 = 'MB_%s' % testname
                figname4 = self.figs[figkey4]

                MBdict = self._get_MBdict(testname)
                MBkwargs = dict(  # suptitle='%s, %s: Master Flat Field' % (stestname, colkey),
                    suptitle='%s: Master "Bias"' % testname, 
                    figname=figname4)

                self.plot_ImgFPA(MBdict, **MBkwargs)

                if self.report is not None:
                    self.addFigure2Report(figname4,
                        figkey=figkey4,
                        caption='%s: Master Bias [ADU].' % testname, 
                        texfraction=0.7)

        # Temporal Stability of Offsets


        if doTemp:

            res_time = self._do_OffRonvsTime('BIAS02')

            if self.report is not None:

                self.report.add_Section(keyword='timestab',
                    Title='Time Stability', level=1)

                self.report.add_Section(keyword='RONtimestab',
                    Title='RON Time Stability', level=2)                

                RonvsTimeText = [
                '\\textbf{Non Repeated} Calibrations, \n',
                'Nr. of measurements (quadrants x tests) = %i\n' % res_time['NOREPCAL']['N'],
                'mean RON: \\textbf{%.3f} ADU\n' % res_time['NOREPCAL']['meanRON'],
                'std RON: \\textbf{%.3f} ADU\n' % res_time['NOREPCAL']['stdRON'],
                '\\newline',
                '\\textbf{Repeated} Calibrations, \n',
                'Nr. of measurements (quadrants x tests) = %i\n' % res_time['REPCAL']['N'],
                'mean RON: \\textbf{%.3f} ADU\n' % res_time['REPCAL']['meanRON'],
                'std RON: \\textbf{%.3f} ADU\n' % res_time['REPCAL']['stdRON'],
                '\\newline'
                ]

                self.report.add_Text(RonvsTimeText)

                self.report.add_Section(keyword='OFFtimestab',
                    Title='Offset Time Stability', level=2)           

                
                OffvsTimeText = [
                '\\textbf{Non Repeated} Calibrations, \n',
                'Nr. of measurements (quadrants x tests) = %i\n' % res_time['NOREPCAL']['N'],
                'std OFF: \\textbf{%.3f} ADU\n' % res_time['NOREPCAL']['stdOFF'],
                '\\newline',
                '\\textbf{Repeated} Calibrations, \n',
                'Nr. of measurements (quadrants x tests) = %i\n' % res_time['REPCAL']['N'],
                'std OFF: \\textbf{%.3f} ADU\n' % res_time['REPCAL']['stdOFF']
                ]

                self.report.add_Text(OffvsTimeText)