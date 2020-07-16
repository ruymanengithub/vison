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
import pandas as pd

from vison.datamodel import cdp
from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal, utils
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files


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
    'offset_ove',
    'deltaoff_pre',
    'deltaoff_ove',
    'flu_med_img',
    'flu_std_img',
    'std_pre',
    'std_ove']


class MetaPTC(MetaCal):
    """ """

    def __init__(self, *args, **kwargs):
        """ """

        super(MetaPTC, self).__init__(*args, **kwargs)

        self.testnames = ['PTC01', 'PTC02_590', 'PTC02_730', 'PTC02_880']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        self.batches_highPRNU = ['14313', '14471']

        self.products['HER_CURVES'] = OrderedDict()

        self.censored = []

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
        #   BLOCK, TEST, REPEAT
        #   wavenm, calibrated HK (voltages),
        #   GAIN, EGAIN, ALPHA, BLOOM_ADU, BLOOM_E
        #   REFERENCES TO CURVES

        CHAMBER = sidd.meta['inputs']['CHAMBER']

        CHAMBER_key = CHAMBER[0]
        chamber_v = np.array([CHAMBER_key])
        sidd.addColumn(chamber_v, 'CHAMBERKEY', IndexS, ix=0)

        ogseobj = ogse.Ogse(CHAMBER=CHAMBER)

        wave = sidd.mx['wave'][0, 0]

        wave_v = np.array([ogseobj.get_wavelength(wave)])
        sidd.addColumn(wave_v, 'WAVENM', IndexS, ix=0)

        block_v = np.array([block])
        sidd.addColumn(block_v, 'BLOCK', IndexS, ix=0)

        test_v = np.array([jrep + 1])
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)

        test_v = np.array([session])
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)

        test_v = np.array([testname])
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)

        gain_mx = sidd.products['gain_mx']
        bloom_mx = sidd.products['bloom_mx']

        HER_dict = sidd.compliances['HER']

        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        gain_v = tmp_v_CQ.copy()
        egain_v = tmp_v_CQ.copy()
        alpha_v = tmp_v_CQ.copy()
        bloom_ADU_v = tmp_v_CQ.copy()
        bloom_e_v = tmp_v_CQ.copy()
        HER_v = tmp_v_CQ.copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                gain_v[0, iCCD, kQ] = gain_mx[CCDk][Q]['gain']
                egain_v[0, iCCD, kQ] = gain_mx[CCDk][Q]['egain']
                alpha_v[0, iCCD, kQ] = gain_mx[CCDk][Q]['alpha']

                bloom_ADU_v[0, iCCD, kQ] = bloom_mx[CCDk][Q]['bloom_ADU']
                bloom_e_v[0, iCCD, kQ] = bloom_mx[CCDk][Q]['bloom_e']

                HER_v[0, iCCD, kQ] = HER_dict[CCDk][Q][1]

        sidd.addColumn(gain_v, 'GAIN', IndexCQ)
        sidd.addColumn(egain_v, 'EGAIN', IndexCQ)
        sidd.addColumn(alpha_v, 'ALPHA', IndexCQ)
        sidd.addColumn(bloom_ADU_v, 'BLOOM_ADU', IndexCQ)
        sidd.addColumn(bloom_e_v, 'BLOOM_E', IndexCQ)
        sidd.addColumn(HER_v, 'HER', IndexCQ)

        productspath = os.path.join(inventoryitem['resroot'], 'products')
        her_pick = os.path.join(productspath, os.path.split(sidd.products['HER_PROFILES'])[-1])
        her_profs = files.cPickleRead(her_pick)['data'].copy()

        herprofkeys_v = np.zeros((1,NCCDs), dtype='S50')

        for iCCD, CCDk in enumerate(CCDkeys):

            herkey = '%s_%s_%s_%i_%s' % (testname, block, session, jrep + 1, CCDk)

            self.products['HER_CURVES'][herkey] = her_profs.copy()

            herprofkeys_v[0, iCCD] = herkey

        sidd.addColumn(herprofkeys_v, 'HERPROF_KEY', IndexC)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _extract_GAIN_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)
        column = 'GAIN_%s_Quad%s' % (CCDk, Q)
        G = PT[column][ixblock][0]
        return G

    def _extract_HER_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)
        column = 'HER_%s_Quad%s' % (CCDk, Q)
        HER = PT[column][ixblock][0]
        return HER

    def _extract_HERprof_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)
        column = 'HERPROF_KEY_%s' % (CCDk, )
        HERprof_key = PT[column][ixblock][0]
        data = self.products['HER_CURVES'][HERprof_key][CCDk][Q]
        ixsel = np.where(data['x']>51+2048)
        HERy = data['y'][ixsel].tolist()

        
        return HERy


    def _extract_BADU_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)

        column = 'BLOOM_ADU_%s_Quad%s' % (CCDk, Q)
        badu = PT[column][ixblock][0]
        if badu > 0:
            return badu
        else:
            return np.nan

    def _extract_BE_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)

        column = 'BLOOM_E_%s_Quad%s' % (CCDk, Q)
        be = PT[column][ixblock][0]
        if be > 0:
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
                        G = PT['GAIN_%s_Quad%s' % (CCDk, Q)].data[ixblock][0]
                        EG = PT['EGAIN_%s_Quad%s' % (CCDk, Q)].data[ixblock][0]
                        gpair = (G, EG)
                        G_MX[block][testname][CCDk][Q] = gpair

        return G_MX

    def _do_GvsTime(self):
        """ """

        def reshape_PT(PT, blocks, CCDs, Quads):

            datadict = OrderedDict()
            cols = ['CCDQID', 'BLOCK', 'CCD', 'time', 'TEST', 'WAVENM', 'GAIN', 'EGAIN']
            for col in cols:
                datadict[col] = []

            ptdf = PT.to_pandas()
            
            for block in blocks:
                ixsel = np.where(ptdf.BLOCK == block)

                for iix in ixsel[0]:

                    for CCD in CCDs:

                        for Q in Quads:

                            ccdqid = '%s_%s' % (ptdf['sn_ccd%i_CCD%i' % (CCD,CCD)][iix],Q)
                            datadict['CCDQID'].append(ccdqid)
                            datadict['BLOCK'].append(ptdf.BLOCK[iix])
                            datadict['CCD'].append(CCD)
                            datadict['time'].append(ptdf['time_CCD%i' % CCD][iix])
                            datadict['TEST'].append(ptdf.TEST[iix])
                            datadict['WAVENM'].append(ptdf.WAVENM[iix])
                            datadict['GAIN'].append(ptdf['GAIN_CCD%i_Quad%s' % (CCD,Q)][iix])
                            datadict['EGAIN'].append(ptdf['EGAIN_CCD%i_Quad%s' % (CCD,Q)][iix])
            for col in cols:
                datadict[col] = np.array(datadict[col])

            optdf = pd.DataFrame(datadict)

            return optdf


        for itest, testname in enumerate(self.testnames):

            t_PT_df = reshape_PT(self.ParsedTable[testname],
                                self.blocks,
                                self.CCDs, self.Quads)
            
            if itest==0:
                PT_df = copy.deepcopy(t_PT_df)
            else:
                PT_df = pd.concat([PT_df, t_PT_df])

        uCCDQIDs = np.unique(PT_df.CCDQID)

        # discriminating CCDs that were calibrated more than once and 
        # those that were not.

        Lreps = []
        Lnoreps = []
        for uCCDQID in uCCDQIDs:
            utime = PT_df.loc[PT_df.CCDQID==uCCDQID].time
            dtime = utime.max()-utime.min()
            if dtime.days > 10:
                Lreps.append(uCCDQID)
            else:
                Lnoreps.append(uCCDQID)
        
        noreps = pd.DataFrame(dict(CCDQID=Lnoreps))
        reps = pd.DataFrame(dict(CCDQID=Lreps))

        PTnorep_df = PT_df[PT_df.CCDQID.isin(noreps.CCDQID)]
        PTrep_df = PT_df[PT_df.CCDQID.isin(reps.CCDQID)]

        # average gain values for each wavelength

        res = PTnorep_df.groupby(['WAVENM']).median()['GAIN']
        avG = res.values.copy()
        waves = res.index.values.copy()
        fcorr = avG/avG[np.where(waves==800)]

        # correcting gains for wavelength effect
        
        for iw,wave in enumerate(waves):

            PTnorep_df.loc[PTnorep_df.WAVENM==wave,
                            PTnorep_df.columns=='GAIN'] /= fcorr[iw]
            PTrep_df.loc[PTrep_df.WAVENM==wave,
                            PTrep_df.columns=='GAIN'] /= fcorr[iw]

        
        def get_avtempGrelstd(df,ccdqids):
            gains = np.array([df.loc[df.CCDQID==iid].GAIN.mean() for iid in ccdqids])
            stdgains = np.array([df.loc[df.CCDQID==iid].GAIN.std() for iid in ccdqids])
            return np.median(stdgains/gains)

        stdGnorep = get_avtempGrelstd(PTnorep_df, Lnoreps)
        stdGrep = get_avtempGrelstd(PTrep_df, Lreps)

        res = dict(NOREPCAL=dict(relstdG=stdGnorep,
                    N=len(Lnoreps)),
                  REPCAL=dict(
                    relstdG=stdGrep,
                    N=len(Lreps)),
                  Gwave=dict(
                    waves=waves,
                    avgG=avG))

        return res




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

            self.figs['BLOOM_ADU_MAP_%s' % testname] = os.path.join(
                self.figspath, 'BLOOM_ADU_MAP_%s.png' % testname)

            self.figs['BLOOM_ELE_MAP_%s' % testname] = os.path.join(
                self.figspath, 'BLOOM_ELE_MAP_%s.png' % testname)

            self.figs['HER_MAP_%s' % testname] = os.path.join(self.figspath,
                                                              'HER_MAP_%s.png' % testname)

            self.figs['HER_curves_%s' % testname] = os.path.join(self.figspath,
                                                                 'HER_curves_%s.png' % testname)

    def init_outcdpnames(self):

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        for testname in self.testnames:
            self.outcdps['GAIN_%s' % testname] = '%s_GAIN_elperadu_MAP.json' % testname
            self.outcdps['HER_%s' % testname] = '%s_HER_profiles_MAP.json' % testname
            self.outcdps['BLOOM_ADU_%s' % testname] = '%s_BLOOM_ADU_MAP.json' % testname
            self.outcdps['BLOOM_ELE_%s' % testname] = '%s_BLOOM_ELE_MAP.json' % testname

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

                        Gcol = 'GAIN_CCD%i_Quad%s' % (CCD, Q)

                        _y.append(PT[Gcol][ixblock][0])

            y.append(np.nanmean(_y))
            ey.append(np.nanstd(_y))

        sort = np.argsort(x)
        x = np.array(x)[sort]
        y = np.array(y)[sort]
        ey = np.array(ey)[sort]

        XYdict = dict(x=x, y=y, ey=ey)

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

            for iCCD, CCD in enumerate(self.CCDs):
                TCCD = (PT['HK_CCD%s_TEMP_T' % CCD] + PT['HK_CCD%s_TEMP_B' % CCD]) / 2.
                for iQ, Q in enumerate(self.Quads):

                    Ycolumn = 'GAIN_CCD%s_Quad%s' % (CCD, Q)

                    _Ydata = PT[Ycolumn].copy()

                    x[testname] += TCCD.tolist()
                    y[testname] += _Ydata.tolist()

            x[testname] = np.array(x[testname])
            y[testname] = np.array(y[testname])

        XYdict = dict(x=x, y=y, labelkeys=labelkeys)

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

            for iCCD, CCD in enumerate(self.CCDs):

                for iQ, Q in enumerate(self.Quads):

                    if Q in ['E', 'F']:
                        Xcolumn = 'HK_CCD%i_OD_T_CAL' % (CCD,)
                    else:
                        Xcolumn = 'HK_CCD%i_OD_B_CAL' % (CCD,)

                    _Xdata = PT[Xcolumn].copy()

                    Ycolumn = 'GAIN_CCD%s_Quad%s' % (CCD, Q)
                    _Ydata = PT[Ycolumn].copy()

                    x[testname] += _Xdata.tolist()
                    y[testname] += _Ydata.tolist()

            x[testname] = np.array(x[testname])
            y[testname] = np.array(y[testname])

        XYdict = dict(x=x, y=y, labelkeys=labelkeys)

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

            for iCCD, CCD in enumerate(self.CCDs):

                for iQ, Q in enumerate(self.Quads):

                    if Q in ['E', 'F']:
                        Xcolumn = 'HK_COMM_RD_T_CAL'
                    else:
                        Xcolumn = 'HK_COMM_RD_B_CAL'

                    _Xdata = PT[Xcolumn].copy()

                    Ycolumn = 'GAIN_CCD%s_Quad%s' % (CCD, Q)
                    _Ydata = PT[Ycolumn].copy()

                    x[testname] += _Xdata.tolist()
                    y[testname] += _Ydata.tolist()

            x[testname] = np.array(x[testname])
            y[testname] = np.array(y[testname])

        XYdict = dict(x=x, y=y, labelkeys=labelkeys)

        return XYdict

    def _get_XYdict_HER(self, testname):
        """ """

        x = dict()
        y = dict()

        PT = self.ParsedTable[testname]

        labelkeys = []

        for block in self.flight_blocks:
            ixsel = np.where(PT['BLOCK'] == block)

            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' % CCD

                herprof_key = PT['HERPROF_KEY_%s' % CCDk][ixsel][0]
                i_her = self.products['HER_CURVES'][herprof_key].copy()

                for kQ, Q in enumerate(self.Quads):

                    pkey = '%s_%s_%s' % (block, CCDk, Q)

                    _x = i_her[CCDk][Q]['x'][10:].copy()
                    _x -= _x.min()
                    _y = i_her[CCDk][Q]['y'][10:].copy()

                    if pkey not in self.censored:

                        x[pkey] = _x.copy()
                        y[pkey] = _y.copy()
                        labelkeys.append(pkey)


        HERdict = dict(x=x, y=y, labelkeys=labelkeys)

        return HERdict


    def dump_aggregated_results(self):
        """ """

        if self.report is not None:
            self.report.add_Section(keyword='dump', 
                Title='Aggregated Results', level=0)

            self.add_DataAlbaran2Report()


        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        outpathroot = self.outpathroot

        doAll = True

        doGainMaps = doAll
        doBloomMaps = doAll
        doHERMaps = doAll
        doHERcurves = doAll
        doGvsWave = doAll
        doGvsT = doAll
        doGvsOD = doAll
        doGvsRD = doAll
        doGvsTime = doAll

        # GAIN matrix (all blocks and test/waves) to dict() saved as pickle

        GAIN_MXdict = self.gen_GAIN_MXdict()

        GAIN_MXpick = os.path.join(outpathroot, 'GAIN_MX_PTC0X.pick')

        files.cPickleDump(GAIN_MXdict, GAIN_MXpick)

        # GAIN maps (all tests/waves)

        if doGainMaps:

            for testname in self.testnames:

                GMAP = self.get_FPAMAP_from_PT(
                    self.ParsedTable[testname],
                    extractor=self._extract_GAIN_fromPT)

                stestname = st.replace(testname, '_', '\_')

                avgG = self.get_stat_from_FPAMAP(GMAP, np.nanmean)
<<<<<<< HEAD
                print(('Average G [%s]: %.2f' % (testname, avgG)))
=======
                avgGtext = 'Average G [%s]: %.2f' % (stestname, avgG)
                if self.report is not None:
                    self.report.add_Text(avgGtext)
                else:
                    print(avgGtext)
>>>>>>> master

                def _count_within_REQ(vals):
                    req = [3.0,3.5]
                    return np.sum(np.array([(item>req[0] and item<=req[1]) for item in vals]))

                NwithinREQ = self.get_stat_from_FPAMAP(GMAP, _count_within_REQ)
                ReqText = 'N-Quadrants within requirement (3.0-3.5 e/ADU): %i' % (NwithinREQ)
                if self.report is not None:
                    self.report.add_Text(ReqText)
                else:
                    print(ReqText)


                figkey1 = 'GAIN_MAP_%s' % testname
                figname1 = self.figs[figkey1]

                self.plot_SimpleMAP(GMAP, **dict(
                    suptitle='%s: GAIN e-/ADU' % stestname,
                    ColorbarText = 'e-/ADU',
                    figname=figname1))

                if self.report is not None:
                    self.addFigure2Report(figname1, 
                        figkey=figkey1, 
                        caption='%s: Gain Map [e-/ADU].' % stestname, 
                        texfraction=0.7)

                g_header = OrderedDict()
                g_header['title'] = 'GAIN MAP'
                g_header['test'] = stestname
                g_header.update(CDP_header)
            
                g_cdp = cdp.Json_CDP(rootname=self.outcdps['GAIN_%s' % testname],
                              path=self.cdpspath)
                g_cdp.ingest_inputs(data=GMAP,
                             header = g_header,
                             meta=dict(units='e-/ADU',
                                        structure='CCDID:Q:gain'))
                g_cdp.savehardcopy()


        if doBloomMaps:

            # BLOOM maps (ADU and e-, from PTC01)

            for testname in self.testnames:

                BADU_MAP = self.get_FPAMAP_from_PT(
                    self.ParsedTable[testname], extractor=self._extract_BADU_fromPT)

                stestname = st.replace(testname, '_', '\_')

                figkey2 = 'BLOOM_ADU_MAP_%s' % testname
                figname2 = self.figs[figkey2]

                self.plot_SimpleMAP(BADU_MAP, **dict(
                    suptitle='%s: BLOOM-ADU [DN]' % stestname,
                    ColorbarText='ADU',
                    figname=figname2))  # ,
                # corekwargs=dict(norm = Normalize(vmin=3e4,vmax=2**16, clip=False))))

                if self.report is not None:
                    self.addFigure2Report(figname2, 
                        figkey=figkey2, 
                        caption='%s: Blooming onset threshold map in ADU.' % stestname, 
                        texfraction=0.7)


                ba_header = OrderedDict()
                ba_header['title'] = 'BLOOMING MAP [ADU]'
                ba_header['test'] = stestname
                ba_header.update(CDP_header)
            
                ba_cdp = cdp.Json_CDP(rootname=self.outcdps['BLOOM_ADU_%s' % testname],
                              path=self.cdpspath)
                ba_cdp.ingest_inputs(data=BADU_MAP,
                             header = ba_header,
                             meta=dict(units='ADU',
                                        structure='CCDID:Q:bloom threshold'))
                ba_cdp.savehardcopy()

            for testname in self.testnames:

                BE_MAP = self.get_FPAMAP_from_PT(
                    self.ParsedTable[testname],
                    extractor=self._extract_BE_fromPT)

                figkey3 = 'BLOOM_ELE_MAP_%s' % testname
                figname3 = self.figs[figkey3]

                stestname = st.replace(testname, '_', '\_')
                self.plot_SimpleMAP(BE_MAP, **dict(
                    suptitle='%s: BLOOM-ELECTRONS' % stestname,
                    ColorbarText='electrons',
                    figname=figname3))  # ,
                # corekwargs=dict(norm = Normalize(vmin=1e5,vmax=2.2E5, clip=False))))

                if self.report is not None:
                    self.addFigure2Report(figname3, 
                        figkey=figkey3, 
                        caption='%s: Blooming Map in electrons.' % stestname, 
                        texfraction=0.7)

                be_header = OrderedDict()
                be_header['title'] = 'BLOOMING MAP [ELECTRONS]'
                be_header['test'] = stestname
                be_header.update(CDP_header)
            
                be_cdp = cdp.Json_CDP(rootname=self.outcdps['BLOOM_ELE_%s' % testname],
                              path=self.cdpspath)
                be_cdp.ingest_inputs(data=BE_MAP,
                             header = be_header,
                             meta=dict(units='electrons',
                                        structure='CCDID:Q:bloom threshold'))
                be_cdp.savehardcopy()


        # HER map

        if doHERMaps:

            for testname in self.testnames:

                HERMAP = self.get_FPAMAP_from_PT(
                    self.ParsedTable[testname],
                    extractor=self._extract_HER_fromPT)

                stestname = st.replace(testname, '_', '\_')

                figkey4 = 'HER_MAP_%s' % testname
                figname4 = self.figs[figkey4]

                self.plot_SimpleMAP(HERMAP, **dict(
                    suptitle='%s: Hard Edge Response Factor' % stestname,
                    ColorbarText='[adim]',
                    figname=figname4))

                if self.report is not None:
                    captemp = '%s: Hard Edge Response Map,'+\
                    ' first pixel relative response.'

                    self.addFigure2Report(figname4, 
                        figkey=figkey4, 
                        caption= captemp % stestname, 
                        texfraction=0.7)

                HERMAPprofs = self.get_FPAMAP_from_PT(
                    self.ParsedTable[testname],
                    extractor=self._extract_HERprof_fromPT)

                
                her_header = OrderedDict()
                her_header['title'] = 'HARD EDGE RESPONSE'
                her_header['test'] = stestname
                her_header.update(CDP_header)
            
                her_cdp = cdp.Json_CDP(rootname=self.outcdps['HER_%s' % testname],
                              path=self.cdpspath)
                her_cdp.ingest_inputs(data=HERMAPprofs,
                             header = her_header,
                             meta=dict(units='adim',
                                structure="CCDID:Q:[serial profile]"))
                her_cdp.savehardcopy()
                

        # HER Curves

        if doHERcurves:

            for testname in self.testnames:

                HERSingledict = self._get_XYdict_HER(testname)

                stestname = st.replace(testname, '_', '\_')

                figkey5 = 'HER_curves_%s' % testname
                figname5 = self.figs[figkey5]

                self.plot_XY(HERSingledict, **dict(
                    title='%s: H.E.R. CURVES' % stestname,
                    doLegend=False,
                    xlabel='Pixel',
                    ylabel='HER [frac]',
                    ylim=[-2.E-4, 5.e-4],
                    xlim=[0, 6],
                    corekwargs=dict(linestyle='-', marker=''),
                    figname=figname5))

                if self.report is not None:

                    captemp = '%s: H.E.R. curves, all 144 quadrant channels.'

                    self.addFigure2Report(figname5, 
                        figkey=figkey5, 
                        caption=captemp % stestname, 
                        texfraction=0.7)

        # GAIN vs. Wavelength

        if doGvsWave:

            GvsWavedict = self._get_XYdict_GvsLAM()

            figkey6 = 'GvsWave'
            figname6 = self.figs[figkey6]

            self.plot_XY(GvsWavedict, **dict(
                title='Gain vs. Wavelength',
                doLegend=False,
                doYErrbars=True,
                xlabel='Wavelength',
                ylabel='Gain [e-/ADU]',
                ylim=[3.3, 3.7],
                corekwargs=dict(linestyle='', marker='o'),
                # figname = ''))
                figname=figname6))

            if self.report is not None:
                self.addFigure2Report(figname6, 
                        figkey=figkey6, 
                        caption='Measured Gain [e-/ADU] as a function of wavelength.', 
                        texfraction=0.7)

        if doGvsTime:
            
            res_GvsTime = self._do_GvsTime()

            if self.report is not None:

                avgGw = ['%.3f' % item for item in res_GvsTime['Gwave']['avgG']]

                GvsTimeText = ['Wavelengths = %s nm\n' % res_GvsTime['Gwave']['waves'].tolist().__repr__(),
                'Avg. Gain in each Wavelength = %s\n' % avgGw.__repr__(),
                '\\textbf{WARNING}: Wavelength dependency of G (through B-F) has been removed in'+\
                ' the following stats.\n',
                '\\textbf{Non Repeated} Calibrations,\n '+\
                    'Nr. of measurements (CCDs x tests) = %i\n' % res_GvsTime['NOREPCAL']['N'],
                'mean rel. std. of G values per quadrant across tests: %.2f \\%%.\n' % \
                            (res_GvsTime['NOREPCAL']['relstdG']*100.,),
                '\\textbf{Repeated} Calibrations,\n '+\
                    'Nr. of measurements (CCDs x tests) = %i\n' % res_GvsTime['REPCAL']['N'],
                'mean rel. std. of G values per quadrant across Tests: %.2f \\%%.\n' % \
                            (res_GvsTime['REPCAL']['relstdG']*100.,),
                ]

                self.report.add_Text(GvsTimeText)

                

        # GAIN vs. detector temperature (PTC01)

        if doGvsT:

            GvsTdict = self._get_XYdict_GvsT()

            figkey7 = 'GvsT'
            figname7 = self.figs[figkey7]

            self.plot_XY(GvsTdict, **dict(
                title='Gain vs. Detector Temperature',
                doLegend=True,
                xlabel='Detector Temperature',
                ylabel='Gain [e-/ADU]',
                corekwargs=dict(linestyle='', marker='.'),
                figname=figname7))

            if self.report is not None:

                self.addFigure2Report(figname7, 
                        figkey=figkey7, 
                        caption='Measured Gain (at 800 nm) vs. detector Temperature,'+\
                                ' all quadrant channels.', 
                        texfraction=0.7)

        # GAIN vs. OD-CAL (PTC01)

        if doGvsOD:

            GvsODdict = self._get_XYdict_GvsOD()

            figkey8 = 'GvsOD'
            figname8 = self.figs[figkey8]

            self.plot_XY(GvsODdict, **dict(
                title='Gain vs. Output Drain',
                doLegend=True,
                xlabel='Output Drain (Calibrated) [V]',
                ylabel='Gain [e-/ADU]',
                corekwargs=dict(linestyle='', marker='.'),
                # figname=''))
                figname=figname8))

            if self.report is not None:
                self.addFigure2Report(figname8, 
                        figkey=figkey8, 
                        caption='Gain vs. Output Drain Voltage,'+\
                        ' all 144 quadrants. The (small) range in OD is given '+\
                        'from differences in voltage outputs / calibration'+\
                        ' across ROEs.', 
                        texfraction=0.7)

        # GAIN vs. RD-CAL (PTC01)

        if doGvsRD:

            figkey9 = 'GvsRD'
            figname9 = self.figs[figkey9]

            GvsRDdict = self._get_XYdict_GvsRD()

            self.plot_XY(GvsRDdict, **dict(
                title='Gain vs. Reset Drain',
                doLegend=True,
                xlabel='Reset Drain (Calibrated) [V]',
                ylabel='Gain [e-/ADU]',
                corekwargs=dict(linestyle='', marker='.'),
                # figname=''))
                figname=figname9))

            if self.report is not None:
                self.addFigure2Report(figname9, 
                        figkey=figkey9, 
                        caption='Gain vs. Reset Drain Voltage,'+\
                        ' all 144 quadrants. The (small) range in RD is given '+\
                        'from differences in voltage outputs / calibration'+\
                        ' across ROEs.', 
                        texfraction=0.7)

        # Save the ParsedTable(s)
