#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Created on Mon Tue 24 16:34:00 2020

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

from vison.support import vcal, utils
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files
from vison.datamodel import ccd as ccdmod
from vison.datamodel import cdp as cdpmod
from vison.plot.baseplotclasses import XYPlot

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
    'HK_RPSU_28V_PRI_I']

class XYPlotPER(XYPlot):

    def populate_axes(self):
        """ """
        super(XYPlotPER,self).populate_axes()        
        self.ax.axvline(x=0., linestyle='--',color='r')
        self.ax.axhline(y=0, linestyle=':', color='b')




class MetaPersist(MetaCal):
    """ """

    def __init__(self, *args, **kwargs):
        """ """

        super(MetaPersist, self).__init__(*args, **kwargs)

        self.testnames = ['PERSIST01']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        # self.blocks = self.blocks[-2:] # TESTS

        self.products['PERSIST'] = OrderedDict()

        self.init_fignames()
        #self.init_outcdpnames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """


        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']

        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        # IndexC = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
        #                             vcore.vIndex('CCD',vals=self.CCDs)])

        IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0]),
                                     vcore.vIndex('CCD', vals=self.CCDs),
                                     vcore.vIndex('Quad', vals=self.Quads)])

        #idd = copy.deepcopy(inventoryitem['dd'])

        sidd = self.parse_single_test_gen(jrep, block, testname, inventoryitem)

        # TEST SPECIFIC


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


        productspath = os.path.join(inventoryitem['resroot'], 'products')

        persist_pick = os.path.join(productspath,
                os.path.split(sidd.products['PERSIST_STATS'])[-1])
        persist = files.cPickleRead(persist_pick)
        #stats = copy.deepcopy(metafitcdp['data']['ANALYSIS'])

        persistkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['PERSIST'][persistkey] = copy.deepcopy(persist)
        persistkey_v = np.array([persistkey])
        sidd.addColumn(persistkey_v, 'PERSIST', IndexS, ix=0)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

        # ADDING Nr. pixels masked: DARK, FLAT, MERGEs

        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        Ndark_v = tmp_v_CQ.copy()
        Nflat_v = tmp_v_CQ.copy()
        Nmerge_v = tmp_v_CQ.copy()
        Ncols_v = tmp_v_CQ.copy()

        productspath = os.path.join(inventoryitem['resroot'], 'products')

        NPIX_pick = os.path.join(productspath, os.path.split(sidd.products['DEF_TB_CDP'])[-1])

        NPIX_dict = files.cPickleRead(NPIX_pick)['data']['DEFECTS'].copy()

        for iCCD, CCDk in enumerate(CCDkeys):
            for kQ, Q in enumerate(self.Quads):

                ixsel = np.where((NPIX_dict['CCD'] == (iCCD + 1)) & (NPIX_dict['Q'] == (kQ + 1)))

                Ndark_v[0, iCCD, kQ] = NPIX_dict['N_DARK'].as_matrix()[ixsel][0]
                Nflat_v[0, iCCD, kQ] = NPIX_dict['N_FLAT'].as_matrix()[ixsel][0]
                Nmerge_v[0, iCCD, kQ] = NPIX_dict['N_MERGE'].as_matrix()[ixsel][0]
                Ncols_v[0, iCCD, kQ] = NPIX_dict['NCOLS_MERGE'].as_matrix()[ixsel][0]

        sidd.addColumn(Ndark_v, 'NDARK', IndexCQ)
        sidd.addColumn(Nflat_v, 'NFLAT', IndexCQ)
        sidd.addColumn(Nmerge_v, 'NMERGE', IndexCQ)
        sidd.addColumn(Ncols_v, 'NCOLS', IndexCQ)

        # Maps of bad pixels

        tmp_v_str = np.zeros((1), dtype='S50')

        for maskkey in self.maskkeys:

            mskkey_v = tmp_v_str.copy()

            all_mask_fits = dict()

            for jCCD, CCDk in enumerate(CCDkeys):

                all_mask_fits[CCDk] = os.path.join(
                    productspath, os.path.split(sidd.products[maskkey][CCDk])[-1])

            _mskkey = '%s_%s_%s_R%i' % (maskkey, block, session, jrep + 1)

            #debug=False
            #if block == 'GUYE' and maskkey=='MERGE':
            #    debug=True

            self.products['MASKSCOOS'][_mskkey] = self._extract_badpix_coordinates(
                all_mask_fits, block) #,debug=debug)

            mskkey_v[0] = _mskkey

            sidd.addColumn(mskkey_v, 'MASKCOOS_%s' % maskkey, IndexS)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit


    def plot_XY(self, XYdict, **kwargs):
        """ """
        plotobj = self._get_plot_gen(XYPlotPER)
        plotobj(self, XYdict, **kwargs)

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)


        self.figs['PERSIST_curves'] = os.path.join(self.figspath,
                        'PERSIST_curves.png')


    def _get_XYdict_PER(self):
        """ """

        x = dict()
        y = dict()
        ey = dict()

        PT = self.ParsedTable['PERSIST01']
        column = 'PERSIST'

        labelkeys = []

        pmeans = []
        pmeds = []
        p25s = []
        p75s = []
        pstds = []

        for iblock, block in enumerate(self.flight_blocks):

            ixblock = np.where(PT['BLOCK'] == block)
            per_key = PT[column][ixblock][0]
            per = self.products['PERSIST'][per_key]

            if iblock == 0:
                dtsec = np.concatenate((
                    per['data']['CCD1']['E']['REF']['deltasec'],
                    per['data']['CCD1']['E']['LAT']['deltasec']))

            for iCCD, CCD in enumerate(self.CCDs):

                CCDk = 'CCD%i' % CCD

                for kQ, Q in enumerate(self.Quads):

                    pmeans.append(np.concatenate(
                        (per['data'][CCDk][Q]['REF']['mean'],
                        per['data'][CCDk][Q]['LAT']['mean']))
                        )
                    pmeds.append(np.concatenate(
                        (per['data'][CCDk][Q]['REF']['p50'],
                        per['data'][CCDk][Q]['LAT']['p50']))
                        )

                    p25s.append(np.concatenate(
                        (per['data'][CCDk][Q]['REF']['p25'],
                        per['data'][CCDk][Q]['LAT']['p25']))
                        )

                    p75s.append(np.concatenate(
                        (per['data'][CCDk][Q]['REF']['p75'],
                        per['data'][CCDk][Q]['LAT']['p75']))
                        )

                    pstds.append(np.concatenate(
                        (per['data'][CCDk][Q]['REF']['std'],
                        per['data'][CCDk][Q]['LAT']['std']))
                        )

        Nprof = len(pmeans)

        avmeans = np.nanmedian(np.stack(pmeans), axis=0)
        avmeds = np.nanmedian(np.stack(pmeds), axis=0)
        avp25 = np.nanmedian(np.stack(p25s), axis=0)
        avp75 = np.nanmedian(np.stack(p75s), axis=0)
        avstd = np.nanmedian(np.stack(pstds), axis=0) / np.sqrt(Nprof)

        x['mean'] = dtsec.copy()
        y['mean'] = avmeans.copy()
        ey['mean'] = avstd.copy()
        labelkeys = ['mean']

        outdict = dict(x=x, y=y, ey=ey,
            labelkeys=labelkeys,
            confidence=dict(
                x = dtsec.copy(),
                yminus = avp25.copy(),
                yplus = avp75.copy()
                )
            )

        return outdict


    def dump_aggregated_results(self):
        """ """


        if self.report is not None:
            self.report.add_Section(keyword='dump', 
                Title='Aggregated Results', level=0)

            self.add_DataAlbaran2Report()

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header['DATE'] = self.get_time_tag()

        # Persistence Main Plot

        PER_singledict = self._get_XYdict_PER()

        figkey1 = 'PERSIST_curves'
        figname1 = self.figs[figkey1]

        PER_kwargs = dict(
            title='Average Persistence Curve',
            doLegend=False,
            doYErrbars=True,
            doConfidence=True,
            xlabel='time [seconds]',
            ylabel='Excess signal [ADUs]',
            figname=figname1,
            corekwargs=dict(
                mean=dict(linestyle='-',marker='o',color='b'))
            )

        PER_kwargs['confidence_kwargs'] = dict(
            color='g',
            alpha=0.2
            )

        self.plot_XY(PER_singledict, **PER_kwargs)

        if self.report is not None:
            self.addFigure2Report(figname1,
                figkey=figkey1,
                caption='PERSIST01: Persistence plot, stacked across all 144 '+\
                    'quadrants in the FPA. The shown signal is the average '+\
                    'within the saturation mask, minus the background, for '+\
                    'the 3 frames acquired before and after the saturation '+\
                    '(at t=0), stacked by the median across the 144 quadrants '+\
                    'in the FPA. The green shaded area marks the average extent '+\
                    'of the 25 to 75 percentile of pixel values within the mask. '+\
                    'The error bars show the (average) standard deviation '+\
                    'divided by sqrt(144), to give an idea of the statistical '+\
                    'error. No significant latent is measured in this test.',
                texfraction=1.)


