#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:32:00 2020

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
from collections import OrderedDict
import string as st
import os
#from scipy import ndimage
#from skimage import exposure
import pandas as pd
from astropy.table import vstack

from vison.datamodel import cdp
from vison.datamodel import ccd as ccdmod
from vison.support import utils
from vison.support import files
from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.metatests.bias import MetaBias
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
    'offset_ove',
    'std_pre',
    'std_ove']

testnames = ['BIAS01',
'BIAS02',
'DARK01',
#'COSMETICS00',
'CHINJ01',
'CHINJ02',
'FLAT01',
'FLAT02_590',
'FLAT02_730',
'FLAT02_880',
'NL02',
'PERSIST01',
'PSF01_590',
'PSF01_730',
'PSF01_800',
'PSF01_880',
'PTC01',
'PTC02_590',
'PTC02_730',
'PTC02_880',
'TP01',
'TP02',
'TP11',
'TP21']


class MetaPano(MetaCal):

    def __init__(self, **kwargs):
        """ """

        super(MetaPano, self).__init__(**kwargs)

        self.testnames = testnames
        self.testtypes = ['BIAS','DARK','CHINJ','FLAT','NL','PERSIST','PSF','PTC','TP']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        #self.products['MB_PROFILES'] = OrderedDict()
        #self.products['MASTERBIAS'] = OrderedDict()

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

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit


    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['doff_colored_bytest'] = os.path.join(self.figspath, 'DELTAOFFSETS_byTEST.png')
        self.figs['dron_colored_bytest'] = os.path.join(self.figspath, 'DELTARON_byTEST.png')

    def init_outcdpnames(self):
        """ """

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

    def stack_PTs(self):
        """ """

        for it, testname in enumerate(self.testnames):

            if it ==0:
                PT = copy.deepcopy(self.ParsedTable[testname])
            else:
                PT = vstack([PT,self.ParsedTable[testname]])
        

        return PT


    def get_pldata_vstime_bytest(self, statistic):
        """ """

        pldata = OrderedDict(labelkeys=[],x=dict(),y=dict())

        PT = self.stack_PTs()

        for block in self.flight_blocks:

            for jCCD, CCD in enumerate(self.CCDs):
                CCDkey = 'CCD%i' % CCD
                for kQ, Q in enumerate(self.Quads):

                    statcol = '%s_%s_Quad%s' % (statistic, CCDkey, Q)
                    tcol = 'time_%s' % CCDkey

                    ixsel = np.where((PT['BLOCK']==block))

                    _data = PT[statcol][ixsel]
                    _data -= np.nanmean(_data)

                    _time = PT[tcol][ixsel]

                    deltatime = np.array([dt.seconds/3600. for dt in _time-_time.min()])

                    tests = PT['TEST'][ixsel]

                    for testname in self.testnames:

                        ixtest = np.where(tests == testname)

                        testtype = [ttype for ttype in self.testtypes if ttype in 'BIAS01'][0]

                        if testtype not in pldata['labelkeys']:
                            pldata['labelkeys'].append(testtype)

                            pldata['x'][testtype] = list(deltatime[ixtest])
                            pldata['y'][testtype] = list(_data[ixtest])
                        else:
                            pldata['x'][testtype] += list(deltatime[ixtest])
                            pldata['y'][testtype] += list(_data[ixtest])

        return pldata



    def dump_aggregated_results(self):
        """ """


        if self.report is not None:
            self.report.add_Section(keyword='dump', 
                Title='Aggregated Results', level=0)

            #self.add_DataAlbaran2Report()

            # Note:
            # Delta-offsets/RONs are referred to the mean value of offset/ron for each quadrant,
            # in each campaign, for tests BIAS02


            # Delta-offset vs. time, color coding by test type

            pldata1= self.get_pldata_vstime_bytest('offset_pre')



            figkey1 = 'doff_colored_bytest'
            figname1 = self.figs[figkey1]

            fig1kwargs = dict(
                title='Delta-Offset vs. time',
                doLegend=True,
                xlabel=r'$\Delta Time\ [hrs]$',
                ylabel=r'$\Delta Offset [ADU]$',
                #ylim=[-3., 7.],
                figname=figname1)


            TESTcolors = cm.rainbow(np.linspace(0, 1, len(self.testtypes)))

            pointcorekwargs1 = dict()
            for jtest, testtype in enumerate(self.testtypes):
                pointcorekwargs1['%s' % (testtype,)] = dict(
                    linestyle='', marker='.', color=TESTcolors[jtest], alpha=0.5)

            fig1kwargs['corekwargs'] = pointcorekwargs1

            self.plot_XY(pldata1, **fig1kwargs)

            if self.report is not None:
                self.addFigure2Report(figname1,
                    figkey=figkey1,
                    caption='',
                    texfraction=0.8)


            # Delta-RON vs. time, color coding by test type

            pldata2= self.get_pldata_vstime_bytest('std_pre')

            figkey2 = 'dron_colored_bytest'
            figname2 = self.figs[figkey2]

            fig2kwargs = dict(
                title='Delta-RON vs. time',
                doLegend=True,
                xlabel=r'$\Delta Time\ [hrs]$',
                ylabel=r'$\Delta RON [ADU]$',
                #ylim=[-3., 7.],
                figname=figname2)

            TESTcolors = cm.rainbow(np.linspace(0, 1, len(self.testtypes)))

            pointcorekwargs2 = dict()
            for jtest, testtype in enumerate(self.testtypes):
                pointcorekwargs2['%s' % (testtype,)] = dict(
                    linestyle='', marker='.', color=TESTcolors[jtest], alpha=0.5)

            fig2kwargs['corekwargs'] = pointcorekwargs2

            self.plot_XY(pldata2, **fig2kwargs)

            if self.report is not None:
                self.addFigure2Report(figname2,
                    figkey=figkey2,
                    caption='',
                    texfraction=0.8)



            # Delta-offset vs. time, color coding by block

            

            # histogram of delta-offsets (prescans)


            # histogram of delta-RONs (prescans)


            


            # pre-scan offset vs. image fluence (requires changing all tests to measure image fluences
            # in check stats!)

            




        


