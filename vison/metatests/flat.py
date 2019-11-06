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
from skimage import exposure
from scipy import ndimage
import gc

from vison.support import vjson
from vison.datamodel import ccd as ccdmod
from vison.support import files
from vison.fpa import fpa as fpamod

from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse


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
    'flu_med_img',
    'flu_std_img',
    'std_pre',
    'std_ove']


class MetaFlat(MetaCal):
    """ """

    def __init__(self, **kwargs):
        """ """

        super(MetaFlat, self).__init__(**kwargs)

        self.testnames = ['FLAT01', 'FLAT02_590', 'FLAT02_730', 'FLAT02_880']

        self.colkeys = dict()
        for test in self.testnames:
            if test == 'FLAT01':
                self.colkeys[test] = ['col001', 'col002', 'col003']
            else:
                self.colkeys[test] = ['col001', 'col002']

        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['MASTERFLATS'] = OrderedDict()

        self.batches_highPRNU = ['14313', '14471']

        self.init_fignames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """

        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']
        Ncols = len(self.colkeys[testname])

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

        # ADDING PRNU AND FLUENCES

        tmp_v_CQ = np.zeros((1, NCCDs, NQuads))

        prnu_pc_vs = OrderedDict()
        prnu_fluadu_vs = OrderedDict()
        prnu_fluele_vs = OrderedDict()
        for icol in range(1, Ncols + 1):
            colkey = 'col%03i' % icol

            prnu_pc_vs[colkey] = tmp_v_CQ.copy()
            prnu_fluadu_vs[colkey] = tmp_v_CQ.copy()
            prnu_fluele_vs[colkey] = tmp_v_CQ.copy()

        productspath = os.path.join(inventoryitem['resroot'], 'products')

        prnucdp_pick = os.path.join(productspath, os.path.split(sidd.products['PRNU_TB_CDP'])[-1])
        prnucdp = files.cPickleRead(prnucdp_pick)

        for icol in range(1, Ncols + 1):

            colkey = 'col%03i' % icol

            prnutb = prnucdp['data']['PRNU_%s' % colkey]

            for jCCD, CCDk in enumerate(CCDkeys):
                for kQ, Q in enumerate(self.Quads):

                    ixsel = np.where((prnutb['CCD'] == (jCCD + 1)) & (prnutb['Q'] == (kQ + 1)))
                    prnu_pc_vs[colkey][0, jCCD, kQ] = prnutb['PRNU_PC'].as_matrix()[ixsel][0]

                    G = self.cdps['GAIN'][block][CCDk][Q][0]

                    prnu_fluadu_vs[colkey][0, jCCD, kQ] = prnutb['AVFLUENCE'].as_matrix()[ixsel][0]

                    prnu_fluele_vs[colkey][0, jCCD, kQ] = G * prnu_fluadu_vs[colkey][0, jCCD, kQ]

        for icol in range(1, Ncols + 1):

            colkey = 'col%03i' % icol

            sidd.addColumn(prnu_fluadu_vs[colkey], 'PRNU_FLUADU_%s' % colkey.upper(), IndexCQ)
            sidd.addColumn(prnu_fluele_vs[colkey], 'PRNU_FLUELE_%s' % colkey.upper(), IndexCQ)
            sidd.addColumn(prnu_pc_vs[colkey], 'PRNU_PC_%s' % colkey.upper(), IndexCQ)

        # ADDING REFERENCES TO MASTER FLAT-FIELDS

        tmp_v_C = np.zeros((1, NCCDs), dtype='S50')

        for icol in range(1, Ncols + 1):
            colkey = 'col%03i' % icol

            mfkey_v = tmp_v_C.copy()

            for jCCD, CCDk in enumerate(CCDkeys):

                mfkey = '%s_%s_%s_%i_%s_%s' % (testname, block, session, jrep + 1, colkey, CCDk)

                mfkey_v[0, jCCD] = mfkey

                self.products['MASTERFLATS'][mfkey] = os.path.split(
                    sidd.products['MasterFFs'][colkey][CCDk])[-1]

            sidd.addColumn(mfkey_v, 'MASTERFLATS_%s' % colkey.upper(), IndexC)

        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit

    def _get_extractor_PRNU_fromPT(self, colkey):
        """ """

        def _extract_PRNU_fromPT(PT, block, CCDk, Q):
            ixblock = self.get_ixblock(PT, block)
            column = 'PRNU_PC_%s_%s_Quad%s' % (colkey.upper(), CCDk, Q)

            PRNU = PT[column][ixblock][0]
            return PRNU

        return _extract_PRNU_fromPT

    def _get_XYdict_PRNUFLU(self, PT, Ncols):
        """ """

        x = dict()
        y = dict()

        labelkeys = []

        for icol in range(1, Ncols + 1):
            colkey = 'col%03i' % icol
            labelkeys.append(colkey)

            x[colkey] = []
            y[colkey] = []

            for iCCD, CCD in enumerate(self.CCDs):
                for iQ, Q in enumerate(self.Quads):
                    Xcolumn = 'PRNU_FLUELE_%s_CCD%i_Quad%s' % (colkey.upper(), CCD, Q)
                    Ycolumn = 'PRNU_PC_%s_CCD%i_Quad%s' % (colkey.upper(), CCD, Q)

                    x[colkey] += (PT[Xcolumn] / 1.E3).tolist()
                    y[colkey] += PT[Ycolumn].tolist()

            x[colkey] = np.array(x[colkey])
            y[colkey] = np.array(y[colkey])

        XYdict = dict(x=x, y=y, labelkeys=labelkeys)

        return XYdict

    def _get_MFdict(self, testname, colkey):
        """ """
        PT = self.ParsedTable[testname]

        MFdict = dict()

        for jY in range(self.NCOLS_FPA):
            for iX in range(self.NSLICES_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                locator = self.fpa.FPA_MAP[Ckey]

                block = locator[0]
                CCDk = locator[1]
                flip = locator[2]

                inventoryitem = self.inventory[block][testname][0]

                productspath = os.path.join(inventoryitem['resroot'], 'products')

                ixblock = np.where(PT['BLOCK'] == block)

                cdpkey = PT['MASTERFLATS_%s_%s' % (colkey.upper(), CCDk)][ixblock][0]

                masterfits = os.path.join(productspath, self.products['MASTERFLATS'][cdpkey])

                ccdobj = ccdmod.CCD(infits=masterfits, getallextensions=True, withpover=True,
                                    overscan=20)

                img = ccdobj.extensions[1].data.transpose().copy()

                # MFdict[Ckey] = dict(img = np.zeros((4172,4238),dtype='float32').copy()) # TESTS
                # MFdict[Ckey] = dict(img=img.copy()) # TESTS
                # continue # TESTS

                simg = ndimage.filters.gaussian_filter(img, sigma=5.,
                                                       mode='constant',
                                                       cval=1.)
                esimg = exposure.equalize_hist(simg, nbins=256)

                MFdict[Ckey] = dict(img=self.fpa.flip_img(esimg, flip))

                #ccdobj = None
                #img = None
                #simg = None
                #esimg = None

        return MFdict

    def _get_XYdict_PRNULAM(self):
        """ """

        x = dict()
        y = dict()
        ey = dict()

        labelkeys = ['low PRNU', 'Hi PRNU']

        for labelkey in labelkeys:
            x[labelkey] = []
            y[labelkey] = []
            ey[labelkey] = []

        for testname in self.testnames:

            if testname == 'FLAT01':
                Fluence = 2
            else:
                Fluence = 1

            colkey = 'col%03i' % Fluence

            PT = self.ParsedTable[testname]

            wave = PT['WAVENM'][0]

            x['low PRNU'].append(wave)
            x['Hi PRNU'].append(wave)

            _y_lo = []
            _y_hi = []

            for block in self.flight_blocks:

                ixblock = np.where(PT['BLOCK'] == block)

                for CCD in self.CCDs:

                    Ckey = self.fpa.get_Ckey_from_BlockCCD(block, CCD)

                    locator = self.fpa.FPA_MAP[Ckey]
                    CCDsn = locator[3]
                    CCDbatch = CCDsn[0:5]

                    for Q in self.Quads:

                        PRNUcol = 'PRNU_PC_%s_CCD%i_Quad%s' % (colkey.upper(), CCD, Q)

                        _PRNU = PT[PRNUcol][ixblock][0]

                        if CCDbatch in self.batches_highPRNU:

                            _y_hi.append(_PRNU)
                        else:
                            _y_lo.append(_PRNU)

            y['low PRNU'].append(np.nanmean(_y_lo))
            ey['low PRNU'].append(np.nanstd(_y_lo))
            y['Hi PRNU'].append(np.nanmean(_y_hi))
            ey['Hi PRNU'].append(np.nanstd(_y_hi))

        for labelkey in labelkeys:
            ixorder = np.argsort(x[labelkey])
            x[labelkey] = np.array(x[labelkey])[ixorder]
            y[labelkey] = np.array(y[labelkey])[ixorder]
            ey[labelkey] = np.array(ey[labelkey])[ixorder]

        XYdict = dict(x=x, y=y, ey=ey, labelkeys=labelkeys)

        return XYdict

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['PRNU_vs_WAVE'] = os.path.join(self.figspath,
                                                 'PRNU_vs_WAVLENGTH.png')

        for testname in self.testnames:
            self.figs['PRNU_vs_FLU_%s' % testname] = os.path.join(self.figspath,
                                                                  'PRNU_vs_FLU_%s.png' % testname)

        for testname in self.testnames:

            for icol, colkey in enumerate(self.colkeys[testname]):

                self.figs['PRNU_MAP_%s_%s' % (testname, colkey)] = os.path.join(
                    self.figspath, 'PRNU_MAP_%s_%s.png' % (testname, colkey))
                self.figs['PRNU_MAP_%s_%s_json' % (testname, colkey)] = os.path.join(
                    self.figspath, 'PRNU_MAP_%s_%s.json' % (testname, colkey))

        for testname in self.testnames:
            for colkey in self.colkeys[testname]:
                self.figs['MF_%s_%s' % (testname, colkey)] = os.path.join(
                    self.figspath, 'MF_%s_%s.png' % (testname, colkey))

    def dump_aggregated_results(self):
        """ """

        # PRNU vs. WAVELENGTH

        doPRNUvsWAVE = True
        doMFs = True
        doPRNUvsFLU = True

        if doPRNUvsWAVE:

            XYPLam = self._get_XYdict_PRNULAM()

            self.plot_XY(XYPLam, **dict(
                title='PRNU vs. Wavelength',
                doLegend=True,
                doYErrbars=True,
                xlabel='Wavelength [nm]',
                ylabel='PRNU',
                figname=self.figs['PRNU_vs_WAVE'],
                corekwargs=dict(linestyle='-', marker='o')))
            gc.collect()

        # MASTER FLATS DISPLAY

        if doMFs:

            for testname in self.testnames:

                stestname = st.replace(testname, '_', '\_')

                for colkey in self.colkeys[testname]:
                    # for colkey in [self.colkeys[testname][2]]:

                    colnum = int(colkey[-1])

                    print('MF: %s %s' % (testname, colkey))

                    if testname == 'FLAT01':
                        suptitle = 'FLAT01 800 nm, Fluence %s' % colnum
                    else:
                        _test, _wave = st.split(testname, '_')
                        suptitle = '%s %s nm, Fluence %i' % (_test, _wave, colnum)
                    print(suptitle)

                    MFdict = self._get_MFdict(testname, colkey)
                    MFkwargs = dict(  # suptitle='%s, %s: Master Flat Field' % (stestname, colkey),
                        suptitle=suptitle,
                        figname=self.figs['MF_%s_%s' % (testname, colkey)])

                    self.plot_ImgFPA(MFdict, **MFkwargs)
                    MFdict = None

        # PRNU vs. FLUENCE

        if doPRNUvsFLU:

            for testname in self.testnames:

                NFluCols = len(self.colkeys[testname])
                stestname = st.replace(testname, '_', '\_')

                XYdict = self._get_XYdict_PRNUFLU(self.ParsedTable[testname], NFluCols)

                self.plot_XY(XYdict, **dict(
                    title='%s: PRNU' % (stestname,),
                    doLegend=True,
                    xlabel='Fluence [ke-]',
                    ylabel='PRNU',
                    figname=self.figs['PRNU_vs_FLU_%s' % testname]))

        # PRNU maps

        for testname in self.testnames:

            stestname = st.replace(testname, '_', '\_')

            for icol, colkey in enumerate(self.colkeys[testname]):

                PRNUMAP = self.get_FPAMAP_from_PT(self.ParsedTable[testname],
                                                  extractor=self._get_extractor_PRNU_fromPT(colkey))

                vjson.save_jsonfile(PRNUMAP, self.figs['PRNU_MAP_%s_%s_json' % (testname, colkey)])

                self.plot_SimpleMAP(PRNUMAP, **dict(
                    suptitle='%s [%s]: PRNU' % (stestname, colkey),
                    figname=self.figs['PRNU_MAP_%s_%s' % (testname, colkey)]
                ))
