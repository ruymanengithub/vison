#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Metatests: DARKS

Created on Mon Feb 10 16:07:00 2020

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
from vison.image import cosmetics

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
    'std_ove',
    'chk_flu_img']


class MetaDark(MetaCal):
    """ """

    def __init__(self, **kwargs):
        """ """

        super(MetaDark, self).__init__(**kwargs)

        self.testnames = ['DARK01']
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()

        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['MD_PROFILES'] = OrderedDict()
        self.products['MASTERDARKS'] = OrderedDict()

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

        dksignal_v = tmp_v_CQ.copy()
        hotpixels_v = tmp_v_CQ.copy()

        productspath = os.path.join(inventoryitem['resroot'], 'products')
        profspath = os.path.join(inventoryitem['resroot'], 'profiles')

        dktbcdp_pick = os.path.join(productspath, os.path.split(sidd.products['DARK_TB_CDP'])[-1])
        dktbcdp = files.cPickleRead(dktbcdp_pick)

        profs_pick = os.path.join(profspath, os.path.split(sidd.products['MD_PROFILES'])[-1])
        profs_dict = files.cPickleRead(profs_pick)

        profkey = '%s_%s_%s_%i' % (testname, block, session, jrep + 1)

        self.products['MD_PROFILES'][profkey] = profs_dict.copy()

        profskeys_v = np.zeros((1), dtype='U50')

        profskeys_v[0] = profkey

        nQ = len(self.Quads)

        for iCCD, CCDk in enumerate(CCDkeys):

            for kQ, Q in enumerate(self.Quads):

                kk = iCCD * nQ + kQ

                dksignal_v[0, iCCD, kQ] = dktbcdp['data']['DARK']['AVSIGNAL'][kk]
                hotpixels_v[0, iCCD, kQ] = dktbcdp['data']['DARK']['N_HOT'][kk]


        sidd.addColumn(dksignal_v, 'DK_SIGNAL', IndexCQ)
        sidd.addColumn(hotpixels_v, 'DK_N_HOT', IndexCQ)

        sidd.addColumn(profskeys_v, 'MDPROFS_KEY', IndexS)

        # ADDING REFERENCES TO MASTER DARKS

        tmp_v_C = np.zeros((1, NCCDs), dtype='U50')

        mdkey_v = tmp_v_C.copy()

        for jCCD, CCDk in enumerate(CCDkeys):

            mdkey = '%s_%s_%s_%i_%s' % (testname, block, session, jrep + 1, CCDk)

            mdkey_v[0, jCCD] = mdkey

            self.products['MASTERDARKS'][mdkey] = os.path.split(
                sidd.products['MasterDKs'][CCDk])[-1]

        sidd.addColumn(mdkey_v, 'MASTERDARKS', IndexC)


        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit


    def _extract_DKSIG_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)
        column = 'DK_SIGNAL_%s_Quad%s' % (CCDk, Q)
        DKSIG = PT[column][ixblock][0]

        return DKSIG




    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['DK_SIGNAL_MAP'] = os.path.join(self.figspath,\
                                        'DK_SIGNAL_MAP.png')

        self.figs['HOTPIX_COUNT_MAP'] = \
                os.path.join(self.figspath, 'HOTPIX_COUNT_MAP.png')

        self.figs['MD'] = os.path.join(self.figspath, \
            'MASTER_DARK.png')

        profdirs = ['hor','ver']

        for profdir in profdirs:
            self.figs['PROFS_%s' % (profdir,)] =  os.path.join(self.figspath,
                                'DK_PROFS_%s.png' % (profdir,))

    def init_outcdpnames(self):
        """ """

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        self.outcdps['DK_SIGNAL'] = 'DARK01_DK_SIGNAL_ADU_MAP.json'

    def _get_MDdict(self):
        """ """

        PT = self.ParsedTable['DARK01']

        MDdict = dict()

        for jY in range(self.NCOLS_FPA):
            for iX in range(self.NSLICES_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                locator = self.fpa.FPA_MAP[Ckey]

                block = locator[0]
                CCDk = locator[1]
                flip = locator[2]

                inventoryitem = self.inventory[block]['DARK01'][0]

                productspath = os.path.join(inventoryitem['resroot'], 'products')

                ixblock = np.where(PT['BLOCK'] == block)

                cdpkey = PT['MASTERDARKS_%s' % CCDk][ixblock][0]

                masterfits = os.path.join(productspath, self.products['MASTERDARKS'][cdpkey])

                ccdobj = ccdmod.CCD(infits=masterfits, getallextensions=True, withpover=True,
                                    overscan=20)

                img = ccdobj.extensions[1].data.transpose().copy()


                # coarse hot pixel masking (just for display)
                medval = np.median(img)
                img[np.where(img>medval*50)] = medval

                # image smoothing

                simg = ndimage.filters.gaussian_filter(img, sigma=5.,
                                                       mode='constant',
                                                       cval=1.)

                esimg = exposure.equalize_hist(simg, nbins=256)

                MDdict[Ckey] = dict(img=self.fpa.flip_img(esimg, flip))

                #ccdobj = None
                #img = None
                #simg = None
                #esimg = None

        return MDdict


    def _get_XYdict_PROFS(self,proftype):
        """ """

        x = dict()
        y = dict()

        labelkeys = []

        PT = self.ParsedTable['DARK01']

        for block in self.flight_blocks:

            ixsel = np.where(PT['BLOCK'] == block)
            prof_key = PT['MDPROFS_KEY'][ixsel][0]

            i_Prof = self.products['MD_PROFILES'][prof_key].copy()

            for iCCD, CCD in enumerate(self.CCDs):
                CCDk = 'CCD%i' % CCD

                for kQ, Q in enumerate(self.Quads):

                    pkey = '%s_%s_%s' % (block, CCDk, Q)

                    _pcq = i_Prof['data'][proftype][CCDk][Q].copy()

                    x[pkey] = _pcq['x'].copy()
                    y[pkey] = _pcq['y'] - np.nanmedian(_pcq['y'])

                    labelkeys.append(pkey)

        Pdict = dict(x=x, y=y, labelkeys=labelkeys)

        return Pdict

    def _extract_logNPHOT_fromPT(self, PT, block, CCDk, Q):
        """ """
        ixblock = self.get_ixblock(PT, block)
        column = 'DK_N_HOT_%s_Quad%s' % (CCDk, Q)
        logNHP = np.log10(max(PT[column][ixblock][0],1.))
        return logNHP

    def dump_aggregated_results(self):
        """ """



        if self.report is not None:
            self.report.add_Section(keyword='dump',
                Title='Aggregated Results', level=0)

            self.add_DataAlbaran2Report()

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        # Dark Signal map, ADUs

        DKSIGMAP = self.get_FPAMAP_from_PT(
            self.ParsedTable['DARK01'],
            extractor=self._extract_DKSIG_fromPT)

        dksig_header = OrderedDict()
        dksig_header['title'] = 'DARK SIGNAL MAP'
        dksig_header.update(CDP_header)

        dksig_cdp = cdp.Json_CDP(rootname=self.outcdps['DK_SIGNAL'],
                              path=self.cdpspath)
        dksig_cdp.ingest_inputs(data=DKSIGMAP,
                             header = dksig_header,
                             meta=dict(units='ADU'))
        dksig_cdp.savehardcopy()

        figkey1 = 'DK_SIGNAL_MAP'
        figname1 = self.figs[figkey1]

        self.plot_SimpleMAP(DKSIGMAP, **dict(
                suptitle="DARK01: 'DARK' SIGNAL [ADU]",
                ColorbarText='ADU',
                figname=figname1))


        if self.report is not None:
            self.addFigure2Report(figname1,
                figkey=figkey1, 
                caption='DARK01: "Dark" Signal (ADU) integrated over 565 seconds. '+\
                'This is mostly stray-light, in fact. See text].', 
                texfraction=0.7)

            dksigcdpdict = dict(
                caption='DARK01: "Dark" Signal (ADU) integrated over 565 seconds. '+\
                '[Mostly stray-light, in fact. See text].',
                valformat='%.2f')

            ignore = self.add_StdQuadsTable2Report( 
                Matrix=DKSIGMAP,
                cdpdict=dksigcdpdict)

        # MASTER DARK DISPLAY 

        figkey2 = 'MD'
        figname2 = self.figs[figkey2]

        MDdict = self._get_MDdict()
        MDkwargs = dict(  # suptitle='%s, %s: Master Flat Field' % (stestname, colkey),
            suptitle='DARK01: Master "Dark" [stray-light]',
            figname=figname2)

        self.plot_ImgFPA(MDdict, **MDkwargs)

        if self.report is not None:
            self.addFigure2Report(figname2,
                figkey=figkey2,
                caption='DARK01: Master Dark [ADU]. Dominated by stray-light. See text.', 
                texfraction=0.7)

        MDdict = None

        # HOT PIXELS MAP [PENDING]

        logNHPMAP = self.get_FPAMAP_from_PT(self.ParsedTable['DARK01'],
                                              extractor=self._extract_logNPHOT_fromPT)

        figkey3 = 'HOTPIX_COUNT_MAP'
        figname3 = self.figs[figkey3]

        self.plot_SimpleMAP(logNHPMAP, **dict(
                suptitle='DARK01: log(NR) of Hot Pixels',
                figname=figname3,
                ColorbarText='log(N)'
            ))

        if self.report is not None:

            captemp = 'log(Number) of hot pixels in each CCD quadrant of the FPA. '

            self.addFigure2Report(figname3, 
                figkey=figkey3, 
                caption=captemp, 
                texfraction=0.7)


        # Vertical and Horizonal AVG PROFILES

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


        for proftype in proftypes:

            XY_profs = self._get_XYdict_PROFS(proftype=proftype)

            figkey4 = 'PROFS_%s' % (proftype,)

            if proftype == 'hor':
                ylim = [-20,50]
            elif proftype == 'ver':
                ylim = [-20,30]

            profkwargs = dict(
                title='DARK01: Avg. profiles, direction: %s' % (proftype,),
                doLegend=False,
                xlabel=xlabels_profs[proftype],
                ylabel=r'$ADU$',
                ylim=ylim,
                xlim=None,
                figname=self.figs[figkey4],
                corekwargs=pointcorekwargs)

            self.plot_XY(XY_profs, **profkwargs)

            if proftype == 'ver':
                captemp = 'DARK01: Stacked profiles of Master Dark quadrant'+\
                ' images in "parallel" direction. Median value of profile has been subtracted for '+\
                'clarity. Each colour corresponds to a different block (each with 3x4 quadrants).'
            elif proftype == 'hor':
                captemp = 'DARK01: Stacked profiles of Master Dark quadrant'+\
                ' images in "serial" direction. Median value of profile has been subtracted for '+\
                'clarity. Each colour corresponds to a'+\
                ' different block (each with 3x4 quadrants).'

            if self.report is not None:

                self.addFigure2Report(self.figs[figkey4], 
                    figkey=figkey4, 
                    caption= captemp,
                    texfraction=0.7)

