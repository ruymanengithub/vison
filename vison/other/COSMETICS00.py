#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: COSMETICS00

Cosmetics masks for the detectors: defects in darkness and PR


Created on Wed Dec 05 17:59:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import unittest
from matplotlib.colors import Normalize
from scipy import interpolate, ndimage
import pandas as pd
import string as st

from vison.other import COSaux
from vison.pipe.task import HKKeys
from vison.image import bits
from vison.datamodel import core
from vison.pipe.task import Task
from vison.support import context
from vison.datamodel import scriptic as sc
from vison.datamodel import ccd
from vison.dark.DarkTask import DarkTask
from vison.datamodel import inputs, cdp
from vison.support import utils
from vison.scripts import vis_cosmetics_masker as vcm
# END IMPORT

isthere = os.path.exists


# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


COS_commvalues = dict(program='CALCAMP', test='COSMETICS00',
                      rdmode='fwd_bas',
                      flushes=7, siflsh=1, siflsh_p=500,
                      swellw=context.sumwell['fwd_bas'][0],
                      swelldly=context.sumwell['fwd_bas'][1],
                      inisweep=1,
                      vstart=0, vend=2086,
                      toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                      chinj=0,
                      s_tpump=0,
                      v_tpump=0,
                      e_shuttr=0,
                      mirr_on=0,
                      wave=4,
                      motr_on=0,
                      source='flat')


class COS_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    # manifesto.update(OrderedDict(sorted([
    #    ('N', ([int], 'Number of Frame Acquisitions.')),
    #])))


class COSMETICS00(DarkTask):
    """ """

    inputsclass = COS_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [
            ('check', self.check_data),
            ('masks', self.do_masks),
            ('meta', self.meta)]

        super(COSMETICS00, self).__init__(inputs=inputs, log=log, drill=drill,
                                          debug=debug, cleanafter=cleanafter)
        self.name = 'COSMETICS00'
        self.type = 'Simple'

        self.HKKeys = HKKeys
        self.figdict = COSaux.get_COS_figs()
        self.CDP_lib = COSaux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds COSMETICS00 script structure dictionary.

        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """

        waveflat = 800
        FW_IDflat = self.ogse.get_FW_ID(waveflat)
        FW_IDflatx = int(FW_IDflat[-1])
        exptimeFL = self.ogse.profile['tFWC_flat']['nm%i' % waveflat] * 0.5

        COS_sdict = OrderedDict()

        COS_sdict['col001'] = dict(frames=5, exptime=0,
                                   shuttr=0, mirr_on=0, motr_on=0,
                                   source='flat',
                                   comments='BIAS')

        COS_sdict['col002'] = dict(frames=5, exptime=exptimeFL,
                                   wave=FW_IDflatx,
                                   shuttr=1, mirr_on=0, motr_on=0,
                                   source='flat',
                                   comments='FLAT')

        Ncols = len(COS_sdict.keys())
        COS_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(COS_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
                diffvalues = diffvalues = dict()

        COS_sdict = sc.update_structdict(
            COS_sdict, commvalues, diffvalues)

        return COS_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(COSMETICS00, self).filterexposures(structure, explog, OBSID_lims,
                                                        colorblind=True,
                                                        wavedkeys=wavedkeys)

    def check_data(self):
        """ """

        kwargs = dict(figkeys=[])

        Task.check_data(self, **kwargs)

    def get_checkstats_ST(self, **kwargs):
        """ """
        pass

    def check_metrics_ST(self, **kwargs):
        """ """

        pass

#    def check_data(self):
#        """ """
#        kwargs = dict(figkeys=[])
#        Task.check_data(self, **kwargs)



    def do_masks(self):
        """

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='masks', Title='COSMETICS00: MASKS', level=0)

        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuads = DDindices.shape
        #Quads = DDindices.get_vals('Quad')
        CCDs = DDindices.get_vals('CCD')

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        productspath = self.inputs['subpaths']['products']

        dpath = self.dd.mx['datapath'][0, 0]

        for key in ['DARK', 'FLAT', 'MERGE']:
            self.dd.products[key] = OrderedDict()

        for jCCD, CCDk in enumerate(CCDs):

            seldk = np.where(self.dd.mx['exptime'][:, 0] == 0.)
            OBSID_list_dk = self.dd.mx['ObsID'][seldk].tolist()

            selfl = np.where(self.dd.mx['exptime'][:, 0] > 0.)
            OBSID_list_fl = self.dd.mx['ObsID'][selfl].tolist()

            inputs = OrderedDict()

            inputs['sn_ccd'] = self.inputs['diffvalues']['sn_%s' % CCDk.lower()]
            inputs['CCD'] = jCCD + 1
            inputs['elvis'] = self.elvis
            inputs['tag'] = self.BLOCKID
            inputs['doDarkMask'] = True
            inputs['doFlatMask'] = True
            inputs['doMerge'] = True
            inputs['outpath'] = productspath
            inputs['inputs_dkmask'] = OrderedDict()
            inputs['inputs_dkmask']['datapath'] = dpath
            inputs['inputs_dkmask']['OBSID_list'] = OBSID_list_dk
            inputs['inputs_dkmask']['thresholds'] = [-1000., 50.0]
            inputs['inputs_flmask'] = OrderedDict()
            inputs['inputs_flmask']['datapath'] = dpath
            inputs['inputs_flmask']['OBSID_list'] = OBSID_list_fl
            inputs['inputs_flmask']['thresholds'] = [0.5, 1.5]

            outputs = vcm.run_maskmaker(inputs)

            for key in outputs.keys():
                self.dd.products[key][CCDk] = outputs[key]

        productsstr = self.dd.products.__str__()
        productsstr = st.replace(productsstr, ',', '\n')

        msg = 'Updated products:'

        if self.log is not None:
            self.log.info(msg)
            self.log.info(productsstr)

        if self.report is not None:
            self.report.add_Text(msg)
            self.report.add_Text(productsstr, verbatim=True)

    def _get_NBadCols(self, maskccdobj, Q, iext=-1):
        """ """
        quaddata = maskccdobj.get_quad(Q, canonical=True, extension=iext)
        bX, bY = np.where(quaddata)
        uCol, uColcounts = np.unique(bX, return_counts=True)
        Ncols = int(np.sum(uColcounts > 200))  # (minor) caveat: this also includes non-adjacent pixels

        return Ncols

    def meta(self):
        """

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='COSMETICS00: Masks Displays and Statistics', level=0)

        maskkeys = ['DARK', 'FLAT', 'MERGE']
        NpixImgQuad = self.ccdcalc.NrowsCCD * self.ccdcalc.NcolsCCD

        indices = copy.deepcopy(self.dd.indices)
        CCDs = indices.get_vals('CCD')
        nC = len(CCDs)
        Quads = self.ccdcalc.Quads
        nQ = len(Quads)

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        NP = nC * nQ

        # Display Masks

        MASKS_2PLOT = OrderedDict()
        for maskkey in maskkeys:
            MASKS_2PLOT[maskkey] = OrderedDict()
            for CCDk in CCDs:
                MASKS_2PLOT[maskkey][CCDk] = OrderedDict()
                for Q in Quads:
                    MASKS_2PLOT[maskkey][CCDk][Q] = OrderedDict()

        DEF_TB = OrderedDict()

        DEF_TB = OrderedDict()
        DEF_TB['CCD'] = np.zeros(NP, dtype='int32')
        DEF_TB['Q'] = np.zeros(NP, dtype='int32')

        for maskkey in maskkeys:
            DEF_TB['N_%s' % maskkey] = np.zeros(NP, dtype='int32')
            DEF_TB['NCOLS_%s' % maskkey] = np.zeros(NP, dtype='int32')

        DEF_TB['PIXLOST'] = np.zeros(NP, dtype='float32')

        normfunction = Normalize(vmin=1. - 0.00001, vmax=1., clip=True)

        if not self.drill:

            for jCCD, CCDk in enumerate(CCDs):
                for jQ, Q in enumerate(Quads):
                    kk = jCCD * nQ + jQ
                    DEF_TB['CCD'][kk] = jCCD + 1
                    DEF_TB['Q'][kk] = jQ + 1

            for maskkey in maskkeys:

                for jCCD, CCDk in enumerate(CCDs):

                    MSKfits = self.dd.products[maskkey][CCDk]

                    MSKccdobj = ccd.CCD(MSKfits)
                    iMext = MSKccdobj.extnames.index('MASK')

                    for jQ, Q in enumerate(Quads):

                        kk = jCCD * nQ + jQ

                        kQ_Ndefects = MSKccdobj.get_stats(Q, sector='img', statkeys=['sum'],
                                                          ignore_pover=True,
                                                          extension=iMext)[0]

                        DEF_TB['N_%s' % maskkey][kk] = kQ_Ndefects
                        DEF_TB['PIXLOST'][kk] = kQ_Ndefects / NpixImgQuad * 100.
                        DEF_TB['NCOLS_%s' %
                               maskkey][kk] = self._get_NBadCols(MSKccdobj, Q, iext=iMext)

                        qdata = MSKccdobj.get_quad(Q, canonical=False, extension=iMext).copy()
                        qdata = 1. - qdata  # inversion
                        sqdata = ndimage.filters.gaussian_filter(qdata, sigma=10.,
                                                                 mode='constant',
                                                                 cval=0.)
                        #if maskkey == 'MERGE':stop()
                        MASKS_2PLOT[maskkey][CCDk][Q]['img'] = sqdata.T.copy()

                self.figdict['Masks_%s' % maskkey][1]['data'] = MASKS_2PLOT[maskkey].copy()
                self.figdict['Masks_%s' % maskkey][1]['meta']['corekwargs']['norm'] = normfunction

        # REPORT THE DEFECTS STATS

        DEF_TB_dddf = OrderedDict()
        DEF_TB_dddf['DEFECTS'] = pd.DataFrame.from_dict(DEF_TB)

        def_tb_cdp = self.CDP_lib['DEF_TB']
        def_tb_cdp.path = self.inputs['subpaths']['products']
        def_tb_cdp.ingest_inputs(
            data=DEF_TB_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy()
        )

        def_tb_cdp.init_wb_and_fillAll(header_title='%s: DEFECTS COUNTS TABLE' %
                                       (self.inputs['test'],))
        self.save_CDP(def_tb_cdp)
        self.pack_CDP_to_dd(def_tb_cdp, 'DEF_TB_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def fi(x): return '%i' % x

            def fpc(x): return '%.2f %%' % x

            cov_formatters = [fccd, fq, fi, fi, fi, fi, fpc]
            columns = ['CCD', 'Q', 'N_DARK', 'N_FLAT', 'N_MERGE', 'NCOLS_MERGE', 'PIXLOST']

            caption = '%s: DEFECTS TABLE.' % \
                (self.inputs['test'],)
            nicecaption = st.replace(caption, '_', '\_')
            Ptex = def_tb_cdp.get_textable(sheet='DEFECTS',
                                           caption=nicecaption,
                                           fitwidth=True,
                                           tiny=True,
                                           formatters=cov_formatters,
                                           columns=columns)
            self.report.add_Text(Ptex)

        # DO THE PLOTS

        if self.report is not None:

            MASKfigkeys = ['Masks_%s' % key for key in maskkeys]
            self.addFigures_ST(figkeys=MASKfigkeys,
                               dobuilddata=False)

        self.canbecleaned = True
