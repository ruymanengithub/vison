#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: DARK01

"Dark Current" analysis script

Created on Tue Aug 29 17:21:00 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
from copy import deepcopy
from collections import OrderedDict
import copy
from scipy import ndimage
from matplotlib.colors import Normalize
import pandas as pd

from vison.pipe.task import HKKeys
from vison.support import context, utils
#from vison.pipe import lib as pilib
#from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import cdp
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#from vison.pipe.task import Task
from DarkTask import DarkTask
from vison.dark import D01aux, darkaux
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#         'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


DARK01_commvalues = dict(program='CALCAMP', test='DARK01',
                         flushes=7, siflsh=1, siflsh_p=500,
                         inisweep=1,
                         vstart=0, vend=2086,
                         shuttr=0, e_shuttr=0,
                         wave=4,
                         source='flat',
                         comments='DARK')


class DARK01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('N', ([int], 'Number of Frame Acquisitions.')),
        ('exptime', ([float], 'Exposure time.')),
    ])))


Flu_lims = OrderedDict(CCD1=OrderedDict(E=[-2., 5.]))
for Q in ['F', 'G', 'H']:
    Flu_lims['CCD1'][Q] = copy.deepcopy(Flu_lims['CCD1']['E'])
for i in [2, 3]:
    Flu_lims['CCD%i' % i] = copy.deepcopy(Flu_lims['CCD1'])


class DARK01(DarkTask):
    """ """

    inputsclass = DARK01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('meta', self.stack_analysis)]
        super(DARK01, self).__init__(inputs=inputs, log=log, drill=drill,
                                     debug=debug, cleanafter=cleanafter)
        self.name = 'DARK01'
        self.type = 'Simple'

        self.HKKeys = HKKeys
        self.figdict = D01aux.get_D01figs()
        self.CDP_lib = D01aux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs',
                                       ccdpickles='ccdpickles',
                                       profiles='profiles',
                                       products='products')

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = dict(N=4, exptime=565)

    def set_perfdefaults(self, **kwargs):
        super(DARK01, self).set_perfdefaults(**kwargs)
        self.perfdefaults['Flu_lims'] = Flu_lims

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds DARK01 script structure dictionary.

        :param diffvalues: dict, opt, differential values.

        """

        N = self.inputs['N']
        exptime = self.inputs['exptime']
        DARK01_sdict = dict(col001=dict(frames=N, exptime=exptime))

        Ncols = len(DARK01_sdict.keys())
        DARK01_sdict['Ncols'] = Ncols

        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(DARK01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
                diffvalues = diffvalues = dict()

        DARK01_sdict = sc.update_structdict(
            DARK01_sdict, commvalues, diffvalues)

        return DARK01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(DARK01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                   wavedkeys=wavedkeys)

    def prep_data(self):
        """

        DARK01: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            [BIAS SUBTRACTION]
            cosmetics masking

        """
        super(DARK01, self).prepare_images(
            doExtract=True, doBadPixels=True,
            doMask=True, doOffset=True, doBias=True, doFF=False)

    def stack_analysis(self):
        """

        **METACODE**

        ::

            f. each CCD:
                f. e. Q:
                    stack all ObsIDs to produce Master Dark
                    produce mask of hot pixels / columns
                    count hot pixels / columns
                    measure average profile along rows
                    measure average profile along cols

            plot average profiles of Master Bias f. each CCD,Q
            show Master Dark (images), include in report
            report stats of defects, include in report
            save name of MasterDark to DataDict, report
            save name of Defects in Darkness Mask to DD, report


        """

        if self.report is not None:
            self.report.add_Section(
                keyword='MasterDK', Title='Master Darks', level=0)

        OnTests = False  # True on TESTS

        settings = dict(
            ID=self.ID,
            BLOCKID=self.BLOCKID,
            CHAMBER=self.CHAMBER
        )

        indices = copy.deepcopy(self.dd.indices)

        CCDs = indices.get_vals('CCD')
        nC = len(CCDs)
        Quads = indices.get_vals('Quad')
        nQ = len(Quads)

        dpath = self.inputs['subpaths']['ccdpickles']
        productspath = self.inputs['subpaths']['products']
        profilespath = self.inputs['subpaths']['profiles']

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        NP = nC * nQ

        DK_TB = OrderedDict()

        DK_TB = OrderedDict()
        DK_TB['CCD'] = np.zeros(NP, dtype='int32')
        DK_TB['Q'] = np.zeros(NP, dtype='int32')
        DK_TB['AVSIGNAL'] = np.zeros(NP, dtype='float32')
        DK_TB['N_HOT'] = np.zeros(NP, dtype='int32')

        MDK_2PLOT = OrderedDict()
        for CCDk in CCDs:
            MDK_2PLOT[CCDk] = OrderedDict()
            for Q in Quads:
                MDK_2PLOT[CCDk][Q] = OrderedDict()
        MDK_p5s = []
        MDK_p95s = []

        # 1D Profiles of Master Dark

        profs1D2plot = OrderedDict()
        profs1D2plot['hor'] = OrderedDict()
        profs1D2plot['ver'] = OrderedDict()

        for CCDk in CCDs:
            for tag in ['hor', 'ver']:
                profs1D2plot[tag][CCDk] = OrderedDict()
            for Q in Quads:
                for tag in ['hor', 'ver']:
                    profs1D2plot[tag][CCDk][Q] = OrderedDict()
                    profs1D2plot[tag][CCDk][Q]['x'] = np.arange(100, dtype='float32')
                    profs1D2plot[tag][CCDk][Q]['y'] = np.zeros(100, dtype='float32')

        if not self.drill:

            self.dd.products['MasterDKs'] = OrderedDict()

            def _pack_profs(CQdict, prof):
                """ """
                _x = prof.data['x'].copy()
                xorder = np.argsort(prof.data['x'])
                _y = prof.data['y'][xorder].copy()
                _x = np.arange(len(_y))

                CQdict['x'] = _x.copy()
                CQdict['y'] = _y.copy()

                return CQdict

            vfullinpath_adder = utils.get_path_decorator(dpath)

            for jCCD, CCDk in enumerate(CCDs):

                DKname = 'EUC_DK_ROE1_%s.fits' % (CCDk,)
                DKpath = os.path.join(productspath, DKname)

                DKlist = np.squeeze(self.dd.mx['ccdobj_name'][:, jCCD]).copy()

                DKlist = vfullinpath_adder(DKlist, extension='pick')

                # MISSING: proper defects and useful area masking
                #   (mask-out pre/over scans)

                if OnTests:
                    DKlist = DKlist[0:2]

                vstart = self.dd.mx['vstart'][0, jCCD]
                vend = self.dd.mx['vend'][0, jCCD]

                jsettings = copy.deepcopy(settings)
                jsettings['CCDSerial'] = self.inputs['diffvalues']['sn_ccd%i' % (jCCD + 1)]

                jsettings['CCDTempTop'] = self.dd.mx['HK_CCD%i_TEMP_T' % (jCCD + 1,)][:].mean()
                jsettings['CCDTempBot'] = self.dd.mx['HK_CCD%i_TEMP_B' % (jCCD + 1,)][:].mean()

                for kQ, Q in enumerate(Quads):
                    jsettings['AVFLU_%s' % Q] = \
                        np.nanmedian(self.dd.mx['chk_flu_img'][:, jCCD, kQ])

                darkaux.produce_MasterDark(DKpath, ccdpickList=DKlist, mask=None,
                                           settings=jsettings)  # COMMENTED ON TESTS

                self.dd.products['MasterDKs'][CCDk] = DKpath

                DK = darkaux.DarkCDP(DKpath)

                iDext = DK.extnames.index('DARK')

                Quads = DK.Quads

                for jQ, Q in enumerate(Quads):

                    kk = jCCD * nQ + jQ

                    DK_TB['CCD'][kk] = jCCD + 1
                    DK_TB['Q'][kk] = jQ + 1
                    DK_TB['N_HOT'][kk] = 0  # PENDING!
                    DK_TB['AVSIGNAL'][kk] = jsettings['AVFLU_%s' % Q]

                    qdata = DK.get_quad(Q, canonical=False, extension=iDext).copy()
                    sqdata = ndimage.filters.gaussian_filter(qdata, sigma=5.,
                                                             mode='constant',
                                                             cval=1.)
                    MDK_2PLOT[CCDk][Q]['img'] = sqdata.transpose()
                    MDK_p5s.append(np.percentile(sqdata, 5))
                    MDK_p95s.append(np.percentile(sqdata, 95))

                    hor1Dprof = DK.get_1Dprofile(Q=Q, orient='hor',
                                                 area='all', stacker='mean',
                                                 vstart=vstart, vend=vend, extension=iDext)

                    profs1D2plot['hor'][CCDk][Q] = _pack_profs(
                        profs1D2plot['hor'][CCDk][Q], hor1Dprof)

                    ver1Dprof = DK.get_1Dprofile(Q=Q, orient='ver',
                                                 area='all', stacker='mean',
                                                 vstart=vstart, vend=vend, extension=iDext)

                    profs1D2plot['ver'][CCDk][Q] = _pack_profs(
                        profs1D2plot['ver'][CCDk][Q], ver1Dprof)

        # PLOTTING 1D PROFILES OF MASTER BIAS

        self.figdict['D01meta_prof1D_hor'][1]['data'] = profs1D2plot['hor'].copy()
        self.figdict['D01meta_prof1D_ver'][1]['data'] = profs1D2plot['ver'].copy()

        if self.report is not None:
            self.addFigures_ST(figkeys=['D01meta_prof1D_hor',
                                        'D01meta_prof1D_ver'],
                               dobuilddata=False)

        # SAVING 1D PROFILES OF MASTER DARK as a CDP

        MD_profiles_cdp = self.CDP_lib['MD_profiles']
        MD_profiles_cdp.header = CDP_header.copy()
        MD_profiles_cdp.path = profilespath
        MD_profiles_cdp.data = profs1D2plot.copy()

        self.save_CDP(MD_profiles_cdp)
        self.pack_CDP_to_dd(MD_profiles_cdp, 'MD_PROFILES')

        # DISPLAYING THE MASTER DARK FRAMES

        self.figdict['D01meta_MasterDark_2D'][1]['data'] = MDK_2PLOT.copy()

        # UPDATING scaling based on data

        if len(MDK_p5s) > 0:
            normfunction = Normalize(vmin=np.min(MDK_p5s),\
                vmax=np.max(MDK_p95s), clip=False)
        else:
            normfunction = False

        self.figdict['D01meta_MasterDark_2D'][1]['meta']['corekwargs']['norm'] = normfunction

        if self.report is not None:
            self.addFigures_ST(figkeys=['D01meta_MasterDark_2D'],
                               dobuilddata=False)

        # REPORT DARK results

        DK_TB_dddf = OrderedDict()
        DK_TB_dddf['DARK'] = pd.DataFrame.from_dict(DK_TB)

        dk_tb_cdp = self.CDP_lib['DARK_TB']
        dk_tb_cdp.path = self.inputs['subpaths']['products']
        dk_tb_cdp.ingest_inputs(
            data=DK_TB_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy()
        )

        dk_tb_cdp.init_wb_and_fillAll(header_title='DARK01 TABLE')
        self.save_CDP(dk_tb_cdp)
        self.pack_CDP_to_dd(dk_tb_cdp, 'DARK_TB_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def fE(x): return '%.2E' % x

            def fi(x): return '%i' % x

            cov_formatters = [fccd, fq, fE, fi]

            caption = 'DARK Results Table.'
            DKtex = dk_tb_cdp.get_textable(sheet='DARK',
                                           caption=caption,
                                           fitwidth=True,
                                           tiny=True,
                                           formatters=cov_formatters)
            self.report.add_Text(DKtex)

        self.canbecleaned = True
