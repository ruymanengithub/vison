# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BF01

Brighter-Fatter Analysis
   Using data from test PTC01 or PTC02

Created on Wed Mar 7 10:57:00 2018

:author: raf

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
import string as st
from collections import OrderedDict
import pandas as pd
import multiprocessing as mp

from vison.pipe.task import HKKeys
from vison.support import context
#from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import HKtools
from vison.datamodel import core, ccd, cdp
from vison.image import calibration
from . import ptc as ptclib
from vison.image import performance
from .FlatTask import FlatTask
from .PTC0X import PTC0X
from vison.pipe.task import Task
from vison.datamodel import inputs
from . import BF01aux
from vison.analysis import Guyonnet15 as G15
from vison.image import covariance as covlib
from sklearn import linear_model

from vison.support import utils
from vison.support.files import cPickleRead, cPickleDumpDictionary
# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


def process_one_fluence_covmaps(q, dd, dpath, CCDs, jCCD, ku, ulabels, Npix,
                                vstart, vend):
    """ """

    labels = dd.mx['label'][:, 0].copy()

    vfullinpath_adder = utils.get_path_decorator(dpath)

    #taskcounter = 0

    ulabel = ulabels[ku]

    six = np.where(labels == ulabel)

    ccdobjNamesList = vfullinpath_adder(dd.mx['ccdobj_name'][six[0], jCCD], 'pick')

    print(('%s: column %i/%i; %i OBSIDs' % (CCDs[jCCD], ku + 1, len(ulabels), len(ccdobjNamesList))))

    ccdobjList = [cPickleRead(item) for item in ccdobjNamesList]

    icovdict = covlib.get_cov_maps(
        ccdobjList, Npix=Npix, vstart=vstart, vend=vend,
        doTest=False, debug=False)
    q.put([jCCD, ku, icovdict])


class BF01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
        ('Npix', ([int], 'Number of Pixels (linear) to consider for Covariance Matrix')),
        ('surrogate', ([str], 'Test to use as surrogate'))
    ])))


def fit_BF(X, Y):
    """ """

    xfit = np.linspace(0, 2**16, 2)

    ixval = np.where(np.isfinite(X) & np.isfinite(Y))

    def_res = dict(xfit=xfit.copy(),
                   yfit=xfit * 0.,
                   intercept=0.,
                   slope=0.,
                   fwhm_hwc=-1.)

    if len(ixval[0]) < 3:
        return def_res

    try:
        ransac = linear_model.RANSACRegressor()
        ransac.fit(np.expand_dims(X[ixval], 1), np.expand_dims(Y[ixval], 1))
    #predictor = ransac.predict
    except BaseException:
        return def_res

    slope = ransac.estimator_.coef_[0][0]
    intercept = ransac.estimator_.intercept_[0]
    fwhm_hwc = ransac.predict(2.**16 / 2.)[0][0]

    yfit = np.squeeze(ransac.predict(np.expand_dims(xfit, 1)))

    res = dict(xfit=xfit.copy(), yfit=yfit.copy(), intercept=intercept,
               slope=slope * 1.e4,
               fwhm_hwc=fwhm_hwc)

    return res


class BF01(PTC0X):
    """ """

    inputsclass = BF01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        #super(BF01, self).__init__(inputs, log, drill, debug)
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract_COV', self.extract_COV),
                         ('extract_BF', self.extract_BF),
                         ('meta', self.meta_analysis)]
        FlatTask.__init__(self, inputs=inputs, log=log, drill=drill, debug=debug,
                          cleanafter=cleanafter)
        #self.inputs['todo_flags'] = self.init_todo_flags()
        # if 'todo_flags' in inputs:
        #    self.inputs['todo_flags'].update(inputs['todo_flags'])

        self.name = 'BF01'
        #self.type = 'Simple'

        self.HKKeys = HKKeys
        self.CDP_lib = BF01aux.get_CDP_lib()
        self.figdict = BF01aux.gt_BF01figs(self.inputs['surrogate'])
        self.inputs['subpaths'] = dict(figs='figs',
                                       ccdpickles='ccdpickles',
                                       covariance='covariance',
                                       kernels='kernels',
                                       products='products')

    def set_inpdefaults(self, **kwargs):
        """ """

        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_inpdefaults(**_kwargs)
        self.inpdefaults['test'] = kwargs['test']
        self.inpdefaults['surrogate'] = kwargs['surrogate']
        self.inpdefaults['Npix'] = 5

    def set_perfdefaults(self, **kwargs):
        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_perfdefaults(**_kwargs)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds PTC0X script structure dictionary.

        #:param exptimes: list of ints [ms], exposure times.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 800 nm.
        :param diffvalues: dict, opt, differential values.

        """
        BF01_sdict = super(BF01, self).build_scriptdict(diffvalues=diffvalues, elvis=elvis)
        Ncols = BF01_sdict['Ncols']
        for i in range(1, Ncols + 1):
            BF01_sdict['col%03i' % i]['test'] = self.inputs['surrogate']
        return BF01_sdict
        # raise NotImplementedError(
        #    "%s: This Task does not build a script, it uses data from another test" % self.name)

    def filterexposures(self, structure, explog, OBSID_lims):
        """

        """
        return Task.filterexposures(self, structure, explog, OBSID_lims,
                                    wavedkeys=['motr_siz'], colorblind=False,
                                    surrogate=self.inputs['surrogate'])

    def check_data(self):

        kwargs = dict(figkeys=['BF01checks_offsets', 'BF01checks_deltaoff',
                               'BF01checks_stds',
                               'BF01checks_flu', 'BF01checks_imgstd'])

        Task.check_data(self, **kwargs)

    def prepare_images(self):
        Task.prepare_images(self, doExtract=True, doBadPixels=True,
                            doMask=True, doOffset=True, doBias=False, doFF=False)

    def extract_COV(self):
        """

        Performs basic analysis of images:
            - extracts COVARIANCE matrix for each fluence

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extractCOV', Title='Covariance-Matrices Extraction', level=0)

        Npix = self.inputs['Npix']

        # labels should be the same accross CCDs. PATCH.
        labels = self.dd.mx['label'][:, 0].copy()
        ulabels = np.unique(labels)
        nL = len(ulabels)

        # vstart and vend should be the same for all OBSIDs in test
        vstart = self.dd.mx['vstart'][0, 0]
        vend = min(self.dd.mx['vend'][0, 0], self.ccdcalc.NrowsCCD)

        indices = copy.deepcopy(self.dd.indices)
        nObs, nC, nQ = indices.shape
        CCDs = indices.get_vals('CCD')
        Quads = indices.get_vals('Quad')

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        dpath = self.inputs['subpaths']['ccdpickles']
        covpath = self.inputs['subpaths']['covariance']
        #prodspath = self.inputs['subpaths']['products']

        profscov_1D = self.CDP_lib['PROFSCOV1D']
        profscov_1D.header = CDP_header.copy()
        profscov_1D.path = covpath
        profscov_1D.data = OrderedDict()

        profscov_1D.data['hor'] = OrderedDict()
        profscov_1D.data['ver'] = OrderedDict()

        for tag in ['hor', 'ver']:
            profscov_1D.data[tag] = OrderedDict()
            for CCDk in CCDs:
                profscov_1D.data[tag][CCDk] = OrderedDict()
                for Q in Quads:
                    profscov_1D.data[tag][CCDk][Q] = OrderedDict()
                    profscov_1D.data[tag][CCDk][Q]['x'] = OrderedDict()
                    profscov_1D.data[tag][CCDk][Q]['y'] = OrderedDict()

        NP = nC * nQ * nL

        COV_dd = OrderedDict()
        COV_dd['CCD'] = np.zeros(NP, dtype='int32')
        COV_dd['Q'] = np.zeros(NP, dtype='int32')
        COV_dd['col'] = np.zeros(NP, dtype='int32')
        COV_dd['av_mu'] = np.zeros(NP, dtype='float32')
        COV_dd['av_var'] = np.zeros(NP, dtype='float32')
        COV_dd['COV_00'] = np.zeros(NP, dtype='float32')
        COV_dd['COV_01'] = np.zeros(NP, dtype='float32')
        COV_dd['COV_10'] = np.zeros(NP, dtype='float32')
        COV_dd['COV_11'] = np.zeros(NP, dtype='float32')

        self.dd.products['COV'] = OrderedDict()
        for CCDk in CCDs:
            self.dd.products['COV'][CCDk] = OrderedDict()

        if not self.drill:

            # doTest=False

            arglist = []

            mgr = mp.Manager()
            queue = mgr.Queue()

            for jCCD, CCDk in enumerate(CCDs):

                for ku, ulabel in enumerate(ulabels):

                    arglist.append([queue, self.dd, dpath, CCDs, jCCD, ku, ulabels,
                                    Npix, vstart, vend])

            pool = mp.Pool(processes=self.processes)

            for ii in range(len(arglist)):
                pool.apply_async(process_one_fluence_covmaps, args=arglist[ii])
            pool.close()
            pool.join()

            replies = []
            while not queue.empty():
                replies.append(queue.get())

            for reply in replies:

                jCCD, ku, icovdict = reply
                CCDk = CCDs[jCCD]
                ulabel = ulabels[ku]

                self.dd.products['COV'][CCDk][ulabel] = copy.deepcopy(
                    icovdict)

                for lQ, Q in enumerate(Quads):
                    jj = jCCD * (nQ * nL) + ku * nQ + lQ

                    COV_dd['CCD'][jj] = jCCD + 1
                    COV_dd['Q'][jj] = lQ + 1
                    COV_dd['col'][jj] = int(st.replace(ulabel, 'col', ''))
                    COV_dd['av_mu'][jj] = icovdict['av_mu'][Q]
                    COV_dd['av_var'][jj] = icovdict['av_var'][Q]
                    COV_dd['COV_00'][jj] = icovdict['av_covmap'][Q][0, 0]
                    COV_dd['COV_01'][jj] = icovdict['av_covmap'][Q][0, 1]
                    COV_dd['COV_10'][jj] = icovdict['av_covmap'][Q][1, 0]
                    COV_dd['COV_11'][jj] = icovdict['av_covmap'][Q][1, 1]

                    profscov_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                        np.arange(Npix - 1)
                    profscov_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                        icovdict['av_covmap'][Q][0, 1:].copy()

                    profscov_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                        np.arange(Npix - 1)
                    profscov_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                        icovdict['av_covmap'][Q][1:, 0].copy()

        else:

            for jCCD, CCDk in enumerate(CCDs):

                for ku, ulabel in enumerate(ulabels):

                    for lQ, Q in enumerate(Quads):

                        profscov_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                            np.arange(Npix - 1)
                        profscov_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                            np.arange(Npix - 1)

                        profscov_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                            np.arange(Npix - 1)
                        profscov_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                            np.arange(Npix - 1)

        # PLOTS

        for tag in ['hor', 'ver']:
            profscov_1D.data[tag]['labelkeys'] = \
                list(profscov_1D.data[tag][CCDs[0]][Quads[0]]['x'].keys())

        for tag in ['ver', 'hor']:

            fdict_C = self.figdict['BF01_COV_%s' % tag][1]
            fdict_C['data'] = profscov_1D.data[tag].copy()
            if self.report is not None:
                self.addFigures_ST(figkeys=['BF01_COV_%s' % tag],
                                   dobuilddata=False)

        # Saving Profiles

        profscov_1D.savetopickle()
        self.dd.products['profscov_1D_name'] = profscov_1D.rootname

        # Table of Results

        COV_dddf = OrderedDict(COV=pd.DataFrame.from_dict(COV_dd))

        covtable_cdp = self.CDP_lib['COVTABLE']
        covtable_cdp.path = self.inputs['subpaths']['products']
        covtable_cdp.ingest_inputs(
            data=COV_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy()
        )

        covtable_cdp.init_wb_and_fillAll(header_title='BF01: COVTABLE')
        self.save_CDP(covtable_cdp)
        self.pack_CDP_to_dd(covtable_cdp, 'COVTABLE_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def fcol(x): return 'col%03i' % x

            def fE(x): return '%.2E' % x

            cov_formatters = [fccd, fq, fcol] + [fE] * 6

            COVtex = covtable_cdp.get_textable(sheet='COV', caption='BF01: COV',
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=cov_formatters)

            self.report.add_Text(COVtex)

        return None

    def extract_BF(self):
        """

        Performs basic analysis of images:
            - extracts BF matrix for each COV matrix

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extractBF', Title='BF-Matrices Extraction', level=0)

        # labels should be the same accross CCDs. PATCH.
        label = self.dd.mx['label'][:, 0].copy()

        indices = copy.deepcopy(self.dd.indices)
        nObs, nC, nQ = indices.shape
        CCDs = np.array(indices.get_vals('CCD'))
        Quads = np.array(indices.get_vals('Quad'))

        ulabels = np.unique(label)
        nL = len(ulabels)

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        # INITIALISATIONS

        figspath = self.inputs['subpaths']['figs']
        kerpath = self.inputs['subpaths']['kernels']

        NP = nC * nQ * nL

        self.dd.products['BF'] = OrderedDict()
        self.dd.products['BF']['ulabels'] = ulabels.copy()
        self.dd.products['BF']['CCDs'] = CCDs.copy()
        self.dd.products['BF']['Quads'] = Quads.copy()

        BF_dd = OrderedDict()
        BF_dd['CCD'] = np.zeros(NP, dtype='int32')
        BF_dd['Q'] = np.zeros(NP, dtype='int32')
        BF_dd['col'] = np.zeros(NP, dtype='S6')
        BF_dd['FWHMx'] = np.zeros(NP, dtype='float32') + np.nan
        BF_dd['FWHMy'] = np.zeros(NP, dtype='float32') + np.nan
        BF_dd['e'] = np.zeros(NP, dtype='float32') + np.nan
        BF_dd['fluence'] = np.zeros(NP, dtype='float32') + np.nan

        profsker_1D = self.CDP_lib['PROFSKER1D']
        profsker_1D.header = CDP_header.copy()
        profsker_1D.path = kerpath
        profsker_1D.data = OrderedDict()

        profsker_1D.data['hor'] = OrderedDict()
        profsker_1D.data['ver'] = OrderedDict()

        for tag in ['hor', 'ver']:
            profsker_1D.data[tag] = OrderedDict()
            for CCDk in CCDs:
                profsker_1D.data[tag][CCDk] = OrderedDict()
                for Q in Quads:
                    profsker_1D.data[tag][CCDk][Q] = OrderedDict()
                    profsker_1D.data[tag][CCDk][Q]['x'] = OrderedDict()
                    profsker_1D.data[tag][CCDk][Q]['y'] = OrderedDict()

        for jCCD, CCDk in enumerate(CCDs):

            self.dd.products['BF'][CCDk] = OrderedDict()

            for kQ, Q in enumerate(Quads):
                self.dd.products['BF'][CCDk][Q] = OrderedDict()

        Npix = 51
        Npixplot = 11

        if not self.drill:

            singlepixmap = np.zeros((Npix, Npix), dtype='float32') + 0.0
            singlepixmap[(Npix - 1) / 2, (Npix - 1) / 2] = 1.

            for jCCD, CCDk in enumerate(CCDs):

                for ix, ulabel in enumerate(ulabels):

                    COV_dict = self.dd.products['COV'][CCDk][ulabel].copy()

                    for kQ, Q in enumerate(Quads):

                        jj = jCCD * (nQ * nL) + ix * nQ + kQ

                        COV_mx = COV_dict['av_covmap'][Q].copy()

                        BF_dd['CCD'][jj] = jCCD + 1
                        BF_dd['Q'][jj] = kQ + 1
                        BF_dd['col'][jj] = ulabel
                        BF_dd['fluence'][jj] = COV_dict['av_mu'][Q].copy()

                        try:

                            Asol_Q, psmooth_Q = G15.solve_for_A_linalg(
                                COV_mx, var=1., mu=1., returnAll=True, doplot=False,
                                verbose=False)

                            kernel_Q = G15.degrade_estatic(singlepixmap, Asol_Q)

                            cross_Q = kernel_Q[Npix / 2 - 1:Npix / 2 +
                                               2, Npix / 2 - 1:Npix / 2 + 2].copy()
                            kerQshape = G15.get_cross_shape_rough(
                                cross_Q, pitch=12.)

                            self.dd.products['BF'][CCDk][Q][ulabel] = OrderedDict(
                                Asol=Asol_Q.copy(),
                                psmooth=copy.deepcopy(psmooth_Q),
                                kernel=kernel_Q.copy(),
                                fluence=BF_dd['fluence'][jj],
                                fwhmx=kerQshape['fwhmx'],
                                fwhmy=kerQshape['fwhmy'],
                                ell=kerQshape['e'])

                            BF_dd['FWHMx'][jj] = kerQshape['fwhmx']
                            BF_dd['FWHMy'][jj] = kerQshape['fwhmy']
                            BF_dd['e'][jj] = kerQshape['e']

                            profsker_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                                np.arange(Npixplot) - Npixplot / 2
                            profsker_1D.data['hor'][CCDk][Q]['y'][ulabel] = kernel_Q[Npix / \
                                2, Npix / 2 - Npixplot / 2:Npix / 2 + Npixplot / 2 + 1].copy()

                            profsker_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                                np.arange(Npixplot) - Npixplot / 2

                            profsker_1D.data['ver'][CCDk][Q]['y'][ulabel] = kernel_Q[Npix / \
                                2 - Npixplot / 2:Npix / 2 + Npixplot / 2 + 1, Npix / 2].copy()

                            # BEWARE, PENDING: dispfig is saved but NOT REPORTED anywhere!

                            dispfig = os.path.join(
                                figspath, 'DISTORT_BF01_%s_%s%s.png' %
                                (ulabel, CCDk, Q))

                            G15.show_disps_CCD273(Asol_Q, stretch=10., peak=BF_dd['fluence'][jj],
                                                  N=13, sigma=1.6,
                                                  title='%s:%s%s' % (ulabel, CCDk, Q),
                                                  figname=dispfig)

                        except BaseException:

                            self.dd.products['BF'][CCDk][Q][ulabel] = OrderedDict()

                            profsker_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                                np.arange(Npixplot) - Npixplot / 2
                            profsker_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                                np.zeros(Npixplot)

                            profsker_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                                np.arange(Npixplot) - Npixplot / 2

                            profsker_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                                np.zeros(Npixplot)
        else:

            for jCCD, CCDk in enumerate(CCDs):

                for ku, ulabel in enumerate(ulabels):

                    for lQ, Q in enumerate(Quads):

                        profsker_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                            np.arange(Npixplot) - Npixplot / 2
                        profsker_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                            np.zeros(Npixplot)

                        profsker_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                            np.arange(Npixplot) - Npixplot / 2
                        profsker_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                            np.zeros(Npixplot)

        # Plots

        for tag in ['hor', 'ver']:
            profsker_1D.data[tag]['labelkeys'] = \
                list(profsker_1D.data[tag][CCDs[0]][Quads[0]]['x'].keys())

        for tag in ['ver', 'hor']:

            fdict_K = self.figdict['BF01_KER_%s' % tag][1]
            fdict_K['data'] = profsker_1D.data[tag].copy()

            if self.report is not None:
                self.addFigures_ST(figkeys=['BF01_KER_%s' % tag],
                                   dobuilddata=False)

        # REPORTING

        BF_dddf = OrderedDict(BF=pd.DataFrame.from_dict(BF_dd))

        bftable_cdp = self.CDP_lib['BFTABLE']
        bftable_cdp.path = self.inputs['subpaths']['products']
        bftable_cdp.ingest_inputs(
            data=BF_dddf.copy(),
            meta=dict(ulabels=ulabels.copy(),
                      CCDs=CCDs.copy(),
                      Quads=Quads.copy()),
            header=CDP_header.copy()
        )

        bftable_cdp.init_wb_and_fillAll(header_title='BF01: G15 TABLE')
        self.save_CDP(bftable_cdp)
        self.pack_CDP_to_dd(bftable_cdp, 'BFTABLE_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def fcol(x): return x

            def fE(x): return '%.2E' % x

            cov_formatters = [fccd, fq, fcol] + [fE] * 4

            BFtex = bftable_cdp.get_textable(sheet='BF', caption='BF01: G15 results',
                                             fitwidth=True,
                                             tiny=True,
                                             formatters=cov_formatters)

            self.report.add_Text(BFtex)

    def meta_analysis(self):
        """

        Analyzes the BF results across fluences.


        """

        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='Meta-Analysis', level=0)

        # INPUTS

        BFTABLE_CDP_pick = self.dd.products['BFTABLE_CDP']
        BFTABLE_CDP = cPickleRead(BFTABLE_CDP_pick)

        BF_df = BFTABLE_CDP['data']['BF']
        CCDv = BF_df['CCD'].as_matrix()
        Qv = BF_df['Q'].as_matrix()
        FWHMx = BF_df['FWHMx'].as_matrix()
        FWHMy = BF_df['FWHMy'].as_matrix()
        #ell = BF_df['e'].as_matrix()
        flu = BF_df['fluence'].as_matrix()

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        # INITIALISATIONS

        indices = copy.deepcopy(self.dd.indices)
        nObs, nC, nQ = indices.shape
        CCDs = np.array(indices.get_vals('CCD'))
        Quads = np.array(indices.get_vals('Quad'))

        NP = nC * nQ

        BFfit_dd = OrderedDict()
        BFfit_dd['CCD'] = np.zeros(NP, dtype='int32')
        BFfit_dd['Q'] = np.zeros(NP, dtype='int32')
        BFfit_dd['FWHMx_HWC'] = np.zeros(NP, dtype='float32') + np.nan
        BFfit_dd['FWHMx_Slope'] = np.zeros(NP, dtype='float32') + np.nan
        BFfit_dd['FWHMy_HWC'] = np.zeros(NP, dtype='float32') + np.nan
        BFfit_dd['FWHMy_Slope'] = np.zeros(NP, dtype='float32') + np.nan
        BFfit_dd['ELL_HWC'] = np.zeros(NP, dtype='float32') + np.nan

        plot_FWHM_dict = OrderedDict()

        for tag in ['fwhmx', 'fwhmy']:
            plot_FWHM_dict[tag] = OrderedDict(labelkeys=['data', 'fit'])
            for CCDk in CCDs:
                plot_FWHM_dict[tag][CCDk] = OrderedDict()
                for Q in Quads:
                    plot_FWHM_dict[tag][CCDk][Q] = OrderedDict()
                    plot_FWHM_dict[tag][CCDk][Q]['x'] = OrderedDict()
                    plot_FWHM_dict[tag][CCDk][Q]['y'] = OrderedDict()

        for iCCD, CCDk in enumerate(CCDs):

            for kQ, Q in enumerate(Quads):

                jj = iCCD * nQ + kQ

                ixsel = np.where((iCCD + 1 == CCDv) & (kQ + 1 == Qv))

                iflu = flu[ixsel]
                ifwhmx = FWHMx[ixsel]
                ifwhmy = FWHMy[ixsel]
                ixorder = iflu.argsort()
                iflu = iflu[ixorder]
                ifwhmx = ifwhmx[ixorder]
                ifwhmy = ifwhmy[ixorder]

                plot_FWHM_dict['fwhmx'][CCDk][Q]['x']['data'] = iflu.copy()
                plot_FWHM_dict['fwhmx'][CCDk][Q]['y']['data'] = ifwhmx.copy()

                # res = dict(xfit=xfit.copy(),yfit=yfit.copy(),intercept=intercept,
                # slope=slope*1.e4)

                resX = fit_BF(iflu, ifwhmx)

                plot_FWHM_dict['fwhmx'][CCDk][Q]['x']['fit'] = resX['xfit'].copy()
                plot_FWHM_dict['fwhmx'][CCDk][Q]['y']['fit'] = resX['yfit'].copy()

                plot_FWHM_dict['fwhmy'][CCDk][Q]['x']['data'] = iflu.copy()
                plot_FWHM_dict['fwhmy'][CCDk][Q]['y']['data'] = ifwhmy.copy()

                resY = fit_BF(iflu, ifwhmy)

                plot_FWHM_dict['fwhmy'][CCDk][Q]['x']['fit'] = resY['xfit'].copy()
                plot_FWHM_dict['fwhmy'][CCDk][Q]['y']['fit'] = resY['yfit'].copy()

                BFfit_dd['CCD'][jj] = iCCD + 1
                BFfit_dd['Q'][jj] = kQ + 1
                BFfit_dd['FWHMx_HWC'][jj] = resX['fwhm_hwc']
                BFfit_dd['FWHMx_Slope'][jj] = resX['slope']
                BFfit_dd['FWHMy_HWC'][jj] = resY['fwhm_hwc']
                BFfit_dd['FWHMy_Slope'][jj] = resY['slope']
                sigmax = resX['fwhm_hwc'] / 2.355
                sigmay = resY['fwhm_hwc'] / 2.355
                ell = np.abs((sigmax**2. - sigmay**2.) / (sigmax**2. + sigmay**2.))
                BFfit_dd['ELL_HWC'][jj] = ell

        for tag in ['fwhmx', 'fwhmy']:

            fdict_FF = self.figdict['BF01_%s_v_flu' % tag][1]
            fdict_FF['data'] = plot_FWHM_dict[tag].copy()

            if self.report is not None:
                self.addFigures_ST(figkeys=['BF01_%s_v_flu' % tag],
                                   dobuilddata=False)

        # REPORTING

        BFfit_dddf = OrderedDict(BFFIT=pd.DataFrame.from_dict(BFfit_dd))

        bfFITtable_cdp = self.CDP_lib['BFfitTABLE']
        bfFITtable_cdp.path = self.inputs['subpaths']['products']
        bfFITtable_cdp.ingest_inputs(
            data=BFfit_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy()
        )

        bfFITtable_cdp.init_wb_and_fillAll(header_title='BF01: FIT G15 TABLE')
        self.save_CDP(bfFITtable_cdp)
        self.pack_CDP_to_dd(bfFITtable_cdp, 'BFfitTABLE_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def ff(x): return '%.3f' % x

            cov_formatters = [fccd, fq] + [ff] * 5

            caption = 'BF01: FIT G15 results. FWHMx\_HWC: FWHM(x)[um] at half ADC range (32768 ADU); ' +\
                      'FWHMx\_slope: FWHM(x) vs. fluence slope in um / 10 kADU; ' +\
                      'FWHMy\_HWC: FWHM(y)[um] at half ADC range (32768 ADU); ' +\
                      'FWHMy\_slope: FWHM(y) vs. fluence slope in um / 10 kADU; ' +\
                      'ELL\_HWC: ellipticity at half ADC range (32768 ADU).'

            BFfittex = bfFITtable_cdp.get_textable(sheet='BFFIT',
                                                   caption=caption,
                                                   fitwidth=True,
                                                   tiny=True,
                                                   formatters=cov_formatters)

            self.report.add_Text(BFfittex)

        self.canbecleaned = True
