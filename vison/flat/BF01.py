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

import matplotlib.cm as cm

# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


def process_one_fluence_covmaps(q, dd, dpath, CCDs, jCCD, ku, ulabels, 
                                **kwargs):
    """ """

    labels = dd.mx['label'][:, 0].copy()

    Npix = kwargs['Npix']
    covfunc = kwargs['covfunc']
    doBiasCorr = kwargs['doBiasCorr']
    central = kwargs['central']
    vstart = kwargs['vstart']
    vend = kwargs['vend']
    clipsigma = kwargs['clipsigma']

    vfullinpath_adder = utils.get_path_decorator(dpath)

    #taskcounter = 0

    ulabel = ulabels[ku]

    six = np.where(labels == ulabel)

    ccdobjNamesList = vfullinpath_adder(dd.mx['ccdobj_name'][six[0], jCCD], 'pick')

    print(('%s: column %i/%i; %i OBSIDs' % (CCDs[jCCD], ku + 1, len(ulabels), len(ccdobjNamesList))))

    ccdobjList = [cPickleRead(item) for item in ccdobjNamesList]

    icovdict = covlib.get_cov_maps(
        ccdobjList, Npix=Npix, vstart=vstart, vend=vend, clipsigma=clipsigma,
        covfunc = covfunc, doBiasCorr=doBiasCorr,central=central,
        doTest=False, debug=False)
    q.put([jCCD, ku, icovdict])


def correct_BFE_one_image(q, dd, inputs, iObs, nObs, CCDs, Quads, picklespath, Asol):
    """ """
    if Asol is None:
        tempccdobj = '%s_proc_bfe'
    else:
        tempccdobj = '%s_proc_bfe_alt'

    for jCCD, CCDkey in enumerate(CCDs):

        inccdobj_name = '%s.pick' % dd.mx['ccdobj_name'][iObs, jCCD]
        full_inccdobj_name = os.path.join(picklespath, inccdobj_name)


        print(('Test %s, OBS %i/%i: correcting BFE in %s...' % (
            inputs['test'], iObs + 1, nObs, inccdobj_name)))

        # loading CCD Object
        ccdobj = copy.deepcopy(cPickleRead(full_inccdobj_name))

        ccdobj_bfe_name = tempccdobj % dd.mx['File_name'][iObs, jCCD]

        fullccdobj_bfe_name = os.path.join(
            picklespath, '%s.pick' % ccdobj_bfe_name)

        colkey = dd.mx['label'][iObs,jCCD]

        hasallAsols = True

        for Q in Quads:

            if Asol is None:
                try:
                    _Asol = dd.products['BF'][CCDkey][Q][colkey]['Asol'].copy()
                except:
                    _Asol = None
                    hasallAsols = False
            else:
                _Asol = Asol[CCDkey][Q].copy()

            if _Asol is not None:
                Qimg = ccdobj.get_quad(Q, canonical=True, extension=-1)
                Qimg = G15.correct_estatic(Qimg, _Asol)
                ccdobj.set_quad(Qimg, Q, canonical=True, extension=-1)

        if hasallAsols:
            cPickleDumpDictionary(ccdobj, fullccdobj_bfe_name)
            q.put([iObs, jCCD, ccdobj_bfe_name])
        else:
            q.put([iObs, jCCD, 'None.pick'])


class BF01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
        ('Npix', ([int], 'Number of Pixels (linear) to consider for Covariance Matrix')),
        ('clipsigma', ([float], 'Nsigma clipping.')),
        ('covfunc', ([str], 'ID of the covariance function.')),
        ('doBiasCorr', ([bool], 'Do correction bias for sigma clipping.')),
        ('central', ([str], 'Central Estimator.')),
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
                         ('correct_BFE_G15', self.correct_BFE_G15),
                         ('extract_PTCs', self.extract_PTCs),
                         ('meta', self.meta_analysis)]
                         #('debugtask', self.debugtask)]
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
        self.window = dict(wpx=300, hpx=300)

    def set_inpdefaults(self, **kwargs):
        """ """

        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_inpdefaults(**_kwargs)
        self.inpdefaults['test'] = kwargs['test']
        self.inpdefaults['surrogate'] = kwargs['surrogate']
        self.inpdefaults['Npix'] = 5
        self.inpdefaults['clipsigma'] = 4.
        self.inpdefaults['covfunc'] = 'ver2'
        self.inpdefaults['doBiasCorr'] = True
        self.inpdefaults['central'] = 'mean'

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
        clipsigma = self.inputs['clipsigma']
        covfunc = self.inputs['covfunc']
        doBiasCorr = self.inputs['doBiasCorr']
        central = self.inputs['central']

        # labels should be the same accross CCDs. PATCH.
        labels = self.dd.mx['label'][:, 0].copy()
        ulabels = np.unique(labels)
        #ulabels = ['col003'] # TEST
        nL = len(ulabels)

        # vstart and vend should be the same for all OBSIDs in test
        vstart = self.dd.mx['vstart'][0, 0]
        vend = min(self.dd.mx['vend'][0, 0], self.ccdcalc.NrowsCCD)

        indices = copy.deepcopy(self.dd.indices)
        nObs, nC, nQ = indices.shape[0:3]
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
            profscov_1D.data[tag]['labelkeys'] = ulabels
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
        COV_dd['CORR_00'] = np.zeros(NP, dtype='float32')
        COV_dd['CORR_01'] = np.zeros(NP, dtype='float32')
        COV_dd['CORR_10'] = np.zeros(NP, dtype='float32')
        COV_dd['CORR_11'] = np.zeros(NP, dtype='float32')

        self.dd.products['COV'] = OrderedDict()
        for CCDk in CCDs:
            self.dd.products['COV'][CCDk] = OrderedDict()

        if not self.drill:

            # doTest=False

            kwargs = dict(Npix=Npix, vstart=vstart, vend=vend, 
                clipsigma=clipsigma, covfunc=covfunc,
                doBiasCorr=doBiasCorr,central=central)

            arglist = []

            mgr = mp.Manager()
            queue = mgr.Queue()

            for jCCD, CCDk in enumerate(CCDs):

                for ku, ulabel in enumerate(ulabels):

                    arglist.append([queue, self.dd, dpath, CCDs, jCCD, ku, ulabels])

            #process_one_fluence_covmaps(*arglist[-3],**kwargs) # TEST
            #stop()

            pool = mp.Pool(processes=self.processes)

            for ii in range(len(arglist)):
                pool.apply_async(process_one_fluence_covmaps, args=arglist[ii],
                                    kwds=kwargs)
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
                    COV_dd['col'][jj] = int(ulabel.replace('col', ''))
                    COV_dd['av_mu'][jj] = icovdict['av_mu'][Q]
                    COV_dd['av_var'][jj] = icovdict['av_var'][Q]
                    COV_dd['CORR_00'][jj] = icovdict['av_corrmap'][Q][0, 0]
                    COV_dd['CORR_01'][jj] = icovdict['av_corrmap'][Q][0, 1]
                    COV_dd['CORR_10'][jj] = icovdict['av_corrmap'][Q][1, 0]
                    COV_dd['CORR_11'][jj] = icovdict['av_corrmap'][Q][1, 1]


                    profscov_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                        np.arange(Npix - 1)
                    profscov_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                        icovdict['av_corrmap'][Q][1:, 0].copy()

                    profscov_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                        np.arange(Npix - 1)
                    profscov_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                        icovdict['av_corrmap'][Q][0, 1:].copy()

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

        #stop() # TESTS

        # PLOTS

        #for tag in ['hor', 'ver']:
        #    profscov_1D.data[tag]['labelkeys'] = \
        #        profscov_1D.data[tag][CCDs[0]][Quads[0]]['x'].keys()

        rbcolors = cm.rainbow(np.linspace(0, 1, len(ulabels)))


        for tag in ['ver', 'hor']:

            fdict_C = self.figdict['BF01_COV_%s' % tag][1]
            fdict_C['data'] = profscov_1D.data[tag].copy()

            for ic, ulabel in enumerate(ulabels):
                fdict_C['meta']['corekwargs'][ulabel]=\
                    dict(marker='.', linestyle='-',color=rbcolors[ic])

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
        nObs, nC, nQ = indices.shape[0:3]
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
            profsker_1D.data[tag] = OrderedDict(labelkeys=ulabels)
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

                        CORR_mx = COV_dict['av_corrmap'][Q].copy()

                        BF_dd['CCD'][jj] = jCCD + 1
                        BF_dd['Q'][jj] = kQ + 1
                        BF_dd['col'][jj] = ulabel
                        fluence = COV_dict['av_mu'][Q]
                        BF_dd['fluence'][jj] = fluence

                        try:

                            Asol_Q, psmooth_Q = G15.solve_for_A_linalg(
                                CORR_mx, var=1., mu=fluence, returnAll=True, doplot=False,
                                verbose=False)


                            kernel_Q = G15.degrade_estatic(singlepixmap*fluence, Asol_Q)

                            cross_Q = kernel_Q[Npix / 2 - 1:Npix / 2 + 2, 
                                               Npix / 2 - 1:Npix / 2 + 2].copy()
                            kerQshape = G15.get_cross_shape_rough(
                                cross_Q, pitch=12.)


                            #kerQshapealt = BF01aux.get_kernel_gauss_shape(kernel_Q,pitch=12)


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
                            profsker_1D.data['hor'][CCDk][Q]['y'][ulabel] = np.log10(kernel_Q[Npix / \
                                2 - Npixplot / 2:Npix / 2 + Npixplot / 2 + 1, Npix / 2].copy())

                            profsker_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                                np.arange(Npixplot) - Npixplot / 2

                            profsker_1D.data['ver'][CCDk][Q]['y'][ulabel] = np.log10(kernel_Q[Npix / \
                                2, Npix / 2 - Npixplot / 2:Npix / 2 + Npixplot / 2 + 1].copy())


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

        #for tag in ['hor', 'ver']:
        #    profsker_1D.data[tag]['labelkeys'] = \
        #        profsker_1D.data[tag][CCDs[0]][Quads[0]]['x'].keys()

        rbcolors = cm.rainbow(np.linspace(0, 1, len(ulabels)))

        for tag in ['ver', 'hor']:

            fdict_K = self.figdict['BF01_KER_%s' % tag][1]
            fdict_K['data'] = profsker_1D.data[tag].copy()

            for ic, ulabel in enumerate(ulabels):
                fdict_K['meta']['corekwargs'][ulabel]=\
                    dict(marker='.', linestyle='-',color=rbcolors[ic])

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

    def f_correct_BFE_G15(self, ccdobjname, fixA=False):
        """Applies BFE solutions from G+15 to images, to later test effectivity 
        through PTC."""


        # Initialize new columns

        Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
        self.dd.initColumn(ccdobjname, Cindices,
                           dtype='S100', valini='None')

        DDindices = copy.deepcopy(self.dd.indices)

        #nObs,nCCD,nQuad = DDindices.shape
        #Quads = DDindices[2].vals

        nObs = DDindices.get_len('ix')
        # nObs = 3  # TESTS!
        #print 'TESTS: task.prepare_images: LIMITTING TO 3 IMAGES!'

        CCDs = DDindices.get_vals('CCD')
        Quads = DDindices.get_vals('Quad')


        if not self.drill:

            picklespath = self.inputs['subpaths']['ccdpickles']

            if fixA:

                fluences = self.dd.mx['flu_med_img'].array.mean(axis=(1,2))
                minflu = np.nanmin(fluences)
                maxflu = np.nanmax(fluences)
                midflu = (minflu+maxflu)/2.
                ixsel = np.argmin(np.abs(fluences-midflu))
                colkey = self.dd.mx['label'][ixsel,0]


                Asol = OrderedDict()
                for CCDk in CCDs:
                    Asol[CCDk] = OrderedDict()
                    for Q in Quads:
                        Asol[CCDk][Q] = self.dd.products['BF'][CCDk][Q][colkey]['Asol'].copy()

            else: 
                Asol =None

            arglist = []

            mgr = mp.Manager()
            queue = mgr.Queue()

            for iObs in range(nObs):
                arglist.append([queue, self.dd, self.inputs,
                                iObs, nObs, CCDs, Quads, picklespath, Asol])

            #correct_BFE_one_image(*arglist[0]) # TEST
            #nobfe = 'results_atCALDATA/DTEST/BF01_730/ccdpickles/EUC_31074_300719D121603T_ROE1_CCD1_proc_bfe.pick'
            #wbfe = 'results_atCALDATA/DTEST/BF01_730/ccdpickles/EUC_31074_300719D121603T_ROE1_CCD1_proc.pick'
            #ccdobj_nobfe = cPickleRead(nobfe)
            #ccdobj_wbfe = cPickleRead(wbfe)
            #arglist = [arglist[0]]


            pool = mp.Pool(processes=self.processes)

            for i in range(len(arglist)):
                pool.apply_async(correct_BFE_one_image, args=arglist[i])
            pool.close()
            pool.join()

            replies = []
            while not queue.empty():
                replies.append(queue.get())

            for reply in replies:
                iObs, jCCD, ccdobj_name = reply
                self.dd.mx[ccdobjname][iObs, jCCD] = ccdobj_name

        return None

    def correct_BFE_G15(self):
        """ """

        if self.report is not None:
            self.report.add_Section(
                keyword='correct_BFE_G15', Title='Correcting BFE (G+15)', level=0)

        self.f_correct_BFE_G15('ccdobj_bfe_fixA_name', fixA=True)
        self.f_correct_BFE_G15('ccdobj_bfe_name', fixA=False)


    def extract_PTCs(self):
        """ """
        # HARDWIRED VALUES
        wpx = self.window['wpx']
        hpx = self.window['hpx']

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='PTC Extraction', level=0)
            self.report.add_Text('Segmenting on %i x %i windows...' % (wpx, hpx))


        medcol = 'sec_med'
        varcol = 'sec_var'
        ccdobjcol = 'ccdobj_name'

        self.f_extract_PTC(ccdobjcol, medcol, varcol)


        medcol_nobfe = 'sec_med_noBFE'
        varcol_nobfe = 'sec_var_noBFE'
        ccdobjcol_nobfe = 'ccdobj_bfe_name'


        self.f_extract_PTC(ccdobjcol_nobfe, medcol_nobfe, varcol_nobfe)


        medcol_nobfealt = 'sec_med_noBFEalt'
        varcol_nobfealt = 'sec_var_noBFEalt'
        ccdobjcol_nobfealt = 'ccdobj_bfe_fixA_name'


        self.f_extract_PTC(ccdobjcol_nobfealt, medcol_nobfealt, varcol_nobfealt)



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
        nObs, nC, nQ = indices.shape[0:3]
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

        hasPTCs = ('sec_var' in self.dd.colnames) and\
         ('sec_var_noBFE' in self.dd.colnames)

        if hasPTCs:

            av_gain = 3.5 # just for display purposes

            has_gain_cal = False
            #try:
            #    wave = self.inputs['wavelength']
            #    gaincdp = self.inputs['inCDPs']['Gain']['nm%i' % wave]
            #    df = cPickleRead(gaincdp)['data']['GAIN_TB']

            #    gaindict = dict()
            #    for iCCD, CCDk in enumerate(CCDs):
            #        gaindict[CCDk] = dict()
            #        for jQ, Q in enumerate(Quads):
            #            _gain = df.loc[df['CCD']==iCCD+1].loc[df['Q']==jQ+1]['gain'].values[0]
            #            gaindict[CCDk][Q] = _gain

            #except:
            #    has_gain_cal = False
            #    if self.log is not None:
            #        self.log.info('Gain matrix not found!')

            medcols = dict(BFE='sec_med',
                            NOBFE='sec_med_noBFE',
                            NOBFEALT='sec_med_noBFEalt')
            varcols = dict(BFE='sec_var',
                            NOBFE='sec_var_noBFE',
                            NOBFEALT='sec_var_noBFEalt')

            labelkeysPTC = ['data','theo']
            curves_cdp = OrderedDict(BFE=OrderedDict(labelkeys=labelkeysPTC),
                NOBFE=OrderedDict(labelkeys=labelkeysPTC),
                NOBFEALT=OrderedDict(labelkeys=labelkeysPTC))

            PTCkeys = ['BFE','NOBFE','NOBFEALT']

            for key in PTCkeys:
                for CCDk in CCDs:
                    curves_cdp[key][CCDk] = OrderedDict()
                    for Q in Quads:
                        curves_cdp[key][CCDk][Q] = OrderedDict()
                        curves_cdp[key][CCDk][Q]['x'] = OrderedDict()
                        curves_cdp[key][CCDk][Q]['y'] = OrderedDict()

            for iCCD, CCDk in enumerate(CCDs):
                for jQ, Q in enumerate(Quads):
                    ixsel = np.where(~np.isnan(self.dd.mx['ObsID_pair'][:]))

                    for key in PTCkeys:
                        medcol = medcols[key]
                        varcol = varcols[key]

                        raw_var = self.dd.mx[varcol][ixsel, iCCD, jQ, :]
                        raw_med = self.dd.mx[medcol][ixsel, iCCD, jQ, :]
                        ixnonan = np.where(~np.isnan(raw_var) & ~np.isnan(raw_med))
                        var = raw_var[ixnonan]
                        med = raw_med[ixnonan]

                        curves_cdp[key][CCDk][Q]['x']['data'] = med.copy()
                        curves_cdp[key][CCDk][Q]['y']['data'] = var.copy()

                        if has_gain_cal:
                            _gain = gaindict[CCDk][Q]
                        else:
                            _gain = av_gain


                        curves_cdp[key][CCDk][Q]['x']['theo'] = med.copy()
                        curves_cdp[key][CCDk][Q]['y']['theo'] = med.copy()/_gain

            for key in PTCkeys:

                fdict_PTC = self.figdict['BF01_PTC_%s' % key][1]
                fdict_PTC['data'] = curves_cdp[key].copy()

            if self.report is not None:
                self.addFigures_ST(figkeys=['BF01_PTC_BFE',
                    'BF01_PTC_NOBFE',
                    'BF01_PTC_NOBFEALT'],
                    dobuilddata=False)

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
