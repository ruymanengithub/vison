#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Generic FPA Test Class

Created on Tue Aug 20 11:25:55 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import os
import datetime
import sys
import string as st
import traceback
import multiprocessing as mp
import numpy as np
import copy
from skimage import exposure
import pandas as pd

from vison.support import vistime
from vison import __version__
from vison.support import files
from vison.support.report import Report
from vison.datamodel import fpa_dm
from vison.datamodel import inputs
from vison.pipe import task
from vison.pipe import lib as pilib
from vison.fpatests import fpatask_lib as fpatasklib
from vison.fpa import fpa as fpamod
from vison.metatests.metacal import MetaCal
from vison.datamodel import cdp as cdpmod

# from vison.metatests
# END IMPORT


class FpaTask(task.Task):
    """ """

    inputsclass = inputs.Inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.ID = None
        if 'ID' in inputs:
            self.ID = inputs['ID']
        self.processes = 1
        if 'processes' in inputs:
            self.processes = inputs['processes']
        if 'elvis' in inputs:
            self.elvis = inputs['elvis']
        else:
            self.elvis = 'FPA'

        self.Model = 'FPA'
        self.internals = dict()
        self.inputs = self.inputsclass()
        self.inpdefaults = dict()
        self.perfdefaults = dict()
        self.log = log
        self.report = None
        self.TestReference = '7-XXX'
        self.name = ''
        self.type = 'Fpa'
        self.HKKeys = []
        self.CDP_lib = dict()
        self.figdict = dict()
        if not hasattr(self, 'subtasks'):
            self.subtasks = [()]
        self.perflimits = dict()
        self.drill = drill
        self.debug = debug
        self.proc_histo = dict(Extract=False)
        self.cleanafter = cleanafter
        self.canbecleaned = False
        #self.subpaths2clean = ['ccdpickles','ccdflats']

        preprocessing = dict()
        preprocessing['offsetkwargs'] = dict(method='row',
                                             scan='pre', trimscan=[25, 5],
                                             ignore_pover=True,
                                             extension=-1)

        self.set_inpdefaults(**inputs)
        _inputs = self.inpdefaults.copy()

        if 'preprocessing' not in _inputs:
            _inputs['preprocessing'] = preprocessing.copy()

        _inputs['todo_flags'] = self.init_todo_flags()

        _inputs.update(inputs)

        self.inputs.update(_inputs)

        self.fpa = fpamod.FPA(self.inputs['FPAdesign'])
        self.NSLICES_FPA = self.fpa.NSLICES
        self.NCOLS_FPA = self.fpa.NCOLS
        self.CCDs = [1, 2, 3]
        self.Quads = ['E', 'F', 'G', 'H']

        self.metacal = MetaCal()

#        self.set_perfdefaults(**inputs)
#        _perfdefaults = self.perfdefaults.copy()
#        self.perflimits.update(_perfdefaults)
#        if 'perflimits' in self.inputs and self.inputs['perflimits'] is not None:
#            self.perflimits.update(self.inputs['perflimits'])
#
#        if 'diffvalues' in inputs and self.inputs['diffvalues'] is not None:
#            diffvalues = inputs['diffvalues'].copy()
#        else:
#            diffvalues = {}


#        self.inputs['structure'] = self.build_scriptdict(
#            diffvalues, elvis=self.elvis)
#
#        images_format = self.get_images_format()

#        NAXIS2withpover = 2*(ccd.NrowsCCD+ccd.voverscan)
#        emptyccd = ccd.CCD(withpover=images_format[1]==NAXIS2withpover)
#        self.ccdcalc = copy.deepcopy(emptyccd)

        self.CDP_header = OrderedDict()

    def set_inpdefaults(self, **kwargs):
        pass

    def init_todo_flags(self):
        init_todo_flags = dict(init=True, check=False, report=False)
        if len(self.subtasks[0]) > 0:
            for v in self.subtasks:
                init_todo_flags[v[0]] = False
        return init_todo_flags

    def __call__(self):
        """Generic test master function."""

        Errors = False

        self.CDP_header = OrderedDict(ID=self.ID,
                                      vison=__version__)

        # INPUTS

        subtasks = self.subtasks

        # inputs loading
        resultspath = self.inputs['resultspath']
        try:
            _paths = self.inputs['subpaths']
        except KeyError:
            _paths = dict()
        testkey = self.inputs['test']
        todo_flags = self.inputs['todo_flags']
        try:
            reportroot = self.inputs['reportroot']
        except KeyError:
            reportroot = '%s_report' % testkey
        try:
            cleantexafter = self.inputs['cleantexafter']
        except KeyError:
            cleantexafter = False

        DataDictFile = os.path.join(resultspath, '%s_DataDict.pick' % testkey)
        reportobjFile = os.path.join(resultspath, '%s_Report.pick' % testkey)

        for _pathkey in _paths:
            subpath = os.path.join(resultspath, _paths[_pathkey])
            _paths[_pathkey] = subpath
        self.inputs['subpaths'] = _paths

        if todo_flags['init']:

            if self.log is not None:
                self.log.info('Initializing: %s' %
                              (self.__module__,))

            if os.path.exists(DataDictFile):
                os.system('rm %s' % DataDictFile)
            if os.path.exists(reportobjFile):
                os.system('rm %s' % reportobjFile)

            # Creating/clearing resultspath
            if not os.path.exists(resultspath):
                os.system('mkdir %s' % resultspath)
            else:
                os.system(
                    'find %s -maxdepth 1 -type f -exec rm -f {} \;' % resultspath)

            # Creating/clearing sub-resultspath
            for _, subpath in self.inputs['subpaths'].items():
                if not os.path.exists(subpath):
                    os.system('mkdir %s' % subpath)
                else:
                    os.system('rm -rf %s/*' % subpath)

            # Initialising Report Object

            if todo_flags['report']:

                self.report = Report(TestName=testkey, Model=self.Model,
                                     Reference=self.TestReference)
                self.report.add_Section(
                    keyword='init', Title='Inputs \& Data Ingestion', level=0)
                self.add_inputs_to_report()
            else:
                self.report = None

            self.ingest_data()

            self.save_progress(DataDictFile, reportobjFile)
        else:
            self.recover_progress(DataDictFile, reportobjFile)

        # DATA-WORK and ANALYSIS

        for subtask in subtasks:

            subtaskname, subtaskmethod = subtask

            if subtaskname not in todo_flags:
                todo_flags[subtaskname] = False

            if todo_flags[subtaskname]:

                if self.log is not None:
                    self.log.info('Executing %s: %s' %
                                  (subtaskname, subtaskmethod.__module__))

                tini = datetime.datetime.now()
                try:
                    subtaskmethod()
                    tend = datetime.datetime.now()
                    dtm = ((tend - tini).seconds) / 60.
                    if self.log is not None:
                        self.log.info(
                            '%.1f minutes in running Sub-task: %s' % (dtm, subtaskname))
                    #print 'saving progress!'
                    self.save_progress(DataDictFile, reportobjFile)
                except BaseException:
                    Errors = True
                    self.dd.flags.add('SUBTASKCRASH')
                    self.catchtraceback()
                    # self.save_progress(DataDictFile,reportobjFile)
                    if not self.debug:
                        if self.log is not None:
                            self.log.info('SUBTASK "%s:%s" FAILED, QUITTING!' % (
                                subtaskname, subtaskmethod.__name__))
                        break
                    else:
                        sys.exit()

                if self.cleanafter and self.canbecleaned:
                    self.cleanaux()
                    self.canbecleaned = False

            else:
                self.recover_progress(DataDictFile, reportobjFile)

        # Write automatic Report of Results

        if todo_flags['report']:
            outfiles = self.report.doreport(
                reportroot,
                cleanafter=cleantexafter,
                silent=True)  # commented on TESTS
            # outfiles = self.report.doreport(reportroot, cleanafter=False, silent=False) # TESTS
            # stop() # TESTS

            for outfile in outfiles:
                os.system('mv %s %s/' % (outfile, resultspath))

        self.save_progress(DataDictFile, reportobjFile)

        if self.log is not None:
            self.log.info('Finished %s' % self.name)

        return Errors

    def catchtraceback(self):
        """ """
        msg_trbk = traceback.format_exc()
        if self.log is not None and not self.debug:
            self.log.info(msg_trbk)
        else:
            print(msg_trbk)

    def ingest_data(self):
        """ """
        testkey = self.inputs['test']
        datapath = self.inputs['datapath']
        OBSID_lims = self.inputs['OBSID_lims']
        #structure = self.inputs['structure']
        explogf = self.inputs['explogf']

        explog = pilib.loadexplogs(explogf, elvis=self.elvis, addpedigree=True,
                                   datapath=datapath)

        explog, checkreport = self.filterexposures(explog, OBSID_lims)

        if self.log is not None:
            self.log.info('%s acquisition consistent with expectations: %s' % (
                testkey, checkreport['checksout']))
            if len(checkreport['failedcols']) > 0:
                self.log.info('%s failed columns: %s' %
                              (testkey, checkreport['failedcols']))
            if len(checkreport['failedkeys']) > 0:
                self.log.info('%s failed keys: %s' %
                              (testkey, checkreport['failedkeys']))
            if len(checkreport['msgs']) > 0:
                self.log.info(['_'] + checkreport['msgs'])

        if self.report is not None:
            ntestkey = st.replace(testkey, '_', '\_')
            nchecksout = ['\\bf{%s}' % checkreport['checksout']]
            nchecksout = [st.replace(
                item, 'False', '$\\textcolor{red}{\\bf{False}}$') for item in nchecksout][0]
            self.report.add_Text(
                '%s acquisition consistent with expectations: %s\\newline' % (ntestkey, nchecksout))

            if (checkreport['failedcols']) > 0:
                nfailedcols = st.replace(
                    checkreport['failedcols'].__repr__(), '_', '\_')
                self.report.add_Text('%s failed columns: %s' %
                                     (ntestkey, nfailedcols))
            if len(checkreport['failedkeys']) > 0:
                nfailedkeys = st.replace(
                    checkreport['failedkeys'].__repr__(), '_', '\_')
                self.report.add_Text('%s failed keys: %s' %
                                     (ntestkey, nfailedkeys))
            if len(checkreport['msgs']) > 0:
                for msg in checkreport['msgs']:
                    nmsg = st.replace(msg, '_', '\_')
                    self.report.add_Text(nmsg)

        # Adding Time Axis

        explog['time'] = np.array(
            list(map(vistime.get_dtobj, explog['date']))).copy()

        # Building DataDict

        self.dd = fpatasklib.DataDict_builder(explog, self.inputs)

        if not checkreport['checksout']:
            self.dd.flags.add('MISSDATA')

    def load_LE1(self, LE1fits):
        """ """
        return fpa_dm.FPA_LE1(LE1fits)

    def check_data(self):
        """ """
        if self.report is not None:
            self.report.add_Section(
                keyword='check_data', Title='Data Validation', level=0)

    def iterate_over_CCDs(self, LE1, method, **kwargs):
        """ """

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                CCDID = 'C_%i%i' % (jY + 1, iX + 1)
                _kwargs = dict(LE1=LE1, CCDID=CCDID)
                _kwargs.update(kwargs)

                method(self, **_kwargs)

    def iterate_over_CCDs_parallel(self, LE1, method, **kwargs):
        """DOES NOT WORK!!"""

        arglist = []

        mgr = mp.Manager()
        queue = mgr.Queue()

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                CCDID = 'C_%i%i' % (jY + 1, iX + 1)
                _kwargs = dict(queue=queue, LE1=LE1, CCDID=CCDID)
                _kwargs.update(kwargs)
                arglist.append(['retrieve', _kwargs, {}])

        pool = mp.Pool(processes=self.processes)

        for i in range(len(arglist)):
            pool.apply_async(method, args=arglist[i])
        pool.close()
        pool.join()

        replies = []
        while not queue.empty():
            replies.append(queue.get())

        for reply in replies:
            method('apply', {}, reply)

    def filterexposures(self, explog, OBSID_lims):
        """Loads a list of Exposure Logs and selects exposures from test 'test'.

        The datapath becomes another column in DataDict. This helps dealing
        with tests that run overnight and for which the input data is in several
        date-folders.


        """

        if len(OBSID_lims) == 0:
            OBSID_lims = [explog['ObsID'][0], explog['ObsID'][-1]]

        testkey = self.inputs['test']

        selbool = (explog['test'] == testkey) & \
            (explog['ObsID'] >= OBSID_lims[0]) & \
            (explog['ObsID'] <= OBSID_lims[1])

        explog = explog[np.where(selbool)]

        # Assess structure - BYPASSED

        checkreport = dict(checksout=True,
                           failedkeys=[],
                           failedcols=[],
                           msgs=[])

        return explog, checkreport

    def get_FPAMAP(self, inData, extractor):

        M = OrderedDict()

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                M[Ckey] = OrderedDict()

                for Q in self.Quads:
                    M[Ckey][Q] = extractor(inData, Ckey, Q)

        return M


#    def doPlot(self, figkey, **kwargs):
#        """ """
#
#        try:
#            figobj = copy.deepcopy(self.figdict[figkey][0]())
#        except:
#            print 'DEBUGGING IN FpaTask.doPlot...'
#            msg_trbk = traceback.format_exc()
#            self.log.info(msg_trbk)
#            self.log.info('%s, %s' %
#                          (figkey, type(self.figdict[figkey][0])))
#            raise RuntimeError
#
#        figobj.configure(**kwargs)
#
#        self.figdict[figkey][0] = copy.deepcopy(figobj)

    def plot_SimpleMAP(self, MAPdict, kwargs):
        self.metacal.plot_SimpleMAP(MAPdict, kwargs)

    #plot_SimpleMAP = classmethod(MetaCal.plot_SimpleMAP.__func__)

    def plot_XY(self, XYdict, kwargs):
        self.metacal.plot_XY(XYdict, kwargs)

    def plot_XYCCD(self, XYCCD, kwargs):
        self.metacal.plot_XYCCD(XYCCD, kwargs)

    def plot_XYMAP(self, XYMAP, kwargs):
        self.metacal.plot_XYMAP(XYMAP, kwargs)

    def plot_ImgFPA(self, BCdict, kwargs):
        self.metacal.plot_ImgFPA(BCdict, kwargs)

    def iter_overCCDs(self, data, assigner, RetDict=None):
        """ """
        if RetDict is None:
            RetDict = dict()

        for jY in range(self.NCOLS_FPA):
            for iX in range(self.NSLICES_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                RetDict = assigner(RetDict, data, Ckey)

        return RetDict

    def get_ImgDictfromLE1(self, LE1, doequalise=False):
        """ """

        def assigner(ImgDict, LE1, Ckey):
            """ """
            locator = self.fpa.FPA_MAP[Ckey]
            flip = locator[2]

            kccdobj = LE1.get_ccdobj(Ckey)
            img = kccdobj.extensions[-1].data.transpose().copy().astype('float32')

            if doequalise:
                img = exposure.equalize_hist(img, nbins=256).astype('float32')

            res = dict(img=self.fpa.flip_img(img, flip))

            ImgDict[Ckey] = res.copy()

            return ImgDict

        ImgDict = self.iter_overCCDs(LE1, assigner)

        return ImgDict

    def add_StandardQuadsTable(self, extractor, cdp=None, cdpdict=None):
        """ """

        if cdpdict is not None:

            assert isinstance(cdpdict, dict)

            TBkey = cdpdict['TBkey']
            meta = cdpdict['meta'].copy()
            CDP_header = cdpdict['CDP_header'].copy()
            header_title = cdpdict['header_title']
            CDP_KEY = cdpdict['CDP_KEY']
            valformat = cdpdict['valformat']
            caption = cdpdict['caption']
        else:
            TBkey = 'TB'
            meta = dict()
            CDP_header = dict()
            header_title = 'generic title'
            CDP_KEY = 'UNK'
            caption = 'CAPTION PENDING'
            valformat = '%.2e'

        NCCDs = len(self.CCDs)
        NBlocks = len(self.fpa.flight_blocks)
        NP = NBlocks * NCCDs

        TB = OrderedDict()
        TB['CCDID'] = np.zeros(NP, dtype='S4')
        TB['BLOCK'] = np.zeros(NP, dtype='S4')
        TB['CCD'] = np.zeros(NP, dtype='S4')
        for Q in self.Quads:
            TB[Q] = np.zeros(NP, dtype='float32')

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]

                indx = jY * self.NCOLS_FPA + iX

                TB['CCDID'][indx] = Ckey
                TB['BLOCK'][indx] = block[0:4]
                TB['CCD'][indx] = CCDk

                for iQ, Q in enumerate(self.Quads):
                    TB[Q][indx] = extractor(self, Ckey, Q)

        TB_ddf = OrderedDict()
        TB_ddf[TBkey] = pd.DataFrame.from_dict(TB)

        if cdp is not None:
            cdp.path = self.inputs['subpaths']['products']
            cdp.ingest_inputs(
                data=TB_ddf.copy(),
                meta=meta.copy(),
                header=CDP_header.copy()
            )

            cdp.init_wb_and_fillAll(header_title=header_title)
            self.save_CDP(cdp)
            self.pack_CDP_to_dd(cdp, CDP_KEY)
        else:
            cdp = cdpmod.Tables_CDP()
            cdp.ingest_inputs(
                data=TB_ddf.copy(),
                meta=meta.copy(),
                header=CDP_header.copy()
            )

        if self.report is not None:

            def fstr(x): return '%s' % x

            def ff(x): return valformat % x

            her_formatters = [fstr, fstr, fstr]
            for iQ, Q in enumerate(self.Quads):
                her_formatters.append(ff)

            nicecaption = st.replace(caption, '_', '\_')
            Ttex = cdp.get_textable(sheet=TBkey, caption=nicecaption,
                                    fitwidth=True,
                                    tiny=True,
                                    formatters=her_formatters)

            self.report.add_Text(Ttex)

        return TB, cdp
