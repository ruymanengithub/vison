#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Generic Task (Test) Class.

Created on Tue Nov 14 14:20:04 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import datetime
import string as st
import copy
from collections import OrderedDict
import sys
import traceback
#from multiprocessing import Pool
import multiprocessing as mp

from vison.support.report import Report
from vison.support import vistime
from vison.support import files
from . import lib as pilib
from vison.support import context
#import task_lib as tlib
from vison.image import performance
from vison.datamodel import ccd
from vison.image import calibration
from vison.ogse import ogse as ogsemod
from vison.support.files import cPickleDumpDictionary
from vison.datamodel import compliance as complimod
from vison import __version__
# END IMPORT

isthere = os.path.exists

HKKeys = [
    'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
    'CCD1_TEMP_B', 'CCD2_TEMP_B', 'CCD3_TEMP_B',
    'CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T',
    'CCD1_OD_B', 'CCD2_OD_B', 'CCD3_OD_B',
    'COMM_RD_T', 'COMM_RD_B', 'VID_PCB_TEMP_T',
    'CCD1_IG1_T', 'CCD2_IG1_T', 'CCD3_IG1_T',
    'CCD1_IG1_B', 'CCD2_IG1_B', 'CCD3_IG1_B',
    'COMM_IG2_T', 'COMM_IG2_B', 'FPGA_BIAS_ID2',
    'VID_PCB_TEMP_T', 'VID_PCB_TEMP_B', 'RPSU_TEMP1',
    'FPGA_PCB_TEMP_T', 'FPGA_PCB_TEMP_B', 'RPSU_TEMP_2',
    'RPSU_28V_PRI_I']


def prepare_one_image(q, dd, ogse, inputs, iObs,
                      nObs, CCDs, picklespath, doBadPixels,
                      doMask, MaskData,
                      doOffset, offsetkwargs,
                      doBias, BiasData,
                      doFF, FFData):

    if doFF:
        FW_ID = dd.mx['wave'][iObs, 0]
        wavelength = ogse.get_wavelength(FW_ID)

    for jCCD, CCDkey in enumerate(CCDs):

        #vstart = self.dd.mx['vstart'][iObs, jCCD]
        #vend = self.dd.mx['vend'][iObs, jCCD]

        ccdobj_name = '%s_proc' % dd.mx['File_name'][iObs, jCCD]

        dpath = dd.mx['datapath'][iObs, jCCD]
        infits = os.path.join(dpath, '%s.fits' %
                              dd.mx['File_name'][iObs, jCCD])

        print('Test %s, OBS %i/%i: preparing %s...' % (
            inputs['test'], iObs + 1, nObs, infits))

        # converting the FITS file into a CCD Object
        ccdobj = ccd.CCD(infits)

        fullccdobj_name = os.path.join(
            picklespath, '%s.pick' % ccdobj_name)

        if doBadPixels:
            imgdata = ccdobj.extensions[-1].data.copy()
            BPmask = np.isclose(imgdata, 2**16 - 1.) | np.isclose(imgdata, 0.0)
            ccdobj.get_mask(BPmask)

        if doMask:
            ccdobj.get_mask(MaskData[CCDkey].extensions[-1].data)

        if doOffset:
            for Quad in ccdobj.Quads:
                # ccdobj.sub_offset(Quad,method='median',scan='pre',trimscan=[5,5],
                #                  ignore_pover=False)
                ccdobj.sub_offset(Quad, **offsetkwargs)

        if doBias:
            ccdobj.sub_bias(BiasData[CCDkey].extensions[-1].data,
                            extension=-1)

        if doFF:
            nmkey = 'nm%i' % wavelength
            if nmkey in FFData:
                FF=FFData[nmkey][CCDkey]
            else:
                FF = FFData[CCDkey]
            ccdobj.divide_by_flatfield(FF.extensions[FF.extnames.index('FLAT')].data,
                                       extension=-1)
            

        # cPickleDumpDictionary(dict(ccdobj=ccdobj),fullccdobj_name)
        cPickleDumpDictionary(ccdobj, fullccdobj_name)
        # ccdobj.writeto(fullccdobj_name,clobber=True)
        # self.dd.mx['ccdobj_name'][iObs, jCCD] = ccdobj_name
        q.put([iObs, jCCD, ccdobj_name])


# def _prepare_one_image(args):
#    prepare_one_image(*args)

class Task(object):
    """ """

    from .task_lib import check_HK, filterexposures, addHKPlotsMatrix, add_labels_to_explog
    from .task_lib import save_CDP, create_mockexplog, get_checkstats_T, check_metrics_T

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.ID = None
        if 'ID' in inputs:
            self.ID = inputs['ID']
        self.BLOCKID = None
        if 'BLOCKID' in inputs:
            self.BLOCKID = inputs['BLOCKID']
        self.CHAMBER = None
        if 'CHAMBER' in inputs:
            self.CHAMBER = inputs['CHAMBER']
        self.processes = 1
        if 'processes' in inputs:
            self.processes = inputs['processes']
        if 'elvis' in inputs:
            self.elvis = inputs['elvis']
        else:
            self.elvis = context.elvis

        self.ogse = ogsemod.Ogse(self.CHAMBER, withpover=True)

        self.Model = 'XM'
        self.internals = dict()
        self.inputs = self.inputsclass()
        self.inpdefaults = dict()
        self.perfdefaults = dict()
        self.log = log
        self.report = None
        self.TestReference = '7-XXX'
        self.name = ''
        self.type = 'Simple'
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
        self.subpaths2clean = ['ccdpickles', 'ccdflats']

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

        self.set_perfdefaults(**inputs)
        _perfdefaults = self.perfdefaults.copy()
        self.perflimits.update(_perfdefaults)
        if 'perflimits' in self.inputs and self.inputs['perflimits'] is not None:
            self.perflimits.update(self.inputs['perflimits'])

        if 'diffvalues' in inputs and self.inputs['diffvalues'] is not None:
            diffvalues = inputs['diffvalues'].copy()
        else:
            diffvalues = {}

        self.inputs['structure'] = self.build_scriptdict(
            diffvalues, elvis=self.elvis)

        images_format = self.get_images_format()

        NAXIS2withpover = 2 * (ccd.NrowsCCD + ccd.voverscan)
        emptyccd = ccd.CCD(withpover=images_format[1] == NAXIS2withpover)
        self.ccdcalc = copy.deepcopy(emptyccd)

        self.CDP_header = OrderedDict()

    def init_todo_flags(self):
        init_todo_flags = dict(init=True, check=False, report=False)
        if len(self.subtasks[0]) > 0:
            for v in self.subtasks:
                init_todo_flags[v[0]] = False
        return init_todo_flags

    def set_inpdefaults(self, **kwargs):
        pass

    def set_perfdefaults(self, **kwargs):
        self.perfdefaults = OrderedDict()
        self.perfdefaults.update(performance.get_perf_rdout(self.BLOCKID))
        self.perfdefaults['SATUR'] = OrderedDict()
        for CCD in ['CCD1', 'CCD2', 'CCD3']:
            self.perfdefaults['SATUR'][CCD] = [0., 0.1]
        # self.perfdefaults.update(kwargs)

    def build_scriptdict(self, diffvalues={}, elvis=context.elvis):
        """ """
        return dict()

    def get_images_format(self):

        strdict = self.inputs['structure']
        Ncols = strdict['Ncols']

        vstarts = []
        vends = []
        for i in range(1, Ncols + 1):
            vstarts.append(strdict['col%03i' % i]['vstart'])
            vends.append(strdict['col%03i' % i]['vend'])

        vstarts = np.array(vstarts)
        vends = np.array(vends)

        assert np.all(vstarts == vstarts[0])
        assert np.all(vends == vends[0])

        Nlines = vends - vstarts

        if Nlines[0] <= ccd.NrowsCCD:
            images_format = (ccd.NAXIS1,
                             ccd.NrowsCCD * 2)
        elif Nlines[0] == ccd.NrowsCCD + ccd.voverscan:
            images_format = (ccd.NAXIS1,
                             (ccd.NrowsCCD + ccd.voverscan) * 2)
        else:
            raise RuntimeError

        return images_format

    def __call__(self):
        """Generic test master function."""

        Errors = False

        self.CDP_header = OrderedDict(ID=self.ID,
                                      BLOCKID=self.BLOCKID,
                                      CHAMBER=self.CHAMBER,
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
            if not isthere(resultspath):
                os.system('mkdir %s' % resultspath)
            else:
                os.system(
                    'find %s -maxdepth 1 -type f -exec rm -f {} \;' % resultspath)

            # Creating/clearing sub-resultspath
            for _, subpath in self.inputs['subpaths'].items():
                if not isthere(subpath):
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

            if self.type == 'Simple':
                self.ingest_data_SimpleTest()
            elif self.type == 'Meta':
                self.ingest_data_MetaTest()

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

    def addHK_2_dd(self):
        """ """
        self.dd = pilib.addHK(self.dd, self.HKKeys, elvis=self.elvis)

    def addmockHK_2_dd(self):
        self.dd = pilib.addmockHK(self.dd, self.HKKeys, elvis=self.elvis)

    def ingest_data_SimpleTest(self):

        testkey = self.inputs['test']
        datapath = self.inputs['datapath']
        OBSID_lims = self.inputs['OBSID_lims']
        structure = self.inputs['structure']
        explogf = self.inputs['explogf']

        #elvis = self.inputs['elvis']

        if self.drill:
            if len(OBSID_lims) > 0:
                OBSID0 = OBSID_lims[0]
            else:
                OBSID0 = 1000
            explog = self.create_mockexplog(OBSID0=OBSID0)
        else:
            explog = pilib.loadexplogs(explogf, elvis=self.elvis, addpedigree=True,
                                       datapath=datapath)

        # META-DATA WORK
        explog, checkreport = self.filterexposures(
            structure, explog, OBSID_lims)

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

        self.dd = pilib.DataDict_builder(explog, self.inputs, structure)

        if not checkreport['checksout']:
            self.dd.flags.add('MISSDATA')

        # Add HK information
        if not self.drill:
            self.addHK_2_dd()
        else:
            self.addmockHK_2_dd()

    def ingest_data_MetaTest(self):
        raise NotImplementedError("Method implemented in child-class")

    def save_progress(self, DataDictFile, reportobjFile):
        """ """
        files.cPickleDumpDictionary(self.dd, DataDictFile)
        files.cPickleDump(self.report, reportobjFile)
        csvFile = st.replace(DataDictFile, '.pick', '.csv')
        self.dd.saveToFile(csvFile)

    def recover_progress(self, DataDictFile, reportobjFile):
        """ """
        self.dd = files.cPickleRead(DataDictFile)
        self.report = files.cPickleRead(reportobjFile)

    def cleanaux(self):
        """ """

        if not self.canbecleaned:
            return

        for subpathkey in self.subpaths2clean:
            if subpathkey in self.inputs['subpaths']:
                subpath = self.inputs['subpaths'][subpathkey]
                execline1 = "find %s/ -type f -name '*.fits' -exec sh -c '%s' {} \;" % (
                    subpath, 'rm "$0"')
                os.system(execline1)
                execline2 = "find %s/ -type f -name '*.pick' -exec sh -c '%s' {} \;" % (
                    subpath, 'rm "$0"')
                os.system(execline2)
                if self.log is not None:
                    self.log.info('\nCleared contents [.fits/.pick] of %s!' % subpath)

    def addFigure2Report(self, figkey):
        """ """
        figobj = self.figdict[figkey][0]
        figname = figobj.figname
        texfraction = figobj.texfraction
        caption = figobj.caption
        assert os.path.exists(figname)
        epsname = '%s.eps' % os.path.splitext(figname)[0]
        os.system('convert %s %s' % (figname, epsname))
        self.report.add_Figure(epsname, texfraction=texfraction,
                               caption=caption, label=figkey)

    def doPlot(self, figkey, **kwargs):
        """ """

        try:

            figobj = copy.deepcopy(self.figdict[figkey][0]())
        except BaseException:
            print('DEBUGGING IN Task.doPlot...')
            msg_trbk = traceback.format_exc()
            self.log.info(msg_trbk)
            self.log.info('%s, %s' %
                          (figkey, type(self.figdict[figkey][0])))
            raise RuntimeError

        figobj.configure(**kwargs)

        if kwargs['dobuilddata']:
            figobj.build_data(self.dd)
            # except:
            # stop()
        else:
            figobj.data = copy.deepcopy(kwargs['data'])
        if 'meta' in kwargs:
            meta = kwargs['meta']
        else:
            meta = {}

        figobj.plot(**meta)

        self.figdict[figkey][0] = copy.deepcopy(figobj)

    def addFigures_ST(self, dobuilddata=True, **kwargs):
        """ """
        try:
            figkeys = kwargs['figkeys']
        except BaseException:
            figkeys = []

        figspath = self.inputs['subpaths']['figs']

        for figkey in figkeys:
            try:
                pmeta = self.figdict[figkey][1]
                pmeta['path'] = figspath
                pmeta['dobuilddata'] = dobuilddata
                self.doPlot(figkey, **pmeta)
                self.addFigure2Report(figkey)
            except BaseException:
                self.catchtraceback()
                nfigkey = 'BS_%s' % figkey
                self.skipMissingPlot(nfigkey, ref=figkey)

    def addComplianceMatrix2Self(self, complidict, label):
        self.dd.compliances[label] = OrderedDict(complidict).copy()

    def addComplianceMatrix2Log(self, complidict, label=''):
        """ """
        #st_compl = complidict.__str__()
        st_compl = st.split(complidict.get_compliance_txt(), '\n')
        self.log.info([label, st_compl])

    def addComplianceMatrix2Report(self, complidict, label='',
                                   caption=''):
        """ """
        nicelabel = st.replace(label, ' ', '\ ')
        #st_compl = complidict.__str__()

        complitex = ['$\\bf{%s}$' % nicelabel]
        #complitex += complimod.gen_compliance_tex(complidict)
        complitex += complidict.get_compliance_tex(caption)
        # complitex = [st.replace(
        #    item, 'False', '$\\textcolor{red}{\\bf{False}}$') for item in complitex]

        self.report.add_Text(complitex)

    def IsComplianceMatrixOK(self, complidict):
        """ """

        def traverse_tree(dictionary, isOK):

            for key, value in list(dictionary.items()):
                #print 'Upper: %s' % key
                if isinstance(value, (dict, OrderedDict)):
                    #print key,value
                    isOK = isOK and traverse_tree(value, isOK)
                else:
                    #print key,value
                    isOK = isOK and value[0]
            return isOK

        isOK = traverse_tree(complidict, True)

        return isOK

    def addFlagsToLog(self):
        """ """
        flagstxt = st.join(self.dd.flags.getFlagsOnList(), ', ')
        self.log.info('FLAGS ON}:\n%s' % flagstxt)

    def addFlagsToReport(self):
        """ """
        niceflagnames = [st.replace(item, '_', '\_')
                         for item in self.dd.flags.getFlagsOnList()]
        flagstxt = st.join(niceflagnames, ', ')
        msgList = ['$\\bf{FLAGS\ ON}$: ', flagstxt]
        self.report.add_Text(msgList)

    def skipMissingPlot(self, key, ref):
        """ """
        self.figdict[key] = copy.deepcopy(self.figdict['BlueScreen'])
        niceref = st.replace(ref, '_', '\_')
        figspath = self.inputs['subpaths']['figs']
        pmeta = dict(path=figspath,
                     caption='$\\bf{MISSING}:$ %s' % niceref,
                     tag=ref,
                     meta=dict(title=niceref),
                     dobuilddata=True)
        self.doPlot(key, **pmeta)
        self.addFigure2Report(key)

    def check_stat_perCCD(self, arr, CCDlims, CCDs=['CCD1', 'CCD2', 'CCD3']):
        """ """
        compliance = complimod.ComplianceMX_CCD(
            CCDs=CCDs, CCDlims=CCDlims.copy())
        compliance.check_stat(arr)
        return compliance

    def check_stat_perCCDandQ(self, arr, CCDQlims, CCDs=['CCD1', 'CCD2', 'CCD3']):
        """ """
        compliance = complimod.ComplianceMX_CCDQ(
            CCDs=CCDs, CCDQlims=CCDQlims.copy())
        compliance.check_stat(arr)
        return compliance

    def check_stat_perCCDandCol(self, arr, lims, CCDs=['CCD1', 'CCD2', 'CCD3']):
        """ """
        colnames = list(lims[CCDs[0]].keys())
        compliance = complimod.ComplianceMX_CCDCol(colnames,
                                                   indexer=self.dd.mx['label'][:, 0],
                                                   CCDs=CCDs, lims=lims.copy())
        compliance.check_stat(arr)
        return compliance

    def check_stat_perCCDQandCol(self, arr, lims, CCDs=['CCD1', 'CCD2', 'CCD3']):
        """ """
        Qs = ['E', 'F', 'G', 'H']
        colnames = list(lims[CCDs[0]][Qs[0]].keys())

        compliance = complimod.ComplianceMX_CCDQCol(colnames,
                                                    indexer=self.dd.mx['label'][:, 0],
                                                    CCDs=CCDs,
                                                    Qs=Qs,
                                                    lims=lims.copy())
        compliance.check_stat(arr)
        return compliance

    def check_data(self, **kwargs):
        """Generic check_data method"""
        if self.report is not None:
            self.report.add_Section(
                keyword='check_data', Title='Data Validation', level=0)
        # INVENTORY OF DATA
        tDict = self.get_data_inventory_table()
        self.dd.meta['data_inventory'] = tDict.copy()
        if self.report is not None:
            self.add_data_inventory_to_report(tDict)
        # CHECK AND CROSS-CHECK HK
        self.check_HK_ST()
        # OBTAIN METRICS FROM IMAGES - TASK
        self.get_checkstats_T()
        # METRICS ASSESSMENT - TASK
        self.check_metrics_T()
        # OBTAIN METRICS FROM IMAGES - SUB-TASK
        self.get_checkstats_ST(**kwargs)
        # METRICS ASSESSMENT - SUB TASK
        self.check_metrics_ST(**kwargs)
        # PLOTs
        if self.report is not None:
            self.report.add_Section(
                keyword='check_plots', Title='Plots', level=1)
            self.addFigures_ST(**kwargs)
            self.addHKPlotsMatrix()
        # Update Report, raise flags, fill-in
        if self.log is not None:
            self.addFlagsToLog()
        if self.report is not None:
            self.addFlagsToReport()

    def check_HK_ST(self):
        """ """
        HKKeys = self.HKKeys
        if self.report is not None:
            self.report.add_Section(keyword='check_HK', Title='HK', level=1)
            HKKeys_tex = st.replace(HKKeys.__repr__(), '_', '\_')
            self.report.add_Text(['Selected HK Keys: %s' % HKKeys_tex, '\\newline'])

        report_HK_perf = self.check_HK(HKKeys, reference='command', limits='P', tag='Performance',
                                       doReport=self.report is not None,
                                       doLog=self.log is not None)
        HK_perf_ok = np.all(
            [value for key, value in report_HK_perf.items()])

        report_HK_safe = self.check_HK(HKKeys, reference='abs', limits='S', tag='Safe',
                                       doReport=self.report is not None,
                                       doLog=self.log is not None)
        HK_safe_ok = np.all(
            [value for ke, value in report_HK_safe.items()])

        if (not HK_perf_ok) or (not HK_safe_ok):
            self.dd.flags.add('HK_OOL')

    def debugtask(self):

        pass

#        if self.report is not None:
#            self.report.add_Section(
#                keyword='debug', Title='Debugging', level=0)
#
#        CCDs = ['CCD1','CCD2','CCD3']
#        Flu_lims = self.perflimits['Flu_lims']  # dict
#
#        _compliance_flu = self.check_stat_perCCDQandCol(
#            self.dd.mx['chk_med_inject'], Flu_lims, CCDs)
#
#        if self.report is not None:
#            self.addComplianceMatrix2Report(
#                _compliance_flu, label='COMPLIANCE FLUENCE:')

    def prepare_images(
            self,
            doExtract=True,
            doBadPixels=False,
            doMask=False,
            doOffset=False,
            doBias=False,
            doFF=False):
        """ """

        if self.report is not None:
            self.report.add_Section(
                keyword='prep_data', Title='Images Pre-Processing', level=0)

        if not doExtract:
            if self.report is not None:
                self.report.add_Text('Not extracting FITS files: Nothing done.')
            self.proc_histo['Extract'] = False
            return
        else:
            if self.report is not None:
                self.report.add_Text('Extracting FITS files to ccd.CCD objects.')

        def _loadCDP(cdpkey, msg):
            CDPData = calibration.load_FITS_CDPs(
                self.inputs['inCDPs'][cdpkey], 
                ccd.CCD,
                getallextensions=True,
                withpover=self.ccdcalc.withpover)
            cdpstr = self.inputs['inCDPs'][cdpkey].__str__()
            cdpstr = st.replace(cdpstr, ',', ',\n')
            if self.log is not None:
                self.log.info(msg)
                self.log.info(cdpstr)
            if self.report is not None:
                self.report.add_Text(msg)
                self.report.add_Text(cdpstr, verbatim=True)
            return CDPData

        def _reportNotFound(reportobj, msg):
            if reportobj is not None:
                reportobj.add_Text(msg)

        if doBadPixels:
            self.proc_histo['BadPixels'] = True
            if self.report is not None:
                self.report.add_Text('Masking out saturated and missing (zero) pixels.')

        if doMask and 'Mask' in self.inputs['inCDPs']:
            # self.inputs['inCDPs']['Mask']['CCD%i']
            MaskData = _loadCDP('Mask', 'Loading and applying Cosmetics Mask...')
            self.proc_histo['Masked'] = True
        elif doMask and 'Mask' not in self.inputs['inCDPs']:
            NotFoundMsg = 'Cosmetics Mask not Found!'
            self.log.info(NotFoundMsg)
            _reportNotFound(self.report, NotFoundMsg)
            doMask = False
            MaskData = None
        else:
            MaskData = None

        if doOffset:
            self.proc_histo['SubOffset'] = True
            offsetkwargs = self.inputs['preprocessing']['offsetkwargs']

            if self.report is not None:
                self.report.add_Text('Subtracting Offset.')
                msg = 'offsetkwargs=%s' % offsetkwargs.__repr__()
                self.report.add_Text(msg, verbatim=True)
        else:
            offsetkwargs = {}

        if doBias and 'bias' in self.inputs['inCDPs']:
            BiasData = _loadCDP('Bias', 'Loading And Subtracting Bias Structure...')
            self.proc_histo['SubBias'] = True
        elif doBias and 'bias' not in self.inputs['inCDPs']:
            NotFoundMsg = 'Bias Structure not found!'
            self.log.info(NotFoundMsg)
            _reportNotFound(self.report, NotFoundMsg)
            doBias = False
            BiasData = None
        else:
            BiasData = None

        if doFF and 'FF' in self.inputs['inCDPs']:
            FFData = _loadCDP('FF', 'Loading Flat-Field Maps...')
            self.proc_histo['FF'] = True
        elif doFF and 'FF' not in self.inputs['inCDPs']:
            NotFoundMsg = 'FFs no found!'
            self.log.info(NotFoundMsg)
            _reportNotFound(self.report, NotFoundMsg)
            doFF = False
            FFData = None
        else:
            FFData = None

        # Initialize new columns

        Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
        self.dd.initColumn('ccdobj_name', Cindices,
                           dtype='S100', valini='None')

        DDindices = copy.deepcopy(self.dd.indices)

        #nObs,nCCD,nQuad = DDindices.shape
        #Quads = DDindices[2].vals

        nObs = DDindices.get_len('ix')
        # nObs = 3  # TESTS!
        #print 'TESTS: task.prepare_images: LIMITTING TO 3 IMAGES!'

        CCDs = DDindices.get_vals('CCD')


        if not self.drill:

            picklespath = self.inputs['subpaths']['ccdpickles']

            arglist = []

            mgr = mp.Manager()
            queue = mgr.Queue()

            for iObs in range(nObs):
                arglist.append([queue, self.dd, self.ogse, self.inputs,
                                iObs, nObs, CCDs, picklespath, doBadPixels,
                                doMask, MaskData,
                                doOffset, offsetkwargs,
                                doBias, BiasData,
                                doFF, FFData])

            
            #prepare_one_image(*arglist[0]) # TEST
            #stop()

            pool = mp.Pool(processes=self.processes)

            for i in range(len(arglist)):
                pool.apply_async(prepare_one_image, args=arglist[i])
            pool.close()
            pool.join()

            replies = []
            while not queue.empty():
                replies.append(queue.get())

            for reply in replies:
                iObs, jCCD, ccdobj_name = reply
                self.dd.mx['ccdobj_name'][iObs, jCCD] = ccdobj_name

        return None

    def add_inputs_to_report(self):
        """ """

        self.report.add_Text('\\textbf{Test Inputs}')

        caption = 'Inputs of Task %s , Test $%s$' % (
            self.name, self.inputs['test'])
        #ncaption = st.replace(caption,'_','\\_')

        names = ['Parameter', 'Value']

        keys = list(self.inputs.manifesto.keys())
        values = []

        for key in keys:
            _val = self.inputs[key]
            # if isinstance(_val, dict):
            if key in ['structure', 'subpaths', 'perflimits', 'inCDPs', 'preprocessing']:
                values.append('Too long dict()')
            else:
                _val = _val.__repr__()
                #n_val = st.replace(_val,'_','\\_')
                n_val = st.replace(_val, '&', '\\&')
                values.append(n_val)

        tDict = OrderedDict(Parameter=keys, Value=values)
        #formats = dict(Parameter='s',Value='char')

        self.report.add_Table(tDict, names=names,
                              caption=caption, col_align='|l|X|')

    def get_data_inventory_table(self):

        tDict = OrderedDict()
        tDict['ObsID'] = self.dd.mx['ObsID'][:].copy()
        tDict['exptime'] = self.dd.mx['exptime'][:, 0].copy()
        tDict['chinj'] = self.dd.mx['chinj'][:, 0].copy()
        tDict['v_tpump'] = self.dd.mx['v_tpump'][:, 0].copy()
        tDict['s_tpump'] = self.dd.mx['s_tpump'][:, 0].copy()
        tDict['source'] = self.dd.mx['source'][:, 0].copy()
        tDict['wave'] = self.dd.mx['wave'][:, 0].copy()

        return tDict

    def add_data_inventory_to_report(self, tDict):
        """ """

        self.report.add_Text('\\textbf{Test Data}')

        caption = 'Data Used by Task %s , Test $%s$. Datapath = "%s"' % (
            self.name, self.inputs['test'], self.inputs['datapath'])
        #ncaption = st.replace(caption,'_','\\_')

        names = list(tDict.keys())

        self.report.add_Table(tDict, names=names,
                              caption=caption, longtable=True)

    def get_time_tag(self):
        from vison.support.vistime import get_time_tag
        return get_time_tag()

    def pack_CDP_to_dd(self, cdp, cdp_key):
        self.dd.products[cdp_key] = os.path.join(
            cdp.path, '%s.pick' % cdp.rootname)
