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


from vison import __version__
from vison.support import files
from vison.support.report import Report
# END IMPORT

class FpaTask(object):
    """ """
    
    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        
        self.ID = None
        if 'ID' in inputs:
            self.ID = inputs['ID']
        self.processes = 1
        if 'processes' in inputs:
            self.processes = inputs['processes']
        
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
        if not hasattr(self,'subtasks'):
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
        
        
        _inputs['todo_flags']= self.init_todo_flags()

        _inputs.update(inputs)
        
        self.inputs.update(_inputs)
        

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
        
    
    def __call__(self):
        """Generic test master function."""
        
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
            for _, subpath in self.inputs['subpaths'].iteritems():
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
                    dtm = ((tend-tini).seconds)/60.
                    if self.log is not None:
                        self.log.info(
                            '%.1f minutes in running Sub-task: %s' % (dtm, subtaskname))
                    #print 'saving progress!'
                    self.save_progress(DataDictFile, reportobjFile)
                except:
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
            outfiles = self.report.doreport(reportroot, cleanafter=cleantexafter, silent=True) # commented on TESTS
            #outfiles = self.report.doreport(reportroot, cleanafter=False, silent=False) # TESTS
            #stop() # TESTS

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
            print msg_trbk
    
    def ingest_data(self):
        """ """
        raise NotImplementedError
        
    def save_progress(self, DataDictFile, reportobjFile):
        """ """
        files.cPickleDumpDictionary(self.dd, DataDictFile)
        files.cPickleDump(self.report, reportobjFile)
        #csvFile = st.replace(DataDictFile,'.pick','.csv')
        #self.dd.saveToFile(csvFile)

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
                execline1 = "find %s/ -type f -name '*.fits' -exec sh -c '%s' {} \;" % (subpath,'rm "$0"')
                os.system(execline1)
                execline2 = "find %s/ -type f -name '*.pick' -exec sh -c '%s' {} \;" % (subpath,'rm "$0"')
                os.system(execline2)
                if self.log is not None:
                    self.log.info('\nCleared contents [.fits/.pick] of %s!' % subpath)
