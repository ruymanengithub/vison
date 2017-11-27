#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Generic Task (Test) Class.

Created on Tue Nov 14 14:20:04 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import datetime
import string as st
from collections import OrderedDict

from vison.support.report import Report
from vison.support import vistime
from vison.support import files
import lib as pilib

# END IMPORT

isthere = os.path.exists

class Task(object):
    """ """
    
    from task_lib import check_HK
    
    def __init__(self,inputs,log=None):
        """ """
        
        self.inputs = self.feeder(inputs)
        self.log = log
        self.name = ''
        self.HKKeys = []
        self.perflimits = dict()
        
    def __call__(self):
        """Generic test master function."""
        
        # INPUTS
        
        subtasks = self.subtasks
        OBSID_lims = self.inputs['OBSID_lims']
        explogf = self.inputs['explogf']
        datapath = self.inputs['datapath']
        resultspath = self.inputs['resultspath']
        try: _paths = self.inputs['subpaths']
        except: _paths = dict()
        elvis = self.inputs['elvis']
        testkey = self.inputs['test']
        
        
        DataDictFile = os.path.join(resultspath,'%s_DataDict.pick' % testkey)
        reportobjFile = os.path.join(resultspath,'%s_Report.pick' % testkey)
        
        if not isthere(resultspath):
            os.system('mkdir %s' % resultspath)
        
        for _pathkey in _paths:
            subpath = os.path.join(resultspath,_paths[_pathkey])            
            if not isthere(subpath): os.system('mkdir %s' % subpath)
            _paths[_pathkey] = subpath
        
        self.inputs['subpaths'] = _paths
        
        
        structure = self.inputs['structure']
            
        try: reportroot = self.inputs['reportroot']
        except KeyError: reportroot = '%s_report' % testkey
        
        try: cleanafter = self.inputs['cleanafter']
        except KeyError: cleanafter = False
        
        todo_flags = self.inputs['todo_flags']
        
        
        if todo_flags['init']:
            
            
            # Let's start from scratch
            
            if os.path.exists(DataDictFile): os.system('rm %s' % DataDictFile)
            if os.path.exists(reportobjFile): os.system('rm %s' % reportobjFile)
            
            os.system('rm %s/*' % self.inputs['resultspath'])
            if 'subpaths' in self.inputs:
                for _pathkey in self.inputs['subpaths'].keys():
                    subpath = self.inputs['subpaths'][_pathkey]
                    os.system('rm %s/*' % subpath)
            
        
            # Initialising Report Object
        
            if todo_flags['report']:
                self.report = Report(TestName=self.name)
            else:
                self.report = None
        
            # META-DATA WORK
            
            
            # Filter Exposures that belong to the test
        

            explog, checkreport = self.filterexposures(structure,explogf,datapath,OBSID_lims,
                                         elvis)
            
    
            if self.log is not None:
                self.log.info('%s acquisition consistent with expectations: %s' % (testkey,checkreport['checksout']))
                if len(checkreport['failedcols'])>0:
                    self.log.info('%s failed columns: %s' % (testkey,checkreport['failedcols']))
                if len(checkreport['failedkeys'])>0:
                    self.log.info('%s failed keys: %s' % (testkey,checkreport['failedkeys']))
            
            # Adding Time Axis            
            
            explog['time'] = np.array(map(vistime.get_dtobj,explog['date'])).copy()

            
            # Building DataDict 
            
            self.dd = pilib.DataDict_builder(explog,self.inputs,structure)
            
            if not checkreport['checksout']:self.dd.flags.add('MISSDATA')
            
            # Add HK information
            self.addHK_2_dd()
            
            #self.save_progress(dd,reportobj,DataDictFile,reportobjFile)
            self.save_progress(DataDictFile,reportobjFile)
            
        else:
            
            #dd, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
            self.recover_progres(DataDictFile,reportobjFile)
        
        
        # DATA-WORK and ANALYSIS        
        
        for subtask in subtasks:
            
            subtaskname, subtaskmethod = subtask
            
            if self.log is not None:
                self.log.info('%s: %s' % (subtaskname,subtaskmethod.__module__))
            
            if todo_flags[subtaskname]:
                
                tini = datetime.datetime.now()
                subtaskmethod()
                tend = datetime.datetime.now()
                dtm = ((tend-tini).seconds)/60.
                if self.log is not None: 
                    self.log.info('%.1f minutes in running Sub-task: %s' % (dtm,subtaskname))
                
                self.save_progress(DataDictFile,reportobjFile)
            else:
                self.recover_progress(DataDictFile,reportobjFile)
        

        # Write automatic Report of Results

        if todo_flags['report']:
            self.report.doreport(reportroot,cleanafter)
            outfiles = self.report.writeto(reportroot,cleanafter)
            
            for outfile in outfiles:
                os.system('mv %s %s/' % (outfile,resultspath))
        

        self.save_progress(DataDictFile,reportobjFile)
        
        if self.log is not None:
            self.log.info('Finished %s' % self.name)

   
    def addHK_2_dd(self):
        """ """
        self.dd = pilib.addHK(self.dd,self.HKKeys,elvis=self.elvis)
       
    def save_progress(self,DataDictFile,reportobjFile):
        """ """
        files.cPickleDumpDictionary(self.dd,DataDictFile)
        files.cPickleDump(self.report,reportobjFile)
        
        
    def recover_progress(self,DataDictFile,reportobjFile):
        """ """
        self.dd = files.cPickleRead(DataDictFile)
        self.report = files.cPickleRead(reportobjFile)
    
    
    def addFigure2Report(self,figkey):
        """ """
        figobj = self.figdict[figkey]
        figname = figobj.figname
        texfraction = figobj.texfraction
        caption = figobj.caption
        
        assert os.path.exists(figname)
        epsname = '%s.eps' % os.path.splitext(figname)[0]
        os.system('convert %s %s' % (figname,epsname))
        self.report.add_Figure(epsname,texfraction=texfraction,
                               caption=caption,label=figkey)
        
    def doPlot(self,figkey,**kwargs):
        """ """
        figobj = self.figdict[figkey]()
        figobj.configure(**kwargs)
        figobj.build_data(self)
        figobj.plot()
        self.figdict[figkey] = figobj
         
    def addComplianceMatrix2Log(self,complidict,label=''):
        """ """
        st_compl = complidict.__str__()
        self.log.info('%s\n%s' % (label,st_compl))
                
    
    def addComplianceMatrix2Report(self,complidict,label=''):
        """ """
        nicelabel = st.replace(label,' ','\ ')
        st_compl = complidict.__str__()
        nice_st_compl = st.replace(st_compl,'False','$\\textcolor{red}{\\bf{False}}$')
        msgList = ['$\\bf{%s}$' % nicelabel,
        '%s' % nice_st_compl,
        '\\']
        self.report.add_Text(msgList)
        
    def IsComplianceMatrixOK(self,complidict):
        """ """

        def traverse_tree(dictionary,isOK):
        
            for key,value in dictionary.items():
                #print 'Upper: %s' % key
                if isinstance(value,(dict,OrderedDict)):
                    #print key,value
                    isOK = isOK and traverse_tree(value,isOK)
                else:
                    #print key,value
                    isOK = isOK and value
            return isOK
        
        isOK = traverse_tree(complidict,True)
        
        return isOK
        
        
    def addFlagsToLog(self):
        """ """
        flagstxt = st.join(self.dd.flags.getFlagsOnList(),', ')
        self.log.info('FLAGS ON}:\n%s' % flagstxt)
    
    def addFlagsToReport(self):
        """ """
        niceflagnames = [st.replace(item,'_','\_') for item in self.dd.flags.getFlagsOnList()]
        flagstxt = st.join(niceflagnames,', ')
        msgList = ['$\\bf{FLAGS\ ON}$: ',flagstxt]
        self.report.add_Text(msgList)