# -*- coding: utf-8 -*-
"""

This is the main script that will orchestrate the analysis of
Euclid-VIS FM Ground Calibration Campaign.

The functions of this module are:
    
    - Take inputs as to what data is to be analyzed, and what analysis scripts
      are to be run on it.
    - Set the variables necessary to process this batch of FM calib. data.
    - Start a log of actions to keep track of what is being done.
    - Provide inputs to scripts, execute the analysis scripts and report location of analysis 
      results.
    
Some Guidelines for Development:

    - Input data is "sacred": read-only.
    - Each execution of Master must have associated a unique ANALYSIS-ID.
    - All the Analysis must be divided in TASKS. TASKS can have SUB-TASKS.
    - All data for each TASK must be under a single directory (TBC).
    - All results from the execution of FMmaster must be under a single directory 
      with subdirectories for each TASK run.
    - A subfolder of this root directory will contain the logging information:
         inputs, outputs, analysis results locations.
    

Created on Wed Jul 27 12:16:40 2016

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import  set_trace as stop
from copy import copy
import os
import numpy as np

from vison import __version__
from vison.support import logger as lg
from vison.support.report import Report
from vison.dark import BIAS01,DARK01
from vison.point import FOCUS00,PSF0X
from vison.flat import NL01, PTC0X, FLAT0X
from vison.inject import CHINJ01,CHINJ02
from vison.pump import TP01, TP02
from vison.other import PERSIST01 as PER01
#from lib import get_time_tag
from vison.pipe import lib as pilib

# END IMPORT

isthere = os.path.exists

defaults = dict(BLOCKID='R00P00CC000000',CHAMBER='A')


class Pipe(object):
    """Master Class of FM-analysis """
    
    from vison.dark import BIAS01
    #from vison.dark.DARK01 import dark01
    #from vison.inject.CHINJ01 import chinj01
    #from vison.inject.CHINJ02 import chinj02
    #from vison.flat.FLAT0X import flat0x
    #from vison.flat.NL01 import nl01
    #from vison.flat.PTC0X import ptc0x
    #from vison.pump.TP01 import tp01
    #from vison.pumpTP02 import tp02
    #from vison.other.PERSIST01 import persist01
    
    Test_dict = dict(BIAS01=BIAS01)
    
    
    def __init__(self,inputdict,dolog=True):
        """ """
        
        
#        self.task_dict = dict(BIAS01=self.bias01)
#        ,DARK1=self.dark01,
#                 CHINJ01=self.chinj01,CHINJ02=self.chinj02,
#                 NL01=self.nl01,                 
#                 PSF01_590=self.psf0x,PSF01_640=self.psf0x,
#                 PSF01_730=self.psf0x,PSF01_800=self.psf0x,
#                 PSF01_880=self.psf0x,
#                 FOCUS00_590=self.focus00,FOCUS00_640=self.focus00,
#                 FOCUS00_730=self.focus00,FOCUS00_800=self.focus00,
#                 FOCUS00_890=self.focus00)
#        
        self.inputs = defaults
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        self.BLOCKID=self.inputs['BLOCKID'] # BLOCK (ROE+RPSU+CCDs) under test
        self.CHAMBER=self.inputs['CHAMBER']
        self.ID = 'FM%s' % pilib.get_time_tag()  # ID of the analysis "session"

        if dolog:
            self.logf = 'Calib_%s.log' % self.ID
            self.log = lg.setUpLogger(self.logf)
            self.log.info(['\n\nStarting FM Calib. Pipeline',
                          'Pipeline ID: %s' % self.ID,
                          'BLOCK ID: %s' % self.BLOCKID,
                          'Chamber: %s\n' % self.CHAMBER,
                          'vison version: %s\n' % __version__,
                          'Tasks: %s\n' % ( ('%s,'*len(self.tasks)) % tuple(self.tasks))[0:-1]])
        else:
            self.log = None
        
         
    def run(self):
        """ """
        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        if not os.path.exists(resultsroot):
            os.system('mkdir %s' % resultsroot)
        
        
        if self.log is not None: self.log.info('\n\nResults will be saved in: %s\n' % resultsroot)
        
        for taskname in tasknames:
            
            taskinputs = self.inputs[taskname]
            
            msg = ['\n\nRunning Task: %s\n' % taskname]
            msg += ['Inputs:']
            for key in taskinputs:
                msg += ['%s = %s' % (key,str(taskinputs[key]))]
            
            if self.log is not None: self.log.info(msg)
            
            taskinputs['resultspath'] = os.path.join(resultsroot,taskinputs['resultspath'])
            
            self.dotask(taskname,taskinputs)
    
    def dotask(self,taskname,inputs):
        """Generic test master function."""
    
        # INPUTS
        
        #todo_flags = dict(init=True,prep=True,basic=True,meta=True,report=True)
        
        #self.Test_dict[taskname].feeder(inputs)
        Test = self.Test_dict[taskname]
        inputs = Test.feeder(inputs)
        
        subtasks = inputs['subtasks']
        OBSID_lims = inputs['OBSID_lims']
        explogf = inputs['explogf']
        datapath = inputs['datapath']
        resultspath = inputs['resultspath']
        elvis = inputs['elvis']
        testkey = inputs['testkey']
        
        DataDictFile = os.path.join(resultspath,'%s_DataDict.pick' % testkey)
        reportobjFile = os.path.join(resultspath,'%s_Report.pick' % testkey)
        
        if not isthere(resultspath):
            os.system('mkdir %s' % resultspath)
        
        structure = inputs['structure']
            
        try: reportroot = inputs['reportroot']
        except KeyError: reportroot = '%s_report' % testkey
        
        try: cleanafter = inputs['cleanafter']
        except KeyError: cleanafter = False
        
        todo_flags = inputs['todo_flags']
        
        
        if todo_flags['init']:
        
            # Initialising Report Object
        
            if todo_flags['report']:
                reportobj = Report(TestName=taskname)
            else:
                reportobj = None
        
            # META-DATA WORK
            
            
            # Filter Exposures that belong to the test
        
            #DataDict, isconsistent = Test.filterexposures(structure,explogf,datapath,OBSID_lims,
            #                             elvis)
            
            explog, checkreport = Test.filterexposures(structure,explogf,datapath,OBSID_lims,
                                         elvis)
    
            if self.log is not None:
                self.log.info('%s acquisition consistent with expectations: %s' % (testkey,checkreport['checksout']))
                if len(checkreport['failedcols'])>0:
                    self.log.info('%s failed columns: %s' % (testkey,checkreport['failedcols']))
                if len(checkreport['failedkeys'])>0:
                    self.log.info('%s failed keys: %s' % (testkey,checkreport['failedkeys']))
            
            # Adding Time Axis            
            
            explog['time'] = np.array(map(pilib.get_dtobj,explog['DATE'])).copy()

            # Filling in the .fits extension
            
            rootFile_name = explog['File_name'].copy()
            File_name  = ['%s.fits' % item for item in rootFile_name]
            explog['Files'] = np.array(File_name).copy()
            
            
            # Building DataDict 
            
            
            DataDict = pilib.DataDict_builder(explog,inputs,structure)
            
            
            HKKeys = Test.HKKeys
            
            # Add HK information
            DataDict = pilib.addHK(DataDict,HKKeys,elvis=elvis)
            
            
            pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
            
        else:
            
            DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        
        # DATA-WORK and ANALYSIS        
        
        for subtask in subtasks:
            
            subtaskname, subtaskfunc = subtask
            
            if todo_flags[subtaskname]:
                DataDict,reportobj = subtaskfunc(DataDict,reportobj,inputs,self.log)
                pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
            else:
                DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        

        # Write automatic Report of Results
        
        if todo_flags['report']:
            
            reportobj.doreport(reportroot,cleanafter)
            outfiles = reportobj.writeto(reportroot,cleanafter)
            
            for outfile in outfiles:
                os.system('mv %s %s/' % (outfile,resultspath))
        
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
        
        if self.log is not None:
            self.log.info('Finished %s' % taskname)
        
        
    
    
            
    def wait_and_run(self):
        """ """
        
