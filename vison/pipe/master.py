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
import copy
import os
import numpy as np
from time import sleep

from vison import __version__
from vison.support import logger as lg
from vison.support.report import Report
from vison.dark import BIAS01,DARK01
from vison.point import FOCUS00,PSF0X
from vison.flat import NL01, PTC0X, FLAT0X
from vison.inject import CHINJ01,CHINJ02
from vison.pump import TP01, TP02
from vison.other import PERSIST01 as PER01
from vison.support import vistime
#from lib import get_time_tag
from vison.pipe import lib as pilib

# END IMPORT

isthere = os.path.exists

defaults = dict(BLOCKID='R00P00CC000000',CHAMBER='A')

waittime = 120 # seconds


class Pipe(object):
    """Master Class of FM-analysis """
    
    from vison.dark import BIAS01, DARK01
    from vison.flat import NL01, FLAT0X, PTC0X
    from vison.inject import CHINJ01, CHINJ02
    from vison.other import PERSIST01
    from vison.point import PSF0X, FOCUS00
    from vison.pump import TP01, TP02
    
    
    Test_dict = dict(BIAS01=BIAS01,DARK01=DARK01,
                     NL01=NL01,FLAT01=FLAT0X,
                     PTC01=PTC0X,
                     CHINJ01=CHINJ01,CHINJ02=CHINJ02,
                     TP01=TP01,TP02=TP02,
                     PERSIST01=PERSIST01)
    
    for wave in [590,640,880]:
        Test_dict['FLAT02_%i' % wave] = FLAT0X
    for wave in [590,640,730,880]:
        Test_dict['PTC02_%i' % wave] = PTC0X
    for temp in [150,156]:
        Test_dict['PTC02_%iK' % temp] = PTC0X
    for wave in [590,640,730,800,880]:
        Test_dict['FOCUS00_%i' % wave] = FOCUS00
    for wave in [590,640,800,880]:
        Test_dict['PSF01_%i' % wave] = PSF0X
    for temp in [150,156]:
        Test_dict['PSF02_%iK' % temp] = PSF0X
    
    def __init__(self,inputdict,dolog=True):
        """ """
        
        self.inputs = defaults
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        self.BLOCKID=self.inputs['BLOCKID'] # BLOCK (ROE+RPSU+CCDs) under test
        self.CHAMBER=self.inputs['CHAMBER']
        self.ID = 'FM%s' % vistime.get_time_tag()  # ID of the analysis "session"
        
        self.inputs['ID'] = self.ID

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
        
    def launchtask(self,taskname):
        """ """
        
        taskinputs = self.inputs[taskname]
        
        inputs =copy.deepcopy(self.inputs)
        alltasks = inputs['tasks']
        inputs.pop('tasks')
        for taskname in alltasks: inputs.pop(taskname)
        taskinputs.update(inputs)        
            
        msg = ['\n\nRunning Task: %s\n' % taskname]
        msg += ['Inputs:']
        for key in taskinputs:
            msg += ['%s = %s' % (key,str(taskinputs[key]))]
            
        if self.log is not None: self.log.info(msg)
        
        self.dotask(taskname,taskinputs)
        
        
    
    def run(self,explogf=None,elvis=None):
        """ """
        
        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        if not os.path.exists(resultsroot):
            os.system('mkdir %s' % resultsroot)
        
        if self.log is not None: self.log.info('\n\nResults will be saved in: %s\n' % resultsroot)
        
        for taskname in tasknames:
            
            taskinputs = self.inputs[taskname]
            taskinputs['resultspath'] = os.path.join(resultsroot,taskinputs['resultspath'])
            
            if explogf is not None:
                taskinputs['explogf'] = explogf
            if elvis is not None:
                taskinputs['elvis'] = elvis
            
            self.inputs[taskname] = taskinputs
            
            self.launchtask(taskname)
    
    def wait_and_run(self,explogf,elvis='6.1.0'):
        """ """
        
        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        if not os.path.exists(resultsroot):
            os.system('mkdir %s' % resultsroot)
        
        if self.log is not None: self.log.info('\n\nResults will be saved in: %s\n' % resultsroot)
        
        # Learn how many ObsIDs will generate each task
        
        tasksequence = []
        
        for taskname in tasknames:
            
            taskinputs = self.inputs[taskname]
            
            testkey = taskinputs['test']
            taskinputs['resultspath'] = os.path.join(resultsroot,taskinputs['resultspath'])
            taskinputs['explogf'] = explogf
            taskinputs['elvis'] = elvis            
            
            Test = self.Test_dict[taskname]
            taskinputs = Test.feeder(taskinputs)            
            
            structure = taskinputs['structure']
            
            Ncols = structure['Ncols']
            
            Nframes = 0
            for ic in range(1,Ncols+1): Nframes += structure['col%i' % ic]['frames']
            
            tasksequence.append((taskname,testkey,Nframes))
            
        # Launching tasks
        
        fahrtig = False
        
        while not fahrtig:
            
            explog = pilib.loadexplogs(explogf,elvis)
            
            for taskitem in tasksequence:
                
                taskname,testkey,Nframes = taskitem
        
                available = pilib.coarsefindTestinExpLog(explog,testkey,Nframes)
                
                if available:
                    self.launchtask(taskname)
                    tasksequence.pop(0)
            
            sleep(waittime)
            
            #print tasksequence
            
            if len(tasksequence) == 0:
                fahrtig = True    
        
        
        return None

    
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
        testkey = inputs['test']
        
        
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
        
            #dd, isconsistent = Test.filterexposures(structure,explogf,datapath,OBSID_lims,
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
            
            explog['time'] = np.array(map(vistime.get_dtobj,explog['date'])).copy()

            
            # Building DataDict 
            
            dd = pilib.DataDict_builder(explog,inputs,structure)
            
            
            HKKeys = Test.HKKeys
            
            # Add HK information
            dd = pilib.addHK(dd,HKKeys,elvis=elvis)
            
            
            pilib.save_progress(dd,reportobj,DataDictFile,reportobjFile)
            
        else:
            
            dd, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        
        # DATA-WORK and ANALYSIS        
        
        for subtask in subtasks:
            
            subtaskname, subtaskfunc = subtask
            
            if todo_flags[subtaskname]:
                dd,reportobj = subtaskfunc(dd,reportobj,inputs,self.log)
                pilib.save_progress(dd,reportobj,DataDictFile,reportobjFile)
            else:
                dd, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        

        # Write automatic Report of Results

        if todo_flags['report']:
            
            reportobj.doreport(reportroot,cleanafter)
            outfiles = reportobj.writeto(reportroot,cleanafter)
            
            for outfile in outfiles:
                os.system('mv %s %s/' % (outfile,resultspath))
        

        pilib.save_progress(dd,reportobj,DataDictFile,reportobjFile)
        
        if self.log is not None:
            self.log.info('Finished %s' % taskname)
        
        

