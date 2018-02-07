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
import datetime
import sys,traceback

from vison import __version__
from vison.support import logger as lg
from vison.support.report import Report
from vison.support import vistime
#from lib import get_time_tag
from vison.pipe import lib as pilib
from vison.support import context
# END IMPORT

isthere = os.path.exists

defaults = dict(BLOCKID='R00P00CC000000',CHAMBER='A')

waittime = 120 # seconds


class Pipe(object):
    """Master Class of FM-analysis """
    
    from vison.dark.BIAS01 import BIAS01
    from vison.dark.DARK01 import DARK01
    from vison.flat.NL01 import NL01
    from vison.flat.FLAT0X import FLAT0X
    from vison.flat.PTC0X import PTC0X
    from vison.inject.CHINJ00 import CHINJ00
    from vison.inject.CHINJ01 import CHINJ01
    from vison.inject.CHINJ02 import CHINJ02
    from vison.other.PERSIST01 import PERSIST01
    from vison.point.PSF0X import PSF0X
    from vison.point.FOCUS00 import FOCUS00
    from vison.point.PSF01_PANCHRO import PSF01_PANCHRO
    from vison.pump.TP00 import TP00
    from vison.pump.TP01 import TP01
    from vison.pump.TP02 import TP02
    
    
    Test_dict = dict(BIAS01=BIAS01,DARK01=DARK01,
                     NL01=NL01,FLAT01=FLAT0X,
                     PTC01=PTC0X,
                     CHINJ00=CHINJ00,CHINJ01=CHINJ01,CHINJ02=CHINJ02,
                     TP00=TP00,
                     TP01=TP01,TP02=TP02,
                     PERSIST01=PERSIST01,
                     PSF01_PANCHRO=PSF01_PANCHRO)
    
    for wave in [0,590,640,730,880]:
        Test_dict['FLAT02_%i' % wave] = FLAT0X
    for wave in [590,640,730,800,880,0]:
        Test_dict['PTC02_%i' % wave] = PTC0X
    for temp in [150,156]:
        Test_dict['PTC02_%iK' % temp] = PTC0X
    for wave in [590,640,730,800,880]:
        Test_dict['FOCUS00_%i' % wave] = FOCUS00
    for wave in [590,640,730,800,880]:
        Test_dict['PSF01_%i' % wave] = PSF0X
    for temp in [150,156]:
        Test_dict['PSF02_%iK' % temp] = PSF0X
    
    def __init__(self,inputdict,dolog=True,drill=False,debug=False,startobsid=0):
        """ """
        
        self.inputs = defaults
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        self.BLOCKID=self.inputs['BLOCKID'] # BLOCK (ROE+RPSU+CCDs) under test
        self.CHAMBER=self.inputs['CHAMBER']
        self.drill = drill
        self.debug = debug
        self.startobsid = startobsid
        
        if self.debug:
            self.ID = 'PipeDebug'
        else:
            self.ID = 'FM%s' % vistime.get_time_tag()  # ID of the analysis "session"
        
        self.inputs['ID'] = self.ID

        if dolog:
            self.logf = 'Calib_%s.log' % self.ID
            
            if os.path.exists(self.logf): os.system('rm %s' % self.logf)
            
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
        
        extinputs = copy.deepcopy(self.inputs)
        alltasks = extinputs['tasks']
        extinputs.pop('tasks')
        for _taskname in alltasks: extinputs.pop(_taskname)
        taskinputs.update(extinputs)
        
        msg = ['\n\nRunning Task: %s\n' % taskname]
        msg += ['Inputs:\n']
        for key in taskinputs:
            msg += ['%s = %s' % (key,str(taskinputs[key]))]
            
        if self.log is not None: self.log.info(msg)

        tini = datetime.datetime.now()
        
        
        self.dotask(taskname,taskinputs,drill=self.drill,debug=self.debug)
        
        tend = datetime.datetime.now()
        dtm = ((tend-tini).seconds)/60.
        if self.log is not None: 
            self.log.info('%.1f minutes in running Task: %s' % (dtm,taskname))
        
    
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
    
    def wait_and_run(self,explogf,elvis=context.elvis):
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
            test = Test(taskinputs)
            taskinputs = copy.deepcopy(test.inputs)
            
            
            structure = taskinputs['structure']
            
            Ncols = structure['Ncols']
            
            Nframes = 0
            for ic in range(1,Ncols+1): Nframes += structure['col%i' % ic]['frames']
            
            tasksequence.append((taskname,testkey,Nframes))
            
        # Launching tasks
        
        fahrtig = False
        
        while not fahrtig:
            
            explog = pilib.loadexplogs(explogf,elvis)
            
            if self.startobsid>0:
                ixstart = np.where(explog['ObsID']==self.startobsid)[0][0]
                explog = explog[ixstart:].copy()
            
            for it,taskitem in enumerate(tasksequence):
                
                taskname,testkey,Nframes = taskitem
                available = pilib.coarsefindTestinExpLog(explog,testkey,Nframes)
                
                if available:
                    #print '%s available, doing nothing!' % taskname # TESTS
                    self.launchtask(taskname)
                    tasksequence.pop(it)
            
            sleep(waittime)
            
            #print tasksequence
            
            if len(tasksequence) == 0:
                fahrtig = True    
        
        
        return None

    
    def dotask(self,taskname,inputs,drill=False,debug=False):
        """Generic test master function."""
        
        try:
            Test = self.Test_dict[taskname](inputs,self.log,drill,debug)
            Test.ID = self.ID
            Test()
        except:
            self.catchtraceback()
            if self.log is not None:
                self.log.info('TASK "%s@%s" FAILED, QUITTING!' % (taskname,self.Test_dict[taskname].__module__))
            else:
                print 'TASK "%s@%s" FAILED, QUITTING!' % (taskname,self.Test_dict[taskname].__module__)
    
    def catchtraceback(self):
        """ """
        msg_trbk = traceback.format_exc()
        if self.log is not None:
            self.log.info(msg_trbk)
        else:
            print msg_trbk

        

