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

from vison import __version__
from vison.support import logger as lg
from vison.point import PSF0X,FOCUS00
from lib import get_time_tag

# END IMPORT

task_dict = dict(PSF0X=PSF0X,FOCUS00=FOCUS00)

defaults = dict(BLOCKID='R00P00CC000000',CHAMBER='A')


class Pipe(object):
    """Master Class of FM-analysis """
    
    
    def __init__(self,inputdict,dolog=True):
        """ """
                
        # inputsdict:
        #    TASKSlist
        
        self.inputs = defaults
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        self.BLOCKID=self.inputs['BLOCKID'] # BLOCK (ROE+RPSU+CCDs) under test
        self.CHAMBER=self.inputs['CHAMBER']
        self.ID = 'FM%s' % get_time_tag()  # ID of the analysis "session"

        if dolog:
            self.logf = 'Calib_%s.log' % self.ID
            self.log = lg.setUpLogger(self.logf)
            self.log.info(['\n\nStarting FM Calib. Pipeline',
                          'Pipeline ID: %s' % self.ID,
                          'BLOCK ID: %s' % self.BLOCKID,
                          'Chamber: %s\n' % self.CHAMBER,
                          'vison version: %s\n' % __version__,
                          'Tasks: %s\n' % ( ('%s,'*len(self.tasks)) % tuple(self.tasks))[0:-1]])
        self.log = None
        
         
    def run(self):
        """ """
        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        if not os.path.exists(resultsroot):
            os.system('mkdir %s' % resultsroot)
        
        
        if self.log is not None: self.log.info('\n\nResults will be saved in: %s\n' % resultsroot)
        
        for taskname in tasknames:
            
            task = task_dict[taskname]
            taskinputs = self.inputs[taskname]
            
            msg = ['\n\nRunning Task: %s\n' % taskname]
            msg += ['Inputs:']
            for key in taskinputs:
                msg += ['%s = %s' % (key,str(taskinputs[key]))]
            
            if self.log is not None: self.log.info(msg)
            
            taskinputs['resultspath'] = os.path.join(resultsroot,taskinputs['resultspath'])
            
            task.run(taskinputs,log=self.log)
            
    def wait_and_run(self):
        """ """
        
