# -*- coding: utf-8 -*-
"""

This is the main script that will orchestrate the analysis of
Euclid-VIS FM Calibration Campaign.

The functions of this module are:
    
    - Take inputs as to what data is to be analyzed, and what analysis scripts
      are to be run on it.
    - Set the variables necessary to process this batch of FM calib. data.
    - Start a log of actions to keep track of what is being done.
    - Provide inputs, execute the analysis scripts and report location of analysis 
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

@author: raf
"""

# IMPORT STUFF
from pdb import  set_trace as stop
from copy import copy

from vissim.support import logger as lg
from vissim.FMcalib import PSF02

import datetime

# END IMPORT

task_dict = dict(PSF02=PSF02)

def get_time_tag():
    """ """
    t = datetime.datetime.now()
    s = t.strftime('%Y%m%d_%H%M%S')
    return s

class FMpipeline(object):
    """Master Class of FM-analysis """
    
    
    def __init__(self,inputdict,dolog=True):
        """ """
                
        # inputsdict:
        #    TASKSlist
        
        self.inputs = inputdict        
        self.tasks = self.inputs['tasks']
        self.BLOCKID=self.inputs['BLOCKID'] # BLOCK (ROE+RPSU+CCDs) under test
        
        self.ID = 'FM%s' % get_time_tag()  # ID of the analysis "session"

        if dolog:
            self.logf = 'Calib_%s.log' % self.ID
            self.log = lg.setUpLogger(self.logf)
            self.log.info('\n\nStarting FM Calib. Pipeline\n\n')
            self.log.info('Pipeline ID: %s' % self.ID)

        
    def run(self):
        """ """
        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        self.log.info('\nResults will be saved in: %s' % resultsroot)
        
        for taskname in tasknames:
            
            self.log.info('\n\nRunning Task: %s' % taskname)
            
            task = task_dict[taskname]
            taskinputs = self.inputs[taskname]
            
            task.run(self.log,taskinputs)
            
            stop()
        