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
    - All data for each TASK must be under a single day-folder.
    - All results from the execution of FMmaster must be under a single directory 
      with subdirectories for each TASK run.
    - A subfolder of this root directory will contain the logging information:
         inputs, outputs, analysis results locations.
    

Created on Wed Jul 27 12:16:40 2016

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import os
import numpy as np
from time import sleep
import datetime
import sys
import traceback
import glob
from collections import OrderedDict
import select
#import inspect

from vison import __version__
from vison.support import logger as lg
#from vison.support.report import Report
from vison.support import vistime
#from lib import get_time_tag
from vison.pipe import lib as pilib
from vison.support import context
from vison.support import utils
from vison.support import flags as flagsmodule
# END IMPORT

isthere = os.path.exists

defaults = dict(BLOCKID='R00P00CC000000', 
                CHAMBER='A_JUN18')

waittime = 120  # seconds
waitTO = 3600.*4. # seconds


class GenPipe(object):
    """Master Pipeline Class."""
    
    def __init__(self, inputdict, dolog=True, drill=False, debug=False,
                 cleanafter=False):
        """ """
        
        self.inputs = defaults.copy()
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        self.drill = drill
        self.debug = debug
        self.cleanafter = cleanafter
        self.completion = OrderedDict()

        if self.debug:
            self.ID = 'PipeDebug%s' % self.tag
        else:
            self.ID = 'FM%s%s' % (vistime.get_time_tag(),self.tag)  # ID of the analysis "session"

        self.inputs['ID'] = self.ID

        if dolog:
            self.logf = 'Calib_%s.log' % self.ID

            if os.path.exists(self.logf):
                os.system('rm %s' % self.logf)

            self.log = lg.setUpLogger(self.logf)
            self.log.info(['_', 'Starting FM Calib. Pipeline'] +
                          self._get_log_header())
        else:
            self.log = None
            
        self.pipe_session = dict(ID=self.ID,
                               processes=self.processes)

    def _get_log_header(self):
        log_header = [
            'Pipeline ID: %s' % self.ID,
            'vison version: %s\n' % __version__,
            'Tasks: %s\n' % self.tasks.__repr__()]
        return log_header

    def get_test(self, taskname, inputs=dict(), log=None, drill=False, debug=False, cleanafter=False):
        """ """
        testclass = self.Test_dict[taskname]
        #pathtotest = os.path.split(inspect.getfile(testclass))[0]
        #os.system('rm %s/*.pyc' % pathtotest) # avoids problems with .pyc files from previous tasks in the run... DIRTY HACK
        test = testclass(inputs=inputs, log=log, drill=drill, debug=debug,
                         cleanafter=cleanafter)
        return test

    def launchtask(self, taskname):
        """ """

        taskinputs = self.inputs[taskname]

        extinputs = copy.deepcopy(self.inputs)
        alltasks = extinputs['tasks']
        extinputs.pop('tasks')
        for _taskname in alltasks:
            extinputs.pop(_taskname)
        taskinputs.update(extinputs)

        taskinputs.update(self.pipe_session.copy())

        msg = ['\n\nRunning Task: %s\n' % taskname]
        msg += ['Inputs:\n']
        for key in taskinputs:
            msg += ['%s = %s' % (key, str(taskinputs[key]))]

        if self.log is not None:
            self.log.info(msg)
        
        print('Running: %s' % taskname)        
        
        taskreport = self.dotask(taskname, taskinputs,
                                 drill=self.drill, debug=self.debug,
                                 cleanafter=self.cleanafter)
        
        timemsg = '%.1f minutes in running Task: %s' %\
                          (taskreport['exectime'], taskname)
        errormsg = 'Task %s exited with Errors: %s' %\
                          (taskname, taskreport['Errors'])
        print(timemsg)
        print(errormsg)
        
        if self.log is not None:
            self.log.info(timemsg)
            self.log.info(errormsg)
        
        self.completion[taskname] = copy.deepcopy(taskreport)

    def run(self, explogf=None, elvis=None):
        """ """

        tini = datetime.datetime.now()

        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        if not os.path.exists(resultsroot):
            os.system('mkdir %s' % resultsroot)

        if self.log is not None:
            self.log.info('\n\nResults will be saved in: %s\n' % resultsroot)

        for taskname in tasknames:

            taskinputs = self.inputs[taskname]
            taskinputs['resultspath'] = os.path.join(
                resultsroot, taskinputs['resultspath'])

            if explogf is not None:
                taskinputs['explogf'] = explogf
            if elvis is not None:
                taskinputs['elvis'] = elvis

            self.inputs[taskname] = taskinputs

            self.launchtask(taskname)
            
            tpart = datetime.datetime.now()
            partdtm = ((tpart-tini).seconds)/60.
            
            partsummary = self.get_execution_summary(exectime=partdtm)
            if self.log is not None:
                self.log.info(partsummary)

        tend = datetime.datetime.now()
        Dtm = ((tend-tini).seconds)/60.

        summary = self.get_execution_summary(exectime=Dtm)
        if self.log is not None:
            self.log.info(summary)

    def get_execution_summary(self, exectime=None):
        """ """
        summary = ['_', '_', '_', '_',
                   '######################################################################'] +\
            self._get_log_header()
        
        for task in self.tasks:
            if task in self.completion:
                _compl = self.completion[task]
                summary += ['_', '%s' % task]
                summary += ['Executed in %.1f minutes' % _compl['exectime']]
                summary += ['OBSIDs range = [%i, %i]' % _compl['ObsID_range']]
                summary += ['Raised Flags = %s' % _compl['flags']]
                if _compl['Errors']:
                    summary += ['Executed with ERROR(s)']
        
        summary += ['_', '_', '_', '_',
                   '######################################################################']
        
        return summary
    
    def dotask(self, taskname, inputs, drill=False, debug=False, cleanafter=False):
        """Generic test master function."""
        
        strip_taskname = utils.remove_iter_tag(taskname,Full=False)

        tini = datetime.datetime.now()

        Errors = False
        
        try:
            test = self.get_test(strip_taskname, inputs=inputs, log=self.log, 
                                 drill=drill, debug=debug, cleanafter=cleanafter)
            
            Errors = test()  # test execution
        except:
            self.catchtraceback()
            Errors = True
            if self.log is not None:
                self.log.info('TASK "%s@%s" FAILED, QUITTING!' %
                              (taskname, self.Test_dict[strip_taskname].__module__))
            else:
                print 'TASK "%s@%s" FAILED, QUITTING!' % (
                    taskname, self.Test_dict[taskname].__module__)
        
        tend = datetime.datetime.now()
        dtm = ((tend-tini).seconds)/60.
        
        execlog = OrderedDict()
        execlog['Task'] = taskname
        execlog['Errors'] = Errors
        execlog['exectime'] = dtm
        try:
            execlog['ObsID_range'] = (test.dd.mx['ObsID'][:].min(),
                                   test.dd.mx['ObsID'][:].max() )
        except:
            execlog['ObsID_range'] = (-1,-1)
        
        try:
            flags = test.dd.flags.getFlagsOnList()
        except:
            flags = ['UNKNOWN']
        execlog['flags'] = flags.__repr__()
        
        test = None
        
        return execlog

    def catchtraceback(self):
        """ """
        msg_trbk = traceback.format_exc()
        if self.log is not None:
            self.log.info(msg_trbk)
        else:
            print msg_trbk


class Pipe(GenPipe):
    """Master Class of FM-analysis at block-level of assembly."""
    from vison.dark.BIAS0X import BIAS0X
    from vison.dark.DARK01 import DARK01
    from vison.flat.NL01 import NL01
    from vison.flat.NL02 import NL02
    from vison.flat.FLAT0X import FLAT0X
    from vison.flat.PTC0X import PTC0X
    from vison.flat.BF01 import BF01
    from vison.inject.CHINJ00 import CHINJ00
    from vison.inject.CHINJ01 import CHINJ01
    from vison.inject.CHINJ02 import CHINJ02
    from vison.other.PERSIST01 import PERSIST01
    from vison.point.PSF0X import PSF0X
    from vison.point.FOCUS00 import FOCUS00
    from vison.point.PSF01_PANCHRO import PSF01_PANCHRO
    from vison.pump.TP00 import TP00
    from vison.pump.TP01 import TP01
    from vison.pump.TP11 import TP11
    from vison.pump.TP02 import TP02
    from vison.pump.TP21 import TP21
    from vison.other.STRAY00 import STRAY00
    from vison.other.MOT_FF import MOT_FF
    from vison.other.MOT_WARM import MOT_WARM
    from vison.other.COSMETICS00 import COSMETICS00

    Test_dict = dict(BIAS01=BIAS0X, 
                     BIAS02=BIAS0X,
                     COSMETICS00=COSMETICS00,
                     DARK01=DARK01,
                     NL01=NL01, NL02=NL02, 
                     FLAT01=FLAT0X,FLAT_STB=FLAT0X,
                     PTC01=PTC0X,BF01=BF01,
                     CHINJ00=CHINJ00, CHINJ01=CHINJ01, CHINJ02=CHINJ02,
                     TP00=TP00,
                     TP01=TP01, TP11=TP11,
                     TP02=TP02, TP21=TP21,
                     PERSIST01=PERSIST01,
                     PSF01_PANCHRO=PSF01_PANCHRO,
                     STRAY00=STRAY00,
                     MOT_FF=MOT_FF,
                     MOT_WARM=MOT_WARM)
    
    for RD in [15.5, 16.0, 16.5]:
        Test_dict['NL02_%iR' % (RD*10,)] = NL02
    for wave in [0, 590, 640, 730, 880]:
        Test_dict['FLAT02_%i' % wave] = FLAT0X
    for wave in [590, 640, 730, 800, 880, 0]:
        Test_dict['PTC02_%i' % wave] = PTC0X
    for wave in [590, 640, 730, 800, 880, 0]:
        Test_dict['BF01_%i' % wave] = BF01
    for temp in [148, 158]:
        Test_dict['PTC02_%iK' % temp] = PTC0X
    for RD in [15.5, 16.0, 16.5]:
        Test_dict['PTC02_%iR' % (RD*10,)] = PTC0X
    for wave in [590, 640, 730, 800, 880]:
        Test_dict['FOCUS00_%i' % wave] = FOCUS00
    for wave in [590, 640, 730, 800, 880]:
        Test_dict['PSF01_%i' % wave] = PSF0X
    for temp in [148, 158]:
        Test_dict['PSF02_%iK' % temp] = PSF0X
    for wave in [590, 640, 2000, 730, 800, 880, 0]:
        Test_dict['PSFLUX00_%i' % wave] = PSF0X
    for wave in [590, 640, 2000, 730, 800, 880, 0]:
        Test_dict['FLATFLUX00_%i' % wave] = PTC0X

    def __init__(self, inputdict, dolog=True, drill=False, debug=False, startobsid=0,
                 processes=1,tag='',cleanafter=False):
        """ """
        
        self.inputs = defaults.copy()
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        # BLOCK (ROE+RPSU+CCDs) under test
        self.BLOCKID = self.inputs['BLOCKID']
        self.CHAMBER = self.inputs['CHAMBER']
        self.drill = drill
        self.debug = debug
        self.startobsid = startobsid
        self.processes = processes
        self.cleanafter = cleanafter
        self.tag = tag
        self.completion = OrderedDict()

        if self.debug:
            self.ID = 'PipeDebug%s' % self.tag
        else:
            self.ID = 'FM%s%s' % (vistime.get_time_tag(),self.tag)  # ID of the analysis "session"

        self.inputs['ID'] = self.ID

        if dolog:
            self.logf = 'Calib_%s.log' % self.ID

            if os.path.exists(self.logf):
                os.system('rm %s' % self.logf)

            self.log = lg.setUpLogger(self.logf)
            self.log.info(['_', 'Starting FM Calib. Pipeline'] +
                          self._get_log_header())
        else:
            self.log = None
            
        self.pipe_session = dict(ID=self.ID,
                               BLOCKID=self.BLOCKID,
                               CHAMBER=self.CHAMBER,
                               processes=self.processes)

    def _get_log_header(self):
        log_header = [
            'Pipeline ID: %s' % self.ID,
            'BLOCK ID: %s' % self.BLOCKID,
            'Chamber: %s\n' % self.CHAMBER,
            'vison version: %s\n' % __version__,
            'Tasks: %s\n' % self.tasks.__repr__()]
        return log_header


    
    def wait_and_run(self, dayfolder, elvis=context.elvis):
        """ """
        
        tasknames = self.tasks
        resultsroot = self.inputs['resultsroot']
        if not os.path.exists(resultsroot):
            os.system('mkdir %s' % resultsroot)
        
        if self.log is not None:
            self.log.info('\n\nResults will be saved in: %s\n' % resultsroot)
        
        # Learn how many ObsIDs will generate each task

        tasksequence = []
        

        for taskname in tasknames:

            taskinputs = self.inputs[taskname]

            testkey = taskinputs['test']
            taskinputs['resultspath'] = os.path.join(
                resultsroot, taskinputs['resultspath'])
            taskinputs['elvis'] = elvis
            taskinputs['datapath'] = dayfolder
                      
                      
            strip_taskname = utils.remove_iter_tag(taskname, Full=False)
            test = self.get_test(strip_taskname, taskinputs, log=self.log)
            taskinputs = copy.deepcopy(test.inputs)
            

            structure = taskinputs['structure']

            Ncols = structure['Ncols']

            Nframes = 0
            for ic in range(1, Ncols+1):
                Nframes += structure['col%03i' % ic]['frames']

            tasksequence.append((taskname, testkey, Nframes))

        # Launching tasks

        fahrtig = False
        tini = datetime.datetime.now()
        tlastavailable = datetime.datetime.now()

        while not fahrtig:
            
            tmpEL = os.path.join(dayfolder, 'EXP_LOG_*.txt')
            explogfs = glob.glob(tmpEL)
            explogfs = pilib.sortbydateexplogfs(explogfs)
            
            
            #t1 = datetime.datetime.now()
            explog = pilib.loadexplogs(explogfs, elvis)
            #t2 = datetime.datetime.now()
            #print '%.1f seconds in loading explogs...' % (t2-t1).seconds
            
                                                          
            if self.startobsid > 0:
                ixstart = np.where(explog['ObsID'] == self.startobsid)[0][0]
                explog = explog[ixstart:].copy()
            

            for it, taskitem in enumerate(tasksequence):

                taskname, testkey, Nframes = taskitem
                available = pilib.coarsefindTestinExpLog(
                    explog, testkey, Nframes)                

                if available:
                    # print '%s available, doing nothing!' % taskname # TESTS
                    self.inputs[taskname]['explogf'] = explogfs
                    if self.startobsid > 0:
                        self.inputs[taskname]['OBSID_lims'] = [
                            self.startobsid, explog['ObsID'][-1]]
                    self.launchtask(taskname)
                    tasksequence.pop(it)
                    
                    tlastavailable = datetime.datetime.now()

                    partdtm = ((tlastavailable-tini).seconds)/60.
            
                    partsummary = self.get_execution_summary(exectime=partdtm)
                    if self.log is not None:
                        self.log.info(partsummary)

            if len(tasksequence) == 0:
                fahrtig = True
            else:
                if self.log is not None:
                    self.log.info(
                        'Pipeline sleeping for %i seconds...' % waittime)
                sleep(waittime)
                justnow = datetime.datetime.now()
                tsincelastavailable = ((justnow-tlastavailable).seconds)/60.
                
                if tsincelastavailable > waitTO:
                    
                    print " Im bored of waiting... do you want to give up? y/n"
                    rlist, _, _ = select([sys.stdin], [], [], 60.) # 1 minute to answer

                    if rlist:
                        ans = sys.stdin.readline().lower()
                    else:
                        ans = 'y'
                        print "No input. Assuming that's a 'y' and hence quitting."
                    
                    
                    if ans == 'y':
                        if self.log is not None:
                            self.log.info('%.1f hours since last test was available... abandoning' % tsincelastavailable/3600.)
                        fahrtig = True
                    

        tend = datetime.datetime.now()
        Dtm = ((tend-tini).seconds)/60.

        summary = self.get_execution_summary(exectime=Dtm)
        if self.log is not None:
            self.log.info(summary)

        return None

    def dotask(self, taskname, inputs, drill=False, debug=False, cleanafter=False):
        """Generic test master function."""
        
        strip_taskname = utils.remove_iter_tag(taskname,Full=False)

        tini = datetime.datetime.now()

        Errors = False
        
        try:
            test = self.get_test(strip_taskname, inputs=inputs, log=self.log, 
                                 drill=drill, debug=debug, cleanafter=cleanafter)
            
            Errors = test()  # test execution
        except:
            self.catchtraceback()
            Errors = True
            if self.log is not None:
                self.log.info('TASK "%s@%s" FAILED, QUITTING!' %
                              (taskname, self.Test_dict[strip_taskname].__module__))
            else:
                print 'TASK "%s@%s" FAILED, QUITTING!' % (
                    taskname, self.Test_dict[taskname].__module__)
        
        tend = datetime.datetime.now()
        dtm = ((tend-tini).seconds)/60.
        
        execlog = OrderedDict()
        execlog['Task'] = taskname
        execlog['Errors'] = Errors
        execlog['exectime'] = dtm
        try:
            execlog['ObsID_range'] = (test.dd.mx['ObsID'][:].min(),
                                   test.dd.mx['ObsID'][:].max() )
        except:
            execlog['ObsID_range'] = (-1,-1)
        
        try:
            flags = test.dd.flags.getFlagsOnList()
        except:
            flags = ['UNKNOWN']
        execlog['flags'] = flags.__repr__()
        
        test = None
        
        return execlog
    