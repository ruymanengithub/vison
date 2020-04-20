#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:07:21 2019

@author: raf
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
# END IMPORT

defaults = dict(BLOCKID='FM',
                CHAMBER='Unk')


class GenPipe(object):
    """Abstract Master Class of FM-analysis, any level of assembly."""

    Test_dict = dict()

    def __init__(self, inputdict, dolog=True, drill=False, debug=False, startobsid=0,
                 processes=1, tag='', cleanafter=False):
        """ """

        self.inputs = defaults.copy()
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
        self.BLOCKID = self.inputs['BLOCKID']
        self.CHAMBER = self.inputs['CHAMBER']
        self.drill = drill
        self.debug = debug
        self.processes = processes
        self.cleanafter = cleanafter
        self.tag = tag
        self.completion = OrderedDict()

        if self.debug:
            self.ID = 'PipeDebug%s' % self.tag
        else:
            self.ID = 'FM%s%s' % (vistime.get_time_tag(), self.tag)  # ID of the analysis "session"

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

    def _get_log_header(self):
        log_header = [
            'Pipeline ID: %s' % self.ID,
            'BLOCK ID: %s' % self.BLOCKID,
            'Chamber: %s\n' % self.CHAMBER,
            'vison version: %s\n' % __version__,
            'Tasks: %s\n' % self.tasks.__repr__()]
        return log_header

    def get_test(
            self,
            taskname,
            inputs=dict(),
            log=None,
            drill=False,
            debug=False,
            cleanafter=False):
        """ """
        testclass = self.Test_dict[taskname]
        #pathtotest = os.path.split(inspect.getfile(testclass))[0]
        # os.system('rm %s/*.pyc' % pathtotest) # avoids problems with .pyc files
        # from previous tasks in the run... DIRTY HACK
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

        taskinputs.update(dict(ID=self.ID,
                               BLOCKID=self.BLOCKID,
                               CHAMBER=self.CHAMBER,
                               processes=self.processes))

        msg = ['\n\nRunning Task: %s\n' % taskname]
        msg += ['Inputs:\n']
        for key in taskinputs:
            msg += ['%s = %s' % (key, str(taskinputs[key]))]

        if self.log is not None:
            self.log.info(msg)

        print(('Running: %s' % taskname))

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
            partdtm = ((tpart - tini).seconds) / 60.

            partsummary = self.get_execution_summary(exectime=partdtm)
            if self.log is not None:
                self.log.info(partsummary)

        tend = datetime.datetime.now()
        Dtm = ((tend - tini).seconds) / 60.

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
            for ic in range(1, Ncols + 1):
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

                    partdtm = ((tlastavailable - tini).seconds) / 60.

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
                tsincelastavailable = ((justnow - tlastavailable).seconds) / 60.

                if tsincelastavailable > waitTO:

                    print(" Im bored of waiting... do you want to give up? y/n")
                    rlist, _, _ = select([sys.stdin], [], [], 60.)  # 1 minute to answer

                    if rlist:
                        ans = sys.stdin.readline().lower()
                    else:
                        ans = 'y'
                        print("No input. Assuming that's a 'y' and hence quitting.")

                    if ans == 'y':
                        if self.log is not None:
                            self.log.info(
                                '%.1f hours since last test was available... abandoning' %
                                tsincelastavailable / 3600.)
                        fahrtig = True

        tend = datetime.datetime.now()
        Dtm = ((tend - tini).seconds) / 60.

        summary = self.get_execution_summary(exectime=Dtm)
        if self.log is not None:
            self.log.info(summary)

        return None

    def dotask(self, taskname, inputs, drill=False, debug=False, cleanafter=False):
        """Generic test master function."""

        strip_taskname = utils.remove_iter_tag(taskname, Full=False)

        tini = datetime.datetime.now()

        Errors = False

        try:
            test = self.get_test(strip_taskname, inputs=inputs, log=self.log,
                                 drill=drill, debug=debug, cleanafter=cleanafter)

            Errors = test()  # test execution
        except BaseException:
            self.catchtraceback()
            Errors = True
            if self.log is not None:
                self.log.info('TASK "%s@%s" FAILED, QUITTING!' %
                              (taskname, self.Test_dict[strip_taskname].__module__))
            else:
                print(('TASK "%s@%s" FAILED, QUITTING!' % (
                    taskname, self.Test_dict[taskname].__module__)))

        tend = datetime.datetime.now()
        dtm = ((tend - tini).seconds) / 60.

        execlog = OrderedDict()
        execlog['Task'] = taskname
        execlog['Errors'] = Errors
        execlog['exectime'] = dtm
        try:
            execlog['ObsID_range'] = (test.dd.mx['ObsID'][:].min(),
                                      test.dd.mx['ObsID'][:].max())
        except BaseException:
            execlog['ObsID_range'] = (-1, -1)

        try:
            flags = test.dd.flags.getFlagsOnList()
        except BaseException:
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
            print(msg_trbk)
