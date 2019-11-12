#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Pipelining for FPA analysis.



Created on Wed Oct  2 16:15:56 2019

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


from vison.pipe import master
from vison import __version__
from vison.support import logger as lg
#from vison.support.report import Report
from vison.support import vistime
#from lib import get_time_tag
#from vison.pipe import lib as pilib
#from vison.support import context
#from vison.support import utils
#from vison.support import flags as flagsmodule
# END IMPORT

isthere = os.path.exists

defaults = dict()


class FpaPipe(master.GenPipe):
    """Master Class of FM-analysis at block-level of assembly."""

    from vison.fpatests.cea_dec19.FWD_WARM import FWD_WARM

    Test_dict = dict(FWD_WARM=FWD_WARM)

    def __init__(self, inputdict, dolog=True, drill=False, debug=False, startobsid=0,
                 processes=1, tag='', cleanafter=False):
        """ """

        self.inputs = defaults.copy()
        self.inputs.update(inputdict)
        self.tasks = self.inputs['tasks']
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
            self.ID = '%s%s' % (vistime.get_time_tag(), self.tag)  # ID of the analysis "session"

        self.inputs['ID'] = self.ID

        if dolog:
            self.logf = 'FPA_%s.log' % self.ID

            if os.path.exists(self.logf):
                os.system('rm %s' % self.logf)

            self.log = lg.setUpLogger(self.logf)
            self.log.info(['_', 'Starting FPA-Analysis Pipeline'] +
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
