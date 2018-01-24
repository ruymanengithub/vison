#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 17:44:04 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import  numpy as np
from pdb import  set_trace as stop
import copy
import os

from vison.inject import extract_injection_lines
from vison.pipe.task import Task
from vison.datamodel import core,ccd
from vison.pipe import lib as pilib
#from vison.pipe.task import Task
from vison.inject.InjTask import InjTask
# END IMPORT

class PumpTask(InjTask):
    
    def __init__(self,*args,**kwargs):
        super(PumpTask,self).__init__(*args,**kwargs)
        
        
    def check_data(self):
        """ """
        test = self.inputs['test']
        if test == 'TP01':
            kwargs = dict(pattern=(2066,0,1))
        elif test == 'TP02':
            kwargs = dict(pattern=(2066,0,1))
        InjTask.check_data(self,**kwargs)

    
    def check_metrics_ST(self,**kwargs):
        """ 
        
        """        
        super(PumpTask,self).check_metrics_ST(**kwargs)