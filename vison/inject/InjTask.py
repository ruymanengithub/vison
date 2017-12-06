#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:56:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from vison.pipe.task import Task
# END IMPORT

class InjTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(InjTask,self).__init__(*args,**kwargs)