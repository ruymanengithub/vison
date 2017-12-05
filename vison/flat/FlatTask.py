#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:00:10 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF

from vison.pipe.task import Task
# END IMPORT

class FlatTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(FlatTask,self).__init__(*args,**kwargs)
    