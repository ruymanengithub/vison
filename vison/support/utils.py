#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

General Purpose Utilities

Created on Tue Apr 10 15:18:07 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import string as st
import os
from functools import reduce
import operator
import tempfile
import shutil
# END IMPORT

credentials_path = os.path.join(os.environ['HOME'], 'CREDENTIALS_VISON')


def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)


def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value


def countdictlayers(d):
    return max(countdictlayers(v) if isinstance(v, dict) else 0 for v in d.values()) + 1


def get_function_module(level=1, reference='vison'):
    import inspect

    stack = inspect.stack()
    stackrow = stack[level]
    funcname = stackrow[3]
    mod_path = stackrow[1]
    mod_path_split = st.split(mod_path, sep=os.path.sep)
    try:
        ixref = mod_path_split.index(reference)
    except ValueError:
        ixref = 0
    module = os.path.join(*mod_path_split[ixref:])

    return funcname, module


def get_path_decorator(dpath):
    def fullinpath_adder(path,extension=''): 
        if len(extension)==0:
            return os.path.join(dpath, path)
        else:
            return os.path.join(dpath,'%s.%s' % (path,extension))
    return np.vectorize(fullinpath_adder)


def remove_iter_tag(taskname,Full=False):
    if '.' in taskname:
        _t =  taskname[0:taskname.find('.')]
        iteration = int(taskname[taskname.find('.')+1:])
        if Full:
            return _t, iteration
        else:
            return _t
    else:
        if Full:
            return taskname, None
        else:
            return taskname

def vec_remove_iter_tag(taskslist,Full=False):
    return np.vectorize(remove_iter_tag)(taskslist,Full).tolist()

def safe_open_file(path,prefix=''):
    _, temp_path = tempfile.mkstemp(prefix=prefix)
    shutil.copy2(path, temp_path)
    f = open(temp_path)
    os.system('rm %s' % temp_path)
    return f