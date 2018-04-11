#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

General Purpose Utilities

Created on Tue Apr 10 15:18:07 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import string as st
import os

# END IMPORT

def get_function_module(level=1,reference='vison'):
    import inspect
    
    stack = inspect.stack()
    stackrow = stack[level]
    funcname = stackrow[3]
    mod_path = stackrow[1]
    mod_path_split = st.split(mod_path,sep=os.path.sep)
    try: ixref = mod_path_split.index(reference)
    except: ixref = 0
    module = os.path.join(*mod_path_split[ixref:])
    
    return funcname,module