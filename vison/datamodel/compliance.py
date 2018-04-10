#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Some functions to present COMPLIANCE MATRICES.

Created on Mon Apr  9 17:32:03 2018

:author: raf
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import pandas as pd
import collections
import string as st
# END IMPORT

def countdictlayers(d):
    return max(countdictlayers(v) if isinstance(v,dict) else 0 for v in d.values()) + 1


def convert_compl_to_nesteditemlist(complidict):
    """ """
    
    def traverse_tree(dictionary,nesteditemlist):   
            Nkeys = len(dictionary.keys())
            if Nkeys >5:
                nesteditemlist +=[
                '\item %s' % dictionary.__repr__()
                        ]
                return nesteditemlist
            
            for key,value in dictionary.items():
                if isinstance(value,(dict,OrderedDict)):
                    nesteditemlist += [
                    '\item %s:' % key, 
                    '\\begin{itemize}']
                    nesteditemlist = traverse_tree(value,nesteditemlist)
                    nesteditemlist += ['\\end{itemize}']
                else:
                    nesteditemlist += [
                    '\item %s : %s' % (key,value)        
                            ]
            return nesteditemlist
    
    tex =  ['\\begingroup'] +\
           ['\\scriptsize'] +\
           ['\\begin{itemize}']+\
           traverse_tree(complidict,[])+\
           ['\\end{itemize}']
           
    return tex


def removescalars_from_dict(indict):
    """ """
    
    def traverse_tree(indict):
        for key,value in indict.iteritems():
            if isinstance(value,(dict,OrderedDict)):
                indict[key] = traverse_tree(value)
            elif not isinstance(value,collections.Sequence):
                indict[key] = [value]
        return indict
        
    outdict = traverse_tree(indict)
        
    return outdict


def gen_compliance_tex(indict):
    """ """
    
    complidict = removescalars_from_dict(indict)
    
    level = countdictlayers(complidict)
    
    if level < 3:
        df = pd.DataFrame.from_dict(complidict)
        tex = df.to_latex(multicolumn=True,multirow=True,longtable=True,index=False)        
    elif level == 3:
        keys = []
        frames = []
        for key,d in complidict.iteritems():
            keys.append(key)
            frames.append(pd.DataFrame.from_dict(d))
        df = pd.concat(frames,keys=keys)
        tex = df.to_latex(multicolumn=True,multirow=True,longtable=True,index=False)
    else:
        tex = convert_compl_to_nesteditemlist(complidict)
    tex = st.split(tex,'\n')
    return tex
    