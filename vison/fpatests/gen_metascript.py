#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Generic Script to direct the meta-processing of CALCAMP results.

Created on Fri Aug 9 14:52:00 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import os

from vison.support.files import cPickleDumpDictionary, cPickleDump, cPickleRead
from vison.fpatests.bias import MetaBias
# END IMPORT 



def run(MetaClass, testkey, testinputs, SkipLoad=False, SkipParse=False):
    """ """
    
    outpathroot = '%s_FPA' % testkey.upper()
    
    if not os.path.exists(outpathroot):
        os.system('mkdir %s' % outpathroot)
    
    testinputs['outpathroot'] = outpathroot
    
    metaobj = MetaClass(**testinputs)
    
    
    dictiopick = os.path.join(outpathroot,'%s_dictionary.pick' % testkey.lower())
    if SkipLoad:
        metaobj.inventory = cPickleRead(dictiopick)
    else:
        metaobj.load_block_results()
        cPickleDumpDictionary(metaobj.inventory,dictiopick)
        
    metapick = os.path.join(outpathroot,'%s_meta.pick' % testkey.lower())
    
    if SkipParse:
        metaobj = cPickleRead(metapick)    
    else: 
        for testname in metaobj.testnames:
            metaobj.parse_test_results(testname)
        cPickleDump(metaobj, metapick)
    
    metaobj.dump_aggregated_results()
