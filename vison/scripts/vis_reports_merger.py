#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Reports Merger - Tests

Created on Wed Dec 19 10:02:02 2018

@author: raf

"""

# IMPORT STUFF
import sys
import os
from pdb import set_trace as stop
#from astropy.io import ascii
import datetime
import numpy as np
from optparse import OptionParser
import string as st

from vison.support.report import Report
from vison.support.files import cPickleRead
from vison.support import vjson
# END IMPORT

def run_merger(infile):
    """ """
    
    indata = vjson.load_jsonfile(infile, useyaml=True)
    
    programme = indata['metadata']['programme']
    block = indata['metadata']['block']
    reference = indata['metadata']['reference']
    try: 
        Issue = indata['metadata']['issue']
    except KeyError: 
        Issue = 0.0
    
    #indata = ascii.read(infile,data_start=0)
    niceblock = st.replace(block,'_.','\_')
    report = Report(TestName='Test Reports, %s' % niceblock,Model=programme,
                    Reference=reference, Issue=Issue)
    
    reportroot = 'EUCL_MSS_TR_%s_%s_v%.1f' % (reference,block,issue)
    
    N2merge = len(indata)
    
    for i in range(N2merge):
        
        ipick = indata.columns[0][i]
        ireport = cPickleRead(ipick)
        
        # TEMPORARY AD-HOC
        path, filename = os.path.split(ipick)
        testname = st.replace(filename,'_Report.pick', '')
        testname = st.replace(testname,'_', '\_')
        path = os.path.normpath(path)
        folders = path.split(os.sep)
        session = folders[1]
        
        report.add_Chapter('%s: %s' % (session,testname))
        report.Contents += ireport.Contents
    
    
    _ = report.doreport(reportroot, cleanafter=False, silent=True)
    

def run_merger_plus(infile):
    """ """
    
    indata = vjson.load_jsonfile(infile, useyaml=True)
    
    programme = indata['metadata']['programme']
    block = indata['metadata']['block']
    reference = indata['metadata']['reference']
    try: 
        Issue = indata['metadata']['issue']
    except KeyError: 
        Issue = 0.0
    

    #indata = ascii.read(infile,data_start=0)
    niceblock = st.replace(block,'_.','\_')
    report = Report(TestName='Test Reports, %s' % niceblock,Model=programme,
                    Reference=reference, Issue=Issue)
    
    reportroot = 'EUCL_MSS_TR_%s_%s' % (reference,block)
    
    parts_names = indata['parts_names']
    parts_dict = indata['parts'].copy()
    
    assert len(set(parts_names)^set(parts_dict.keys()))==0, "parts names do not match parts keys!"
    
    for ipart, part in enumerate(parts_names):
        parts_list = parts_dict[part]
        
        report.add_Part(part)
    
        jN2merge = len(parts_list)
        
        for i in range(jN2merge):
            
            ipick = parts_list[i]
            ireport = cPickleRead(ipick)
            
            # TEMPORARY AD-HOC
            path, filename = os.path.split(ipick)
            testname = st.replace(filename,'_Report.pick', '')
            testname = st.replace(testname,'_', '\_')
            path = os.path.normpath(path)
            folders = path.split(os.sep)
            session = folders[1]
                        
            report.add_Chapter('%s: %s' % (session,testname))
            report.Contents += ireport.Contents
    
    
    _ = report.doreport(reportroot, cleanafter=False, silent=True)
    

if __name__ == '__main__':

    parser = OptionParser()
    
    parser.add_option("-i", "--infile", dest="infile",
                      default='', help="Inputs file with list of reports to be merged.")
    #parser.add_option("-p", "--programme", dest="programme", 
    #                  default="FM", help="Calibration programme (QM/FQM/FM).")
    #parser.add_option("-b", "--block", dest="block",
    #                  default="BLOCK", help="Name of the block being calibrated.")
    #parser.add_option("-r", "--reference", dest="reference",
    #                   default="6_666", help="")

    (options, args) = parser.parse_args()

    infile = options.infile
    #programme = options.programme
    #block = options.block
    #reference = options.reference
    

    if infile == '':
        parser.print_help()
        sys.exit()

    header = '\n\n#########################################################\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#       running Reports Merger                          #\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#########################################################\n'

    print header

    #run_merger(infile, programme, block, reference)
    run_merger_plus(infile)
    
