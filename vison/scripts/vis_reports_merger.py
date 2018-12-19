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
from astropy.io import ascii
import datetime
import numpy as np
from optparse import OptionParser
import string as st

from vison.support.report import Report
from vison.support.files import cPickleRead

# END IMPORT

def run_merger(infile, programme='FM', block='BLOCK', reference='6_666'):
    """ """
    
    indata = ascii.read(infile)
    
    report = Report(TestName='All Calibration Tests',Model=programme,
                    Reference=reference)
    
    reportroot = 'EUCL_MSS_TR_%s_%s' % (reference,block)
    
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
    
    

if __name__ == '__main__':

    parser = OptionParser()
    
    parser.add_option("-i", "--infile", dest="infile",
                      default='', help="Inputs file with list of reports to be merged.")
    parser.add_option("-p", "--programme", dest="programme", 
                      default="FM", help="Calibration programme (QM/FQM/FM).")
    parser.add_option("-b", "--block", dest="block",
                      default="BLOCK", help="Name of the block being calibrated.")
    parser.add_options("-r", "--reference", dest="reference",
                       default="6_666", help="")

    (options, args) = parser.parse_args()

    infile = options.infile
    programme = options.programme
    block = options.block
    reference = options.reference
    

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

    run_merger(infile, programme, block, reference)
    
