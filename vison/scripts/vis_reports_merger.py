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

def run_merger(infile):
    """ """
    
    indata = ascii.read(infile)
    
    report = Report(TestName='All',Model='FQM',
                    Reference='6-666')
    
    reportroot = 'All_report'
    
    N2merge = len(indata)
    
    for i in range(N2merge):
        
        ipick = indata.columns[0][i]
        ireport = cPickleRead(ipick)
        
        # TEMPORARY AD-HOC
        
        path, filename = os.path.split(ipick)
        testname = st.replace(filename,'_Report.pick', '')
        path = os.path.normpath(path)
        folders = path.split(os.sep)
        session = folders[1]
        
        report.add_Chapter('%s: %s' % (session,testname))
        report.Contents += ireport.Contents
    
    
    outfiles = report.doreport(reportroot, cleanafter=True, silent=True)
    
    

if __name__ == '__main__':

    parser = OptionParser()
    
    parser.add_option("-i", "--infile", dest="infile",
                      default='', help="Inputs file with list of reports to be merged.")

    (options, args) = parser.parse_args()

    infile = options.infile

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

    run_merger(infile)
    
