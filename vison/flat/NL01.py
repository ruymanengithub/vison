# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: NL01

End-To-End Non-Linearity Curve


Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract offset level.
    - Divide by Flat-field.
    - Synoptic analysis:
        fluence ratios vs. extime ratios >> non-linearity curve
    - extract: Non-Linearity curve for each CCD and quadrant
    - produce synoptic figures
    - Save results.


Created on Mon Apr  3 17:38:00 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from vison.pipe import lib as plib
from vison.pipe import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
# END IMPORT

isthere = os.path.exists

HKKeys_NL01 = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

NL01_structure = dict(col1=dict(N=5,Exptime=0),
                          col2=dict(N=2,Exptime=1.),
                          col3=dict(N=2,Exptime=5.),
                          col4=dict(N=2,Exptime=10.),
                          col5=dict(N=2,Exptime=15.),
                          col6=dict(N=2,Exptime=18.),
                   Ncols=6)


def filterexposures_NLC01(inwavelength,explogf,datapath,OBSID_lims,structure=NL01_structure,elvis='5.7.04'):
    """Loads a list of Exposure Logs and selects exposures from test PSF0X.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    """


def prep_data_NL01(DataDict,RepDict,inputs,log=None):
    """Takes Raw Data and prepares it for further analysis. 
    Also checks that data quality is enough."""
    
def basic_analysis_NL01(DataDict,RepDict,inputs,log=None):
    """Performs basic analysis:
         - builds Non-Lin curves
         - fits Non-Linearity curves
    """



def save_progress(DataDict,reportobj,DataDictFile,reportobjFile):        
    files.cPickleDumpDictionary(DataDict,DataDictFile)
    files.cPickleDump(reportobj,reportobjFile)


def recover_progress(DataDictFile,reportobjFile):
    DataDict = files.cPickleRead(DataDictFile)
    reportobj = files.cPickleRead(reportobjFile)
    return DataDict,reportobj
    
def run(inputs,log=None):
    """Test NL01 master function."""
    
