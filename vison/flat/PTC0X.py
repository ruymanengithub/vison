# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PTC_0X

Photon-Transfer-Curve Analysis
   PTC01 - nominal temperature
   PTC02 - alternative temperatures

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract pairs of exposures with equal fluence
    - Synoptic analysis:
        variance vs. fluence
        variance(binned difference-frames) vs. fluence
    - extract: RON, gain, gain(fluence)
    - produce synoptic figures
    - Save results.



Created on Mon Apr  3 17:00:24 2017

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

HKKeys_PTC0X = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

PTC0X_structure = dict(col1=dict(N=5,Exptime=0),
                          col2=dict(N=2,Exptime=1.),
                          col3=dict(N=2,Exptime=5.),
                          col4=dict(N=2,Exptime=10.),
                          col5=dict(N=2,Exptime=15.),
                          col6=dict(N=2,Exptime=18.),
                   Ncols=6)


def filterexposures_PTC0X(inwavelength,explogf,datapath,OBSID_lims,structure=PTC0X_structure,elvis='5.7.04'):
    """Loads a list of Exposure Logs and selects exposures from test PSF0X.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    
    """


def prep_data_PTC0X(DataDict,RepDict,inputs,log=None):
    """Takes Raw Data and prepares it for further analysis. 
    Also checks that data quality is enough."""
    
def basic_analysis_PTC0X(DataDict,RepDict,inputs,log=None):
    """Performs basic analysis of images:
         - builds PTC curves: both on non-binned and binned images
    """

def meta_analysis_PSF0X(DataDict,RepDict,inputs,log=None):
    """
    
    Analyzes the variance and fluence:
       gain, and gain(fluence)
    
    """

def save_progress(DataDict,reportobj,DataDictFile,reportobjFile):        
    files.cPickleDumpDictionary(DataDict,DataDictFile)
    files.cPickleDump(reportobj,reportobjFile)


def recover_progress(DataDictFile,reportobjFile):
    DataDict = files.cPickleRead(DataDictFile)
    reportobj = files.cPickleRead(reportobjFile)
    return DataDict,reportobj
    
def run(inputs,log=None):
    """Test PTC0X master function."""
    
