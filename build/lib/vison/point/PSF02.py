# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PSF02

PSF vs. Fluence, and Wavelength

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract offset level.
    - Divide by Flat-field.
    - Crop stamps of the sources on each CCD/Quadrant.
       - save snapshot figures of sources.
    - for each source:
       - measure shape using weighted moments
       - measure shape using Gaussian Fit
       - Bayesian Forward Modelling the optomechanic+detector PSF
    - Produce synoptic figures.
    - Save results.

Created on Thu Dec 29 15:01:07 2016

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
# END IMPORT

HKKeys_PSF02 = []

def filterexposures_PSF02(explogfs,datapaths):
    """Loads a list of Exposure Logs and selects exposures from test PSF02."""
    
    


def run(inputs,log=None):
    """Test PSF02"""
    
    OBSID_lims = inputs['OBSID_lims']
    explogfs = inputs['explogfs']
    datapaths = inputs['datapaths']
    resultspath = inputs['resultspath']
    
    
    DataDict = filterexposures_PSF02(explogfs,datapaths,OBSID_lims)
    
    DataDict = addHK(DataDict,HKKeys_PSF02)
    
    
    