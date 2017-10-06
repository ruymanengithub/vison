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
#from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
from copy import deepcopy
# END IMPORT

isthere = os.path.exists

HKKeys_NL01 = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS


NL01_commvalues = dict(program='CALCAMP',
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,exptime=0.,shuttr=1,
  siflsh=1,siflsh_p=500,
  wave=6,
  source='flat',
  comments='')


def build_NL01_scriptdict(expts,exptinter,frames,wavelength=0,
                          diffvalues=dict(),elvis='6.3.0'):
    """Builds NL01 script structure dictionary.
    
    :param expts: list of ints [ms], exposure times.
    :param exptinter: int, ms, exposure time of interleaved source-stability exposures.
    :param frames: list of ints, number of frames for each exposure time.
    :param wavelength: int, wavelength. Default: 0 (Neutral Density Filter)
    :param diffvalues: dict, opt, differential values.
    """

    assert  len(expts) == len(frames)
    
    
    FW_ID = ogse.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    
    NL01_commvalues['wave'] = FW_IDX
    
    NL01_sdict = dict()
    
    NL01_sdict['col1'] = dict(frames=4,exptime=0,comment='BGD')
    NL01_sdict['col2'] = dict(frames=1,exptime=exptinter,comment='STAB')
    
    for ix,ifra in enumerate(frames):

        iexp = expts[ix]
        
        colkeyFlu = 'col%i' % (ix*2+2,)
    
        NL01_sdict[colkeyFlu] = dict(frames=ifra,exptime=iexp,comment='Fluence%i' % (ix+1,))
        
        colkeySta = 'col%i' % (ix*2+2+1,)
    
        NL01_sdict[colkeySta] = dict(frames=1,exptime=exptinter,comment='STA')
    

    Ncols = len(NL01_sdict.keys())    
    NL01_sdict['Ncols'] = Ncols
              
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(NL01_commvalues)
    
    NL01_sdict = sc.update_structdict(NL01_sdict,commvalues,diffvalues)
    
    return NL01_sdict


def filterexposures_NLC01():
    """Loads a list of Exposure Logs and selects exposures from test PSF0X.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    """


def check_data_NL01(DataDict,RepDict,inputs,log=None):
    """
    
    NL01: Checks that data quality is good enough.
    
    **METACODE**
    
    ::
        
        Check common HK values are within safe / nominal margins
        Check voltages in HK match commanded voltages, within margins
    
        f.e.ObsID:
            f.e.CCD:
                f.e.Q.:
                    measure offsets/means in pre-, img-, over-
                    measure std in pre-, img-, over-
        assess std in pre- is within allocated margins
        (assess offsets in pre- and over- are equal, within allocated  margins)
        assess image-fluences are within allocated margins for each exposure time
        
        plot fluences vs. exposure time
        plot std-pre vs. time
    
        issue any warnings to log
        issue update to report
    
    """

def prep_data_NL01(DataDict,RepDict,inputs,log=None):
    """
    
    Takes Raw Data and prepares it for further analysis. 
    
    **METACODE**
    
    ::
        
        f.e. ObsID:
            f.e.CCD:
                f.e.Q:
                    subtract offset
                    opt: [sub bias frame] 
    
    """
    
def extract_stats(DataDict,RepDict,inputs,log=None):
    """
    
    Performs basic analysis: extracts statistics from 
    image regions to later build NLC.
    
    **METACODE**
    
    ::
        
        create segmentation map given grid parameters    
    
        f.e. ObsID:
            f.e.CCD:
                f.e.Q:
                    f.e. "img-segment": (done elsewhere)
                        measure central value
                        measure variance

    """


def produce_NLCs(DataDict,RepDict,inputs,log=None):
    """ 
    
    **METACODE**
    
    ::
    
        Obtains Best-Fit Non-Linearity Curve
    
        f.e. CCD:
            f.e. Q:
        
                [opt] apply correction for source variability (interspersed exposure 
                  with constant exptime)
                Build NL Curve (NLC) - use stats and exptimes
                fit poly. shape to NL curve
    
        plot NL curves for each CCD, Q
        report max. values of NL (table)
    
    """
    
def do_satCTE(DataDict,RepDict,inputs,log=None):
    """
    
    **METACODE**
    
    ::
    
        select ObsIDs with fluence(exptime) >~ 0.5 FWC
        
        f.e. ObsID: 
            CCD: 
                Q:
                    measure CTE from amount of charge in over-scan relative to fluence
    
        f.e. CCD: 
            Q:
                get curve of CTE vs. fluence
                measure FWC from curve in ADU
        
        report FWCs in electrons [via gain in inputs] f.e. CCD, Q (table)
    
    """
    
    
    
def feeder(inputs,elvis='6.3.0'):
    """ """
    
    subtasks = [('check',check_data_NL01),('prep',prep_data_NL01),
                ('extract',extract_stats),
                ('NL',produce_NLCs),
                ('satCTE',do_satCTE)]
    
    
    expts = inputs['exptimes']
    exptinter = inputs['exptinter']
    frames = inputs['frames']
    wavelength = inputs['wavelength']
    if 'elvis' in inputs:
        elvis = inputs['elvis']
    if 'diffvalues' in inputs:
        diffvalues = inputs['diffvalues']
    else:
        diffvalues = {}
        
    
    scriptdict = build_NL01_scriptdict(expts,exptinter,frames,
                        wavelength=wavelength,
                          diffvalues=diffvalues,elvis=elvis)
    

    
    inputs['structure'] = scriptdict
    inputs['subtasks'] = subtasks
    
    return inputs

    
