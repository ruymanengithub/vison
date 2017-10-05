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
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=1,chinj=0,tpump=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 6',pos_cal_mirror=polib.mirror_nom['Filter4'],
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
    
    NL01_commvalues['wavelength'] = 'Filter %i' % FW_IDX
    
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
    METACODE
    
    Checks that data quality is good enough.
    
    # check common HK values are within safe / nominal margins
    # check voltages in HK match commanded voltages, within margins
    
    # f.e.ObsID, f.e.CCD, f.e.Q.:
        # measure offsets/means in pre-, img-, over-
        # measure std in pre-, img-, over-
    # assess std in pre- is within allocated margins
    # (assess offsets in pre- and over- are equal, within allocated  margins)
    # assess image-fluences are within allocated margins for each exposure time
    
    # plot fluences vs. exposure time
    # plot std-pre vs. time
    
    # issue any warnings to log
    # issue update to report    
    
    
    """

def prep_data_NL01(DataDict,RepDict,inputs,log=None):
    """
    
    METACODE
    
    Takes Raw Data and prepares it for further analysis. 
        
    # f.e. ObsID, f.e.CCD, f.e.Q:
        # subtract offset
        # opt: [sub bias frame]
           
    
    """
    
def extract_stats(DataDict,RepDict,inputs,log=None):
    """
    
    METACODE
    
    Performs basic analysis: extracts statistics from 
    image regions to later build NLC.
    
    # create segmentation map given grid parameters    
    
    
    # f.e. ObsID, f.e.CCD, f.e.Q:
        
        # f.e. segment:
            # measure central value
            # measure variance
    
    """


def produce_NLCs(DataDict,RepDict,inputs,log=None):
    """ 
    
    METACODE
    
    Obtains Best-Fit Non-Linearity Curve
    
    # f.e. CCD, f.e. Q:
        
        # [opt] apply correction for source variability (interspersed exposure 
        #      with constant exptime)
        # Build NL Curve (NLC) - use stats and exptimes
        # fit poly. shape to NL curve
    
    # plot NL curves for each CCD, Q
    # report max. values of NL (table)
    
    
    """
    
def do_satCTE(DataDict,RepDict,inputs,log=None):
    """
    METACODE 
    
    # select ObsIDs with fluence(exptime) >~ 0.5 FWC
    
    # f.e. ObsID, CCD, Q:
        - measure CTE from amount of charge in over-scan relative to fluence
    
    # f.e. CCD, Q:
        - get curve of CTE vs. fluence
        - measure FWC from curve in ADU
        
    report FWCs in electrons [via gain in inputs] 
    f.e. CCD, Q (table)
    
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

    
