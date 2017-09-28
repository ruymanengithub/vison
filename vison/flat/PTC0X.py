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
#from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from copy import deepcopy
# END IMPORT

isthere = os.path.exists

HKKeys_PTC0X = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS



PTC0X_commvalues = dict(program='CALCAMP',
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=1,chinj=0,tpump=0,motor=0,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength = 'Filter 4', pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='')
  


def build_PTC0X_scriptdict(exptimes,frames,wavelength=800,diffvalues=dict(),
                           elvis='6.0.0'):
    """Builds PTC0X script structure dictionary.
    
    :param exptimes: list of ints [ms], exposure times.
    :param frames: list of ints, number of frames for each exposure time.
    :param wavelength: int, wavelength. Default: 800 nm.
    :param diffvalues: dict, opt, differential values.   
        
    """
    
    assert  len(exptimes) == len(frames)
    
    FW_ID = ogse.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    
    if wavelength == 800: subtest = '01'
    else: subtest = '02'
    
    PTC0X_commvalues['test'] = 'PTC%s_%i' % (subtest,wavelength)
    PTC0X_commvalues['wavelength'] = 'Filter %i' % FW_IDX

    PTC0X_sdict = dict()
    
    for ix,ifra in enumerate(frames):
        iexp = exptimes[ix]
        
        colkey= 'col%i' % (ix+1,)
    
        PTC0X_sdict[colkey] = dict(frames=ifra,exptime=iexp)

    Ncols = len(PTC0X_sdict.keys())    
    PTC0X_sdict['Ncols'] = Ncols
    
               
               
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(PTC0X_commvalues)
    
    PTC0X_sdict = sc.update_structdict(PTC0X_sdict,commvalues,diffvalues)
    
    return PTC0X_sdict




def filterexposures_PTC0X():
    """Loads a list of Exposure Logs and selects exposures from test PSF0X.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    
    """


def check_data(DataDict,RepDict,inputs,log=None):
    """
    
    Checks quality of ingested data.
    
    **METACODE**
    
    ::
        
        check common HK values are within safe / nominal margins
        check voltages in HK match commanded voltages, within margins
    
        f.e.ObsID:
            f.e.CCD:
                f.e.Q.:
                    measure offsets/means in pre-, img-, over-
                    measure std in pre-, img-, over-
        assess std in pre- is within allocated margins
        assess offsets in pre- and over- are equal, within allocated  margins
        assess image-fluences are within allocated margins
    
        plot fluences vs. exposure time
        plot std-pre vs. time
    
        issue any warnings to log
        issue update to report
    
    """
    


def extract_PTC(DataDict,RepDict,inputs,log=None):
    """
    
    Performs basic analysis of images:
        - builds PTC curves: both on non-binned and binned images
    
    **METACODE**
    
    ::
    
        create list of OBSID pairs
    
        create segmentation map given grid parameters
    
        f.e. OBSID pair:
            CCD:
                Q:
                    [apply defects mask if available]
                    subtract CCD images
                    f.e. segment:
                        measure central value
                        measure variance
    
    """



def meta_analysis(DataDict,RepDict,inputs,log=None):
    """

    Analyzes the variance and fluence:
    gain, and gain(fluence)

    METACODE
    
    ::
    
        f.e. CCD:
            Q:
                (using stats across segments:)
                fit PTC to quadratic model
                solve for gain
                solve for alpha (pixel-correls, Guyonnet+15)
                solve for blooming limit
    
        plot PTC curves with best-fit f.e. CCD, Q
        report on gain estimates f. e. CCD, Q (table)
        report on blooming limits (table)
    
    """


def feeder(inputs,elvis='6.1.0'):
    """ """
    
    subtasks = [('check',check_data),('extract',extract_PTC),
                ('meta',meta_analysis)]
    
    
    exptimes = inputs['exptimes']
    frames = inputs['frames']
    wavelength = inputs['wavelength']
    if 'elvis' in inputs:
        elvis = inputs['elvis']
    if 'diffvalues' in inputs:
        diffvalues = inputs['diffvalues']
    else:
        diffvalues = {}
        
    scriptdict = build_PTC0X_scriptdict(exptimes,frames,wavelength=wavelength,
                           diffvalues=diffvalues,
                           elvis=elvis)
    
    inputs['structure'] = scriptdict
    inputs['subtasks'] = subtasks
    
    return inputs
    
