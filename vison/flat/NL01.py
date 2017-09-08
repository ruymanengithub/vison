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
from vison.pipe import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
from copy import deepcopy
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
                          diffvalues=dict(),elvis='6.0.0'):
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



    
def run(inputs,log=None):
    """Test NL01 master function."""
    
