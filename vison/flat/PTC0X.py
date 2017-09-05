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
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
# END IMPORT

isthere = os.path.exists

HKKeys_PTC0X = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

#PTC0X_structure = dict(col1=dict(N=5,Exptime=0),
#                          col2=dict(N=2,Exptime=1.),
#                          col3=dict(N=2,Exptime=5.),
#                          col4=dict(N=2,Exptime=10.),
#                          col5=dict(N=2,Exptime=15.),
#                          col6=dict(N=2,Exptime=18.),
#                   Ncols=6)

PTC0X_commvalues = dict(program='CALCAMP',
  IDL=13000,IDH=18000,IG1=5000,IG2=5000,
  OD_1=26000,RD_1=16000,
  OD_2=26000,RD_2=16000,
  OD_3=26000,RD_3=16000,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=0,chinj_rows_on=20,
  chinj_rows_off=20,chinj_repeat=1,id_width=100,
  id_delay=100,tpump=0,ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength = 'Filter 4', pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  


def build_PTC0X_scriptdict(exptimes,frames,wavelength=800,diffvalues=dict()):
    """ 
   
    PTC01:
        11 av. fluences: 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 
        120% x FWC. 10 frames per fluence (except for last 3 fluences, for which it's only 4 frames).
    
    PTC02 (alt. wavelengths):    
        6 av. fluences: 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
        
    """
    
    assert  len(exptimes) == len(frames)
    
    FW_ID = pilib.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    
    if wavelength == 800: subtest = '01'
    else: subtest = '02'
    
    PTC0X_commvalues['test'] = 'PTC%s_%i' % (subtest,wavelength)
    PTC0X_commvalues['wavelength'] = 'Filter %i' % FW_IDX

    PTC0X_sdict = dict()
    
    for ix,ifra in enumerate(frames):
        iexp = exptimes[ix]
        
        colkey= ['col%i' % (ix+1,)]
    
        PTC0X_sdict[colkey] = dict(frames=ifra,exptime=iexp)

    Ncols = len(PTC0X_sdict.keys())    
    PTC0X_sdict['Ncols'] = Ncols
    
    PTC0X_sdict = sc.update_structdict(PTC0X_sdict,PTC0X_commvalues,diffvalues)
    
    return PTC0X_sdict




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

def run(inputs,log=None):
    """Test PTC0X master function."""
    
