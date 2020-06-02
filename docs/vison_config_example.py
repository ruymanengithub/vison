#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

EUCLID VIS CALIBRATION CAMPAIGN
"vison" CONFIGURATION FILE

BLOCKID: EINSTEIN   # NEEDS USER INPUTS.
Filled in by: RAF   # NEEDS USER INPUTS.

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import sys
import os
from vison.campaign import campaign, devel
import copy
# END IMPORT


# ALL to_do dictionaries:
#
# MOT_WARM: [None,None],dict(init=True,check=True,basic=True)
# COSMETICS00: [None,None],dict(init=True,check=True,masks=True,meta=True)
# ...
                 
def get_config_dict(testgenerator, datapath, resroot, 
                    BLOCKID,CHAMBER,diffvalues,
                    inCDPs,perflimits,tasks2execute=[],
                    doReport=True,elvis='7.5.X'):
    """Produces input dictionaries to run the pipeline on a series of Tasks/Tests."""
    
    return inputdict



def add_RUN_specifics(inputdict,RUN,dataroot=''):
    """Adds specifics to the input dictionaries of each Task/Test. This is where
    we taylor execution for each test (e.g. choose what to do within the Task, 
    or tell the pipeline where to find the data and what OBSIDs to consider for 
    each Test).
    
    """
    import os
    import copy
    import glob
    from vison.pipe.lib import sortbydateexplogfs
    
    # How to fill-in the specifics for each test within a RUN:
    # Code for the test spec lists:
    #    ['TEST-NAME', [START-OBSID(integer), END-OBSID(integer)],todo(dictionary)]
    #    All to-do dictionaries provided at the beginning of this script.
    #    If you want to by-pass a test within a RUN, comment the line with a 'hash'.
    
    # NEEDS USER INPUTS: Tests lists, ObsID limits, sub-tasts todo booleans
    
    test_specifics = dict(
    
    D00 = dict(
            
            schedule=[
                    ['MOT_WARM',[30332, 30339],dict(init=True,check=True,basic=True)]
                    ],
            datapath=os.path.join(dataroot,'25_Jul_19')
        )
    )
    
    # END USER INPUTS
    
    # ...
    
    
    return inputdict

if __name__ == '__main__':
    
    # MASTER
    
    # NEEDS USER INPUTS:
    dataroot = '../atCALDATA/data'
    datapath = os.path.join(dataroot,'NN_Mmm_19')
    cdppath = 'calproducts'
    resroot = 'results_atCALDATA/'
    BLOCKID = 'EINSTEIN'
    CHAMBER = 'B_EINSTEIN' # example: 'A_MAX'
    
    diffvalues = dict(operator = 'unk',
            sn_ccd1 = '14471-19-01', # example: '14183-18-01',
            sn_ccd2 = '15081-15-02', # example: '14311-04-01',
            sn_ccd3 = '14471-10-02', # example: '14173-14-02',
            sn_roe= 'FM14', # example: FM01
            sn_rpsu = 'FM14') # example: FM02

    inCDPs = dict(Mask=dict(
            CCD1=os.path.join(cdppath,'masks/EUC_MASK_%s_CCD1_SN_%s.fits' % \
                              (BLOCKID,diffvalues['sn_ccd1'])),
            CCD2=os.path.join(cdppath,'masks/EUC_MASK_%s_CCD2_SN_%s.fits' % \
                              (BLOCKID,diffvalues['sn_ccd2'])),
            CCD3=os.path.join(cdppath,'masks/EUC_MASK_%s_CCD3_SN_%s.fits' % \
                              (BLOCKID,diffvalues['sn_ccd3']))),
        Gain=dict(
            nm730=os.path.join(cdppath, 'gain/nm730/PTC02_730_GAIN_TB.pick')
            ),
        FF=dict(
            nm730=dict(
                CCD1=os.path.join(cdppath,'flats/nm730/EUC_FF_730nm_col001_ROE1_CCD1.fits'),
                CCD2=os.path.join(cdppath,'flats/nm730/EUC_FF_730nm_col001_ROE1_CCD2.fits'),
                CCD3=os.path.join(cdppath,'flats/nm730/EUC_FF_730nm_col001_ROE1_CCD3.fits')))
        )

    perflimits = dict()
    doReport = True
    elvis='7.5.X'
    # END USER INPUTS.
    
    # NEEDS USER INPUTS.
    
    tasks2execute_dict=dict(
                            D00=['MOT_WARM'],
                            D11=['COSMETICS00', 'FOCUS00', 'PSFLUX00', 'FLATFLUX00', 'FLAT_STB', 'PTC02WAVE'],
                            D12=['BIAS02','BIAS01','CHINJ01','CHINJ02','TP11','TP21',
                                 'FLATFLUX00', 'PTC01','NL02','FLAT01','DARK01','PTC02WAVE',
                                 'PERSIST01','PSFLUX00'],
                            D12E=['BF01','BF01WAVE'],
                            D21=['BIAS02','PSF01','PTC02WAVE'],
                            D21E=['BF01WAVE'],
                            D22=['BIAS02','PSF01','FLAT02'],
                            DTEST=['BF01WAVE']#'BF01WAVE']
                        )
    # END USER INPUTS.
    
    testgenerator = campaign.generate_test_sequence
    
    inputdict = dict()
    
    if np.any([RUN == '' for RUN in RUNs]):
        sys.exit('Hey, what RUN do you want to process?')
    
    for RUN in RUNs:
        
        tasks2execute = tasks2execute_dict[RUN]
        
            # ...
            
            
            inputdict.update(_inputdict)
            

