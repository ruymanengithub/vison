# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PSF0X

PSF vs. Fluence, and Wavelength
   PSF01 - nominal temperature
   PSF02 - alternative temperatures

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
import os
from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.datamodel import scriptic as sc
from vison.flat import FlatFielding as FFing
from vison.point import lib as polib
from vison.support.report import Report
from vison.support import files
from copy import deepcopy
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

PSF0X_commvalues = dict(program='CALCAMP',
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,exptime=0.,shuttr=1,
  siflsh=1,siflush_p=500,
  motr_on=1,
  motr_cnt=2,
  source='point',
  mirr_on=1,
  wave=4,mirr_pos=polib.mirror_nom['F4'],
  comments='')

testdefaults = dict(PSF01=dict(waves=[590,640,800,880],
                               exptimes=dict(),
                               frames=[20,15,10,4,3]))

for w in testdefaults['PSF01']['waves']:
    testdefaults['PSF01']['exptimes']['nm%i' % w] = np.array([5.,25.,50.,75.,90.])/100.*ogse.tFWC_point['nm%i' % w]



def build_PSF0X_scriptdict(exptimes,frames,wavelength=800,
        diffvalues=dict(),elvis='6.3.0'):
    """ 
    
    Builds PSF0X script structure dictionary.
    
    :param exptimes: list of ints, [ms], exposure times.
    :param frames: list of frame numbers. Same length as exptimes.
    :param wavelength: int, [nm], wavelength.
    :param diffvalues: dict, opt, differential values.
    
    
    """
    
    assert len(exptimes) == len(frames)
    
    FW_ID = ogse.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    
    PSF0X_commvalues['wave'] = FW_IDX
    PSF0X_commvalues['mirr_pos'] = polib.mirror_nom['F%i' % FW_IDX]
    
    ncols = len(exptimes)
    
    PSF0X_sdict = dict()
    
    for ic in range(ncols):
        colid = 'col%i' % (ic+1,)
    
        PSF0X_sdict[colid]=dict(frames=frames[ic],exptime=exptimes[ic])

    Ncols = len(PSF0X_sdict.keys())    
    PSF0X_sdict['Ncols'] = Ncols
    
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(PSF0X_commvalues)               
               
    PSF0X_sdict = sc.update_structdict(PSF0X_sdict,commvalues,diffvalues)
    
    
    return PSF0X_sdict


def filterexposures(structure,explogf,datapath,OBSID_lims,elvis='6.3.0'):
    """Loads a list of Exposure Logs and selects exposures from test PSF0X.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    
    """
    
    # load exposure log(s)
    explog = pilib.loadexplogs(explogf,elvis=elvis,addpedigree=True,
                               datapath=datapath)
    Ncols = structure['Ncols']
    
    Filters = [structure['col%i' % i]['wave'] for i in range(1,Ncols+1)]
    Filter = Filters[0]
    assert np.all(np.array(Filters) == Filter)
    
    
    testkey = structure['col1']['test']
    
    selbool = (explog['test'] == testkey) & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) & \
        (explog['wave'] == Filter) # TESTS
        
    explog = explog[selbool]

    # Assess structure
        
    checkreport = pilib.check_test_structure(explog,structure,CCDs=[1,2,3],
                                           wavedkeys=[])
    
    # Labeling of exposures
    explog['label'] = np.array(['None']*len(explog))
    
    frcounter = 0
    for ic in range(1,Ncols+1):
        _frames = structure['col%i' % ic]['frames']
        #print frcounter,frcounter+_frames*3
        explog['label'][frcounter:frcounter+_frames*3] = 'col%i' % ic
        frcounter += _frames*3
    
    
    return explog, checkreport
    

def check_data(dd,report,inputs,log=None):
    """ 
    
    
    
    """
    
    raise NotImplementedError
    


def prep_data(dd,report,inputs,log=None):
    """Takes Raw Data and prepares it for further analysis. Also checks that
    data quality is enough."""
    
    raise NotImplementedError
    
    # Inputs un-packing
    
    FFs = inputs['FFs'] # flat-fields for each CCD (1,2,3)
    
    # Load Flat-Field data for each CCD
    
    FFdata = dict()
    for CCDindex in range(1,4):
        CCDkey = 'CCD%i' % CCDindex
        if CCDkey in FFs.keys():
            FFdata[CCDkey] = FFing.FlatField(fitsfile=FFs[CCDkey])
        
    RepDict['prepPSF0X'] = dict(Title='Data Pre-Processing and QA',items=[])
    
    
    # Loop over CCDs
    
    for CCDindex in range(1,4):
        
        CCDkey = 'CCD%i' % CCDindex
        
        if CCDkey in DataDict:
            
            # Loop over ObsIDs
            
            ObsIDs = DataDict[CCDkey]['ObsID'].copy()
            label = DataDict[CCDkey]['label'].copy()
            Nobs = len(ObsID)
            
            
            for iObs, ObsID in enumerate(ObsIDs):
                
                ilabel = label[iObs]
                
                stop()
                
                # log object being analyzed: ObsID
                
                # retrieve FITS file and open it
                
                # subtract offset and add to DataDict
                # measure STD in pre and over-scan and add to DataDict
                
                # Divide by flat-field
                
                for spotID in spotIDs[CCDkey]:
                    
                    if log is not None: log.info('ObsID - spotID = %s-%s' % (ObsID,spotID))
                    
                    # get coordinates of spotID
                    
                    # Cut-out stamp of the spot
                    
                    # Measure background locally and add to DataDict
                    
                    # if ilabel != fluence_0:
                    #   do basic measurements on each spot and add to DataDict
                    #     peak fluence, peak position (CCD coordinates), FWHM
                    #     quadrupole moments, ellipticity, R2
                    #
                
                    # save spot-data to a hard-file and add path to DataDict
    
            # Data Quality Assessment:
            
            # plot peak fluence, fwhm, ellipticity, vs. Exposure time
            #    save plots as hard-files, keep file name in RepDict
            #
            # Check all parameters are within expected ranges:
            #    offsets, STDs, peak fluence, fwhm, ellipticity
            #    save reports to RepDict
    
       
    return DataDict, RepDict


def basic_analysis(dd,report,inputs,log=None):
    """Performs basic analysis on spots:
         - Gaussian Fits: peak, position, width_x, width_y
    """
    
    raise NotImplementedError
    
    return dd,report

    
def bayes_analysis(dd,report,inputs,log=None):
    """ 
    Performs bayesian decomposition of the spot images:
        - optomechanic PSF and detector PSF.
    Also measures R2, fwhm and ellipticity of "extracted" detector PSF.
    """
    
    raise NotImplementedError
    
    return dd,report

def meta_analysis(dd,report,inputs,log=None):
    """
    
    Analyzes the relation between detector PSF and fluence.
    
    """
    
    raise NotImplementedError
    
    return dd, report
    

def feeder(inputs,elvis='6.3.0'):
    """ """
    
    subtasks = [('check',check_data),('prep',prep_data),
                ('basic',basic_analysis),('bayes',bayes_analysis),
                ('meta',meta_analysis)]
    
    wavelength = inputs['wavelength']
    
    
    testkey = inputs['test']
    if 'PSF01' in testkey: _testkey = 'PSF01'
    
    if 'exptimes' in inputs:
        exptimes = inputs['exptimes']
    else:
        exptimes = testdefaults[_testkey]['exptimes']['nm%i' % wavelength]
    if 'frames' in inputs:
        frames = inputs['frames']
    else:
        frames = testdefaults[_testkey]['frames']
        
    if 'elvis' in inputs:
        elvis = inputs['elvis']
    if 'diffvalues' in inputs:
        diffvalues = inputs['diffvalues']
    else:
        diffvalues = {}
    
    diffvalues['test'] = testkey    
    scriptdict = build_PSF0X_scriptdict(exptimes,frames,wavelength,
                                        diffvalues,elvis=elvis)
    
    inputs['structure'] = scriptdict
    inputs['subtasks'] = subtasks
    
    return inputs
