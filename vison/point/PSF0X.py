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

HKKeys_PSF0X = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

PSF0X_structure = dict(col1=dict(N=5,Exptime=0),
                          col2=dict(N=20,Exptime=1.),
                          col3=dict(N=18,Exptime=5.),
                          col4=dict(N=10,Exptime=10.),
                          col5=dict(N=4,Exptime=15.),
                          col6=dict(N=3,Exptime=18.),
                   Ncols=6)

PSF0X_commvalues = dict(program='CALCAMP',
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=1,sinvflushp=500,
  chinj=0,tpump=0,motor=0,
  matrix_size=2,
  add_h_overscan=0,add_v_overscan=0,
  toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  comments='')



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
    
    PSF0X_commvalues['wavelength'] = 'Filter %i' % FW_IDX
    PSF0X_commvalues['pos_cal_mirror'] = polib.mirror_nom['Filter%i' % FW_IDX]
    
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



def filterexposures_PSF0X(inwavelength,explogf,datapath,OBSID_lims,structure=PSF0X_structure,elvis='5.7.04'):
    """Loads a list of Exposure Logs and selects exposures from test PSF0X.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    
    """
    
    # load exposure log(s)
    explog = pilib.loadexplogs(explogf,elvis=elvis,addpedigree=True)
    
    # add datapath(s)
    
    if isinstance(datapath,list):
        longestdatapathname = max([len(item) for item in datapath])
        explog['datapath'] = np.zeros(len(explog),dtype='S%i' % longestdatapathname)
        explognumber=explog['explognumber']
        for idata in range(len(datapath)):
            explog['datapath'][explognumber == idata] = datapath[idata]
    
    #OBSIDs = np.array(explog['ObsID'].copy())
    #Exptime = np.array(explog['Exptime'].copy())
    rootFile_name = explog['File_name'].copy()
    #TEST = np.array(explog['TEST'].copy())
    #Wavelength = np.array(explog['Wavelength'].copy())
    #CCD = np.array(explog['CCD'].copy())
    
    
    DataDict = {}
        
    
    for key in pilib.FW.keys():
        if pilib.FW[key] == '%inm' % inwavelength: Filter = key
    
    Filter = '%i' % inwavelength # TESTS
    
    selbool = (['PSF02' in item for item in explog['TEST']]) & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) & \
        (explog['Wavelength'] == Filter) # TESTS
    
    
    # Assess structure
    
    isconsistent = pilib.check_test_structure(explog,selbool,structure)
    
    
    # Build DataDict
    
    DataDict = dict(meta = dict(inwavelength=inwavelength,structure=structure))
    
    for CCDindex in [1,2,3]:
        
        CCDkey = 'CCD%i' % CCDindex
        
        DataDict[CCDkey] = dict()
        
        CCDselbool = selbool & (explog['CCD'] == CCDkey)
        
        if len(np.where(CCDselbool)[0]) == 0:
            continue
        
        Nobs = len(np.where(CCDselbool)[0])
        
        for key in explog.colnames:
            DataDict[CCDkey][key] = explog[key][CCDselbool]
        
        
        Exptime = DataDict[CCDkey]['Exptime'].copy()
        
        label = np.zeros(Nobs,dtype='40str')
        
        uexptimes = np.sort(np.unique(Exptime))
        
        
        for ixflu,uexptime in enumerate(uexptimes):
            
            ixsel = np.where(Exptime == uexptime)
            
            label[ixsel] = 'fluence_%i' % ixflu
        
        DataDict[CCDkey]['label'] = label.copy()
        
        
        rootFile_name = DataDict[CCDkey]['File_name'].copy()
        
        File_name  = ['%s.fits' % item for item in rootFile_name]
        
        DataDict[CCDkey]['Files'] = np.array(File_name).copy()
        
    
    return DataDict, isconsistent
    
    
def prep_data_PSF0X(DataDict,RepDict,inputs,log=None):
    """Takes Raw Data and prepares it for further analysis. Also checks that
    data quality is enough."""
    
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


def basic_analysis_PSF0X(DataDict,RepDict,inputs,log=None):
    """Performs basic analysis on spots:
         - Gaussian Fits: peak, position, width_x, width_y
    """
    
def bayes_analysis_PSF0X(DataDict,RepDict,inputs,log=None):
    """ 
    Performs bayesian decomposition of the spot images:
        - optomechanic PSF and detector PSF.
    Also measures R2, fwhm and ellipticity of "extracted" detector PSF.
    """
    

def meta_analysis_PSF0X(DataDict,RepDict,inputs,log=None):
    """
    
    Analyzes the relation between detector PSF and fluence.
    
    """
    


def run(inputs,log=None):
    """Test PSF0X master function."""
    
    
    # INPUTS
    
    todo_flags = dict(init=True,prep=True,basic=True,bayes=True,meta=True,report=True)
    
    OBSID_lims = inputs['OBSID_lims']
    explogf = inputs['explogf']
    datapath = inputs['datapath']
    resultspath = inputs['resultspath']
    wavelength = inputs['wavelength']
    elvis = inputs['elvis']
    
    DataDictFile = os.path.join(resultspath,'PSF0X_%snm_DataDict.pick' % wavelength)
    reportobjFile = os.path.join(resultspath,'PSF0X_%snm_Report.pick' % wavelength)
    
    if not isthere(resultspath):
        os.system('mkdir %s' % resultspath)
    
    try: 
        structure = inputs['structure']
    except: 
        structure = PSF0X_structure
        
    try: reportroot = inputs['reportroot']
    except KeyError: reportroot = 'PSF0X_%inm_report' % wavelength
    
    try: cleanafter = inputs['cleanafter']
    except KeyError: cleanafter = False
    
    if 'todo_flags' in inputs: todo_flags.update(inputs['todo_flags'])
    
    doIndivs = []
    for key in ['basic','bayes']:
        if todo_flags[key] == True: doIndivs.append(True)
    
    resindivpath = '%s_indiv' % resultspath
    
    if not isthere(resindivpath) and np.any(doIndivs):
        os.system('mkdir %s' % resindivpath)
    
    
    if todo_flags['init']:
    
        # Initialising Report Object
    
        if todo_flags['report']:
            reportobj = Report(TestName='PSF0X: %s nm' % wavelength)
        else:
            reportobj = None
    
        # META-DATA WORK
        
    
        # Filter Exposures that belong to the test
    
        DataDict, isconsistent = filterexposures_PSF0X(wavelength,explogf,datapath,OBSID_lims,
                                     structure,elvis)
    
        if log is not None:
            log.info('PSF0X acquisition is consistent with expectations: %s' % isconsistent)
    
   
        # Add HK information
        DataDict = pilib.addHK(DataDict,HKKeys_PSF0X,elvis=elvis)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
        
    else:
        
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        
    # DATA-WORK
    
    # Prepare Data for further analysis (subtract offsets, divide by FFF, trim snapshots). 
    # Check Data has enough quality:
    #     median levels in pre-scan, image-area, overscan
    #     fluences and spot-sizes (coarse measure) match expectations for all spots
    #     
    
    if todo_flags['prep']:
    
        DataDict, reportobj = prep_data_PSF0X(DataDict,reportobj,inputs,log)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)        
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
    
    # Optional
    # Perform Basic Analysis : Gaussian fits and Moments
    
    if todo_flags['basic']:
        
        DataDict, reportobj = basic_analysis_PSF0X(DataDict,reportobj,inputs,log)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)        
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
    
    # Optional
    # Perform Bayesian Analysis
    
    
    if todo_flags['bayes']:
        DataDict, reportobj = bayes_analysis_PSF0X(DataDict,reportobj,inputs,log)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)  
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        
    
    # Optional
    # Produce Summary Figures and Tables
    
    if todo_flags['meta']:
        DataDict, reportobj = meta_analysis_PSF0X(DataDict,reportobj,inputs,log)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)  
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
    
    # Write automatic Report of Results
    
    if todo_flags['report']:
        reportobj.generate_Texbody()
        outfiles = reportobj.writeto(reportroot,cleanafter)
        
        for outfile in outfiles:
            os.system('mv %s %s/' % (outfile,resultspath))
    
    pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
    
    if log is not None:
        log.info('Finished PSF0X')
    
