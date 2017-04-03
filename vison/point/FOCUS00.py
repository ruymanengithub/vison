# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FOCUS00

Focus analysis script

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check quality of data (integrated fluxes are roughly constant, matching expected level).
    - Subtract offset level.
    - Divide by Flat-field.
    - Crop stamps of the sources on each CCD/Quadrant.
       - save snapshot figures of sources.
    - for each source (5 x Nquadrants):
       - measure shape using Gaussian Fit
    - Find position of mirror that minimizes PSF sizes
    - Produce synoptic figures:
        source size and ellipticity across combined FOV (of 3 CCDs)
    - Save results.

Created on Mon Apr 03 16:21:00 2017

:author: Ruyman Azzollini
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

HKKeys_FOCUS00 = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

FOCUS00_structure = dict(col1=dict(N=5,Exptime=0,mirror=50),
                          col2=dict(N=2,Exptime=10.,mirror=45),
                          col3=dict(N=2,Exptime=10.,mirror=46),
                          col4=dict(N=2,Exptime=10.,mirror=47),
                          col5=dict(N=2,Exptime=10.,mirror=48),
                          col6=dict(N=2,Exptime=10.,mirror=49),
                          col7=dict(N=2,Exptime=10.,mirror=50),
                          col8=dict(N=2,Exptime=10.,mirror=51),
                          col9=dict(N=2,Exptime=10.,mirror=52),
                          col10=dict(N=2,Exptime=10.,mirror=53),
                          col11=dict(N=2,Exptime=10.,mirror=54),
                          col12=dict(N=2,Exptime=10.,mirror=55),
                   Ncols=12)


def filterexposures_FOCUS00(inwavelength,explogf,datapath,OBSID_lims,structure=FOCUS00_structure,elvis='5.7.04'):
    """Loads a list of Exposure Logs and selects exposures from test FOCUS00.
    
    The filtering takes into account an expected structure for the 
    acquisition script.

    The datapath becomes another column in DataDict. This helps dealing
    with tests that run overnight and for which the input data is in several
    date-folders.

    
    """
    
    # load exposure log(s)
    explog = plib.loadexplogs(explogf,elvis=elvis,addpedigree=True)
    
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
        
    
    for key in plib.FW.keys():
        if plib.FW[key] == '%inm' % inwavelength: Filter = key
    
    Filter = '%i' % inwavelength # TESTS
    
    selbool = (['PSF02' in item for item in explog['TEST']]) & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) & \
        (explog['Wavelength'] == Filter) # TESTS
    
    
    # Assess structure
    
    isconsistent = plib.check_test_structure(explog,selbool,structure)
    
    
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
    
    
def prep_data_FOCUS00(DataDict,RepDict,inputs,log=None):
    """Takes Raw Data and prepares it for further analysis. 
    Also checks that data quality is enough."""
    
    # Inputs un-packing
    
    FFs = inputs['FFs'] # flat-fields for each CCD (1,2,3)
    
    
    # Load Flat-Field data for each CCD
    
    FFdata = dict()
    for CCDindex in range(1,4):
        CCDkey = 'CCD%i' % CCDindex
        if CCDkey in FFs.keys():
            FFdata[CCDkey] = FFing.FlatField(fitsfile=FFs[CCDkey])
        
    RepDict['prepFOCUS00'] = dict(Title='Data Pre-Processing and QA',items=[])
    
    
    # Loop over CCDs
    
    for CCDindex in range(1,4):
        
        CCDkey = 'CCD%i' % CCDindex
        
        if CCDkey in DataDict:
            
            # Loop over ObsIDs
            
            ObsIDs = DataDict[CCDkey]['ObsID'].copy()
            label = DataDict[CCDkey]['label'].copy()
            #Nobs = len(ObsIDs)
            
            
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


def basic_analysis_FOCUS00(DataDict,RepDict,inputs,log=None):
    """Performs basic analysis on spots:
         - Gaussian Fits: peak, position, width_x, width_y
    """
    


def meta_analysis_FOCUS00(DataDict,RepDict,inputs,log=None):
    """
    
    Analyzes the relation between PSF shape and mirror position.
    
    """
    

def save_progress(DataDict,reportobj,DataDictFile,reportobjFile):        
    files.cPickleDumpDictionary(DataDict,DataDictFile)
    files.cPickleDump(reportobj,reportobjFile)


def recover_progress(DataDictFile,reportobjFile):
    DataDict = files.cPickleRead(DataDictFile)
    reportobj = files.cPickleRead(reportobjFile)
    return DataDict,reportobj

def run(inputs,log=None):
    """Test FOCUS00 master function."""
    
    
    # INPUTS
    
    todo_flags = dict(init=True,prep=True,basic=True,meta=True,report=True)
    
    OBSID_lims = inputs['OBSID_lims']
    explogf = inputs['explogf']
    datapath = inputs['datapath']
    resultspath = inputs['resultspath']
    wavelength = inputs['wavelength']
    elvis = inputs['elvis']
    
    DataDictFile = os.path.join(resultspath,'FOCUS00_%snm_DataDict.pick' % wavelength)
    reportobjFile = os.path.join(resultspath,'FOCUS00_%snm_Report.pick' % wavelength)
    
    if not isthere(resultspath):
        os.system('mkdir %s' % resultspath)
    
    try: 
        structure = inputs['structure']
    except: 
        structure = FOCUS00_structure
        
    try: reportroot = inputs['reportroot']
    except KeyError: reportroot = 'FOCUS00_%inm_report' % wavelength
    
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
            reportobj = Report(TestName='FOCUS00: %s nm' % wavelength)
        else:
            reportobj = None
    
        # META-DATA WORK
        
    
        # Filter Exposures that belong to the test
    
        DataDict, isconsistent = filterexposures_FOCUS00(wavelength,explogf,datapath,OBSID_lims,
                                     structure,elvis)
    
        if log is not None:
            log.info('FOCUS00 acquisition is consistent with expectations: %s' % isconsistent)
    
   
        # Add HK information
        DataDict = plib.addHK(DataDict,HKKeys_FOCUS00,elvis=elvis)
        save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
        
    else:
        
        DataDict, reportobj = recover_progress(DataDictFile,reportobjFile)
        
    # DATA-WORK
    
    # Prepare Data for further analysis (subtract offsets, divide by FFF, trim snapshots). 
    # Check Data has enough quality:
    #     median levels in pre-scan, image-area, overscan
    #     fluences and spot-sizes (coarse measure) match expectations for all spots
    #     
    
    if todo_flags['prep']:
    
        DataDict, reportobj = prep_data_FOCUS00(DataDict,reportobj,inputs,log)
        save_progress(DataDict,reportobj,DataDictFile,reportobjFile)        
    else:
        DataDict, reportobj = recover_progress(DataDictFile,reportobjFile)
    
    # Optional
    # Perform Basic Analysis : Gaussian fits and Moments
    
    if todo_flags['basic']:
        
        DataDict, reportobj = basic_analysis_FOCUS00(DataDict,reportobj,inputs,log)
        save_progress(DataDict,reportobj,DataDictFile,reportobjFile)        
    else:
        DataDict, reportobj = recover_progress(DataDictFile,reportobjFile)
    
    
    # Optional
    # Produce Summary Figures and Tables
    
    if todo_flags['meta']:
        DataDict, reportobj = meta_analysis_FOCUS00(DataDict,reportobj,inputs,log)
        save_progress(DataDict,reportobj,DataDictFile,reportobjFile)  
    else:
        DataDict, reportobj = recover_progress(DataDictFile,reportobjFile)
    
    # Write automatic Report of Results
    
    if todo_flags['report']:
        reportobj.generate_Texbody()
        outfiles = reportobj.writeto(reportroot,cleanafter)
        
        for outfile in outfiles:
            os.system('mv %s %s/' % (outfile,resultspath))
    
    save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
    
    if log is not None:
        log.info('Finished FOCUS00')
    