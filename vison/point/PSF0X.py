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
from vison.datamodel import core
from vison.datamodel import ccd
from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.datamodel import HKtools
from vison.datamodel import scriptic as sc
from vison.flat import FlatFielding as FFing
from vison.point import lib as polib
from vison.support.report import Report
from vison.support import files
from vison.image import calibration
import copy
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

stampw = polib.stampw


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
    
    commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
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

    PSF0X: Checks quality of ingested data.
    
    **METACODE**
    
    
    ::
      
      check common HK values are within safe / nominal margins
      check voltages in HK match commanded voltages, within margins
    
      f.e.ObsID:
          f.e.CCD: 
              f.e.Q.:
                  measure offsets in pre-,and over-
                  measure std in pre-, over-
                  measure median in img area, excluding spots (bgd)
                  measure (bgd-sub'd-) fluence of spots
                  measure size of spots
                  
      assess std in pre- (~RON) is within allocated margins
      assess offsets in pre-, img-, over- are equal, within allocated  margins
      assess offsets are within allocated margins
      assess median-img is within allocated margins
      assess fluences of spots are within allocated margins
      assess sizes of spots are within allocated margins
    
      plot offsets vs. time
      plot std vs. time
      plot spot fluences vs. time
      plot spot size vs. time
          
      issue any warnings to log
      issue update to report
      update flags as needed
    
    """
    
    if log is not None:
        log.info('PSF0X.check_data')
    
    report.add_Section(keyword='check_data',Title='Data Validation',level=0)

    
    bypass = True # TESTS
    
    
    # CHECK AND CROSS-CHECK HK
    
    #print 'HK-perf' # TESTS
    report_HK_perf = HKtools.check_HK_vs_command(HKKeys,dd,limits='P',elvis=inputs['elvis'])
    #print 'HK-safe' # TESTS
    report_HK_safe = HKtools.check_HK_abs(HKKeys,dd,limits='S',elvis=inputs['elvis'])
    

    # Initialize new columns

    Qindices = copy.deepcopy(dd.indices)
    
    if 'Quad' not in Qindices.names:
        Qindices.append(core.vIndex('Quad',vals=pilib.Quads))

    
    Sindices = copy.deepcopy(Qindices)
    if 'Spot' not in Sindices.names:
        Sindices.append(core.vIndex('Spot',vals=polib.Point_CooNom['names']))
    

    newcolnames_off = ['offset_pre','offset_ove']
    for newcolname_off in newcolnames_off:
        dd.initColumn(newcolname_off,Qindices,dtype='float32',valini=np.nan)
    
    newcolnames_std = ['std_pre','std_ove']
    for newcolname_std in newcolnames_std:
        dd.initColumn(newcolname_std,Qindices,dtype='float32',valini=np.nan)
    

    dd.initColumn('bgd_img',Qindices,dtype='float32',valini=np.nan)
    
    SpotNames = Sindices[3].vals
    nSpot = len(SpotNames)
    
    dd.initColumn('chk_x',Sindices,dtype='float32',valini=np.nan)
    dd.initColumn('chk_y',Sindices,dtype='float32',valini=np.nan)
    dd.initColumn('chk_peak',Sindices,dtype='float32',valini=np.nan)
    dd.initColumn('chk_fluence',Sindices,dtype='float32',valini=np.nan)
    dd.initColumn('chk_fwhmx',Sindices,dtype='float32',valini=np.nan)
    dd.initColumn('chk_fwhmy',Sindices,dtype='float32',valini=np.nan)
    
    chkkeycorr = dict(chk_x='x',chk_y='y',chk_peak='peak',chk_fwhmx='fwhmx',chk_fwhmy='fwhmy')
    
    nObs,nCCD,nQuad = Qindices.shape
    Quads = Qindices[2].vals
    CCDs = Qindices[1].vals
    
    # Get statistics in different regions
    
    if not bypass:
    
        for iObs in range(nObs):
            
            for jCCD in range(nCCD):
                
                CCD = CCDs[jCCD]
                
                dpath = dd.mx['datapath'][iObs,jCCD]
                ffits = os.path.join(dpath,'%s.fits' % dd.mx['File_name'][iObs,jCCD])
                
                ccdobj = ccd.CCD(ffits)
                
                for kQ in range(nQuad):
                    Quad = Quads[kQ]
                    
                    for reg in ['pre', 'ove']:
                        altstats = ccdobj.get_stats(Quad,sector=reg,statkeys=['mean','std'],trimscan=[5,5],
                                ignore_pover=True,extension=-1)
                        dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = altstats[0]
                        dd.mx['std_%s' % reg][iObs,jCCD,kQ] = altstats[1]
                    
                    
                    # To measure the background we mask out the sources
                    
                    alt_ccdobj = copy.deepcopy(ccdobj)
                    
                    mask_sources = polib.gen_point_mask(CCD,Quad,width=stampw,sources='all')
                    
                    alt_ccdobj.get_mask(mask_sources)
                    
                    imgstats = alt_ccdobj.get_stats(Quad,sector='img',statkeys=['median'],trimscan=[5,5],
                                ignore_pover=True,extension=-1)
                    
                    dd.mx['bgd_img'][iObs,jCCD,kQ] = imgstats[0]
                    
                    alt_ccdobj = None
                    
                    for xSpot in range(nSpot):
                        SpotName = SpotNames[xSpot]
                        
                        coo = polib.Point_CooNom['CCD%i' % CCD][Quad][SpotName]

                        spot = polib.extract_spot(ccdobj,Quad,coo,log=log,
                                                  stampw=stampw)
                        
                        res_bas = spot.measure_basic(rin=10,rap=10,rout=-1)
                        
                        for chkkey in chkkeycorr:
                            dd.mx[chkkey][iObs,jCCD,kQ,xSpot] = res_bas[chkkeycorr[chkkey]]
                        

                    stop()
                    
                
    # Assess metrics are within allocated boundaries
    
    

#    offset_lims = [1000.,3500.] # TESTS, should be in common limits file
#    offset_diffs = dict(img=[-1.,+5.],ove=[-1.,+6.]) # TESTS, should be in common limits file
#    
#    std_lims = [0.5,2.] # TESTS, should be in common limits file
#    
#    # absolute value of offsets
#    
#    compliance_offsets = dict()
#    for reg in ['pre','img','ove']:
#        
#        test = ((dd.mx['offset_%s' % reg][:] <= offset_lims[0]) |\
#                      (dd.mx['offset_%s' % reg][:] >= offset_lims[1]))
#        compliance_offsets[reg] = np.any(test,axis=(1,2)).sum()
#    
#    
#    # cross-check of offsets
#        
#    xcheck_offsets = dict()
#    for reg in ['img','ove']:
#        
#        test = dd.mx['offset_%s' % reg][:]-dd.mx['offset_pre'][:]
#        
#        testBool = (test <= offset_diffs[reg][0]) | \
#               (test >= offset_diffs[reg][1])
#        xcheck_offsets['img'] = np.any(testBool,axis=(1,2)).sum()
#    
#    # absolute value of std
#        
#    compliance_std = dict()
#    for reg in ['pre']:
#        
#        test = ((dd.mx['std_%s' % reg][:] <= std_lims[0]) |\
#                      (dd.mx['std_%s' % reg][:] >= std_lims[1]))
#        compliance_std[reg] = np.any(test,axis=(1,2)).sum()
#    
#    # Do some Plots
#    
#    # offsets vs. time
#    # std vs. time
#    
#    if log is not None:
#        log.info('Plotting MISSING in BIAS01.check_data')
#    
#    # Update Report, raise flags, fill-in log
#
    if log is not None:
        log.info('Reporting and Flagging MISSING in PSF0X.check_data')    
    
    
    return dd, report
    
    


def prep_data(dd,report,inputs,log=None):
    """
    
    PSF0X: Preparation of data for further analysis.
    
    **METACODE**
    
    ::
        f.e. ObsID:
            f.e.CCD:                
                apply cosmetic mask, if available
                f.e.Q:
                    subtract offset
                subtract superbias, if available
                divide by flat-field, if available
                
                save image as a datamodel.ccd.CCD object.                
                cuts-out and save stamps of pre-processed spots for further analysis.
    
    """
    
    report.add_Section(keyword='prep_data',Title='Data Pre-Processing',level=0)
    
    bypass = True # TESTS
    
    # Inputs un-packing
    
    doMask = False
    doBias = False
    doOffset = True
    doFlats = False
    
    if 'mask' in inputs['inCDPs']:
        Maskdata = calibration.load_CDPs(inputs['inCDPs']['Mask'],ccd.CCD)
        doMask = True
        
        
    if 'bias' in inputs['inCDPs']:
        Biasdata = calibration.load_CDPs(inputs['inCDPs']['Bias'],ccd.CCD)
        doBias = True
    
    
    if 'FF' in inputs['inCDPs']:
        FFdata = calibration.load_CDPs(inputs['inCDPs']['FF'],FFing.FlatField)
        doFlats = True
    
    # INDEX HANDLING
    
    dIndices = copy.deepcopy(dd.indices)
    
    CCDs = dIndices[dIndices.names.index('CCD')].vals
    #nCCD = len(CCDs)
    Quads = dIndices[dIndices.names.index('Quad')].vals
    nQuads = len(Quads)
    nObs = dIndices[dIndices.names.index('ix')].len
    SpotNames = dIndices[dIndices.names.index('Spot')].vals
    nSpots = len(SpotNames)
    
    CIndices = copy.deepcopy(dIndices)
    CIndices.pop(CIndices.names.index('Quad'))
    CIndices.pop(CIndices.names.index('Spot'))
    
    # INITIALIZATIONS
    
    dd.initColumn('ccdobj_name',CIndices,dtype='S100',valini='None')
    dd.initColumn('spots_name',CIndices,dtype='S100',valini='None')
    
    if not bypass:
        
        rpath = dd.meta['inputs']['resultspath']      
        
        for iObs in range(nObs):
            
            for jCCD,CCD in enumerate(CCDs):
                
                CCDkey = 'CCD%i' % CCD
                
                dd.mx['ccdobj_name'][iObs,jCCD] = '%s_proc' % dd.mx['File_name'][iObs,jCCD]
                
                dpath = dd.mx['datapath'][iObs,jCCD]
                infits = os.path.join(dpath,'%s.fits' % dd.mx['File_name'][iObs,jCCD])
                
                ccdobj = ccd.CCD(infits)
                
                fullccdobj_name = os.path.join(rpath,'%s.pick' % dd.mx['ccdobj_name'][iObs,jCCD]) 
                
                if doMask:
                    ccdobj.get_mask(Maskdata[CCDkey].extensions[-1])
                
                if doOffset:
                    for Quad in Quads:
                        ccdobj.sub_offset(Quad,method='median',scan='pre',trims=[5,5],
                                          ignore_pover=False)
                if doBias:
                    ccdobj.sub_bias(Biasdata[CCDkey].extensions[-1],extension=-1)
                
                if doFlats:
                    ccdobj.divide_by_flatfield(FFdata[CCDkey].extensions[1],extension=-1)
                
                ccdobj.writeto(fullccdobj_name,clobber=True)
                
                
                # Cut-out "spots"
                
                spots_array = np.zeros((nQuads,nSpots),dtype=object)
                
                for kQ,Quad in enumerate(Quads):
                    
                    for lS,SpotName in enumerate(SpotNames):
                    
                        coo = polib.Point_CooNom[CCDkey][Quad][SpotName]
                        lSpot = polib.extract_spot(ccdobj,Quad,coo,stampw=stampw)
                        
                        spots_array[kQ,lS] = copy.deepcopy(lSpot)
                
                
                # save "spots" to a separate file and keep name in dd
                
                dd.mx['spots_name'][iObs,jCCD] = '%s_spots' % dd.mx['File_name'][iObs,jCCD]
                
                fullspots_name = os.path.join(rpath,'%s.pick' % dd.mx['spots_name'][iObs,jCCD])
                
                spotsdict = dict(spots=spots_array)
                
                files.cPickleDumpDictionary(spotsdict,fullspots_name)
                
                
                
    
    return dd,report




def basic_analysis(dd,report,inputs,log=None):
    """Performs basic analysis on spots:
         - shape from moments
         - Gaussian fitting: peak intensity, position, width_x, width_y
         
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
