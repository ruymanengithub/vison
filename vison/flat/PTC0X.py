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
import warnings

from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from copy import deepcopy
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

PTC0X_commvalues = dict(program='CALCAMP',
  flushes=7,exptime=0.,shuttr=1,
  siflsh=1,siflsh_p=500,
  wave = 4,
  source='flat',
  comments='')
  
PTC01_exptimes = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100.*ogse.tFWC_flat['nm800'] # ms
PTC02waves = [590,640,730,880]
PTC02TEMP_exptimes = exptimes=np.array([10.,30.,50.,70.,80.,90.])/100.*ogse.tFWC_flat['nm800']


testdefaults = dict(PTC01=dict(exptimes=PTC01_exptimes,
                         frames=[10,10,10,10,10,10,10,10,4,4,4]),
                    PTC02WAVE=dict(waves=PTC02waves,
                                   frames=[4,4,4,4,4,4],
                                   exptimes=dict()),
                    PTC02TEMP=dict(frames=[4,4,4,4,4,4]),
                                   exptimes=PTC02TEMP_exptimes)
                    

for w in testdefaults['PTC02WAVE']['waves']:
    testdefaults['PTC02WAVE']['exptimes']['nm%i' % w] = np.array([10.,30.,50.,70.,80.,90.])/100.*ogse.tFWC_flat['nm%i' % w]
    


def build_PTC0X_scriptdict(testkey,exptimes,frames,wavelength=800,diffvalues=dict(),
                           elvis='6.3.0'):
    """Builds PTC0X script structure dictionary.
    
    :param exptimes: list of ints [ms], exposure times.
    :param frames: list of ints, number of frames for each exposure time.
    :param wavelength: int, wavelength. Default: 800 nm.
    :param diffvalues: dict, opt, differential values.   
        
    """
    
    assert  len(exptimes) == len(frames)
    
    FW_ID = ogse.get_FW_ID(wavelength)
    FW_IDX = int(FW_ID[-1])
    

    PTC0X_commvalues['test'] = testkey
    PTC0X_commvalues['wave'] = FW_IDX

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




def filterexposures(structure,explogf,datapath,OBSID_lims,elvis='6.3.0'):
    """
    
    """    
    wavedkeys = []
    return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
                          wavedkeys=wavedkeys,elvis=elvis)
    


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
    
    if report is not None: report.add_Section(keyword='check_data',Title='Data Validation',level=0)
    
    
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
    if 'Sector' not in Sindices.names:
        Sindices.append(core.vIndex('Sector',vals=Sectors))
    

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
    
    warnings.warn('NOT FINISHED')
    
    


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
    
    raise NotImplementedError


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
    
    raise NotImplementedError


def feeder(inputs,elvis='6.3.0'):
    """ """
    
    subtasks = [('check',check_data),('extract',extract_PTC),
                ('meta',meta_analysis)]
    
    wavelength = inputs['wavelength']
    testkey = inputs['test']
    
    if testkey == 'PTC01': 
        _testkey = 'PTC01'
    if 'PTC02' in testkey:
        if testkey[-1] == 'K': _testkey = 'PTCTSEMP'
        else: _testkey = 'PTCWAVE'
    
    
    if 'exptimes' in inputs: 
        exptimes = inputs['exptimes']
    else:
        if _testkey == 'PTC02WAVE':
            exptimes = testdefaults[_testkey]['exptimes']['%nm%i' % wavelength]
        else:
            exptimes = testdefaults[_testkey]['exptimes']
    
    if 'frames' in inputs: frames = inputs['frames']
    else: frames = testdefaults[_testkey]['frames']
    
    
    if 'elvis' in inputs:
        elvis = inputs['elvis']
    if 'diffvalues' in inputs:
        diffvalues = inputs['diffvalues']
    else:
        diffvalues = {}
        
    scriptdict = build_PTC0X_scriptdict(testkey,exptimes,frames,wavelength=wavelength,
                           diffvalues=diffvalues,
                           elvis=elvis)
    
    inputs['structure'] = scriptdict
    inputs['subtasks'] = subtasks
    
    return inputs
    
