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
import copy
import string as st
from collections import OrderedDict

from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import HKtools
from vison.datamodel import core, ccd
from vison.image import calibration
import ptc as ptclib
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
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
               
    commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(PTC0X_commvalues)
    
    PTC0X_sdict = sc.update_structdict(PTC0X_sdict,commvalues,diffvalues)
    
    return PTC0X_sdict




def filterexposures(structure,explogf,datapath,OBSID_lims,elvis='6.3.0'):
    """
    
    """    
    wavedkeys = []
    return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
                          wavedkeys=wavedkeys,elvis=elvis)
    


def check_data(dd,report,inputs,log=None):
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
    

    newcolnames_off = ['offset_pre','offset_ove']
    for newcolname_off in newcolnames_off:
        dd.initColumn(newcolname_off,Qindices,dtype='float32',valini=np.nan)
    
    newcolnames_std = ['std_pre','std_img','std_ove']
    for newcolname_std in newcolnames_std:
        dd.initColumn(newcolname_std,Qindices,dtype='float32',valini=np.nan)
    

    dd.initColumn('chk_flu_img',Qindices,dtype='float32',valini=np.nan)
    
        
    
    nObs,nCCD,nQuad = Qindices.shape
    Quads = Qindices[2].vals
    CCDs = Qindices[1].vals
        
    # Get statistics in different regions
    
    if not bypass:
    
        for iObs in range(nObs):
            
            for jCCD in range(nCCD):
                
                #CCD = CCDs[jCCD]
                
                dpath = dd.mx['datapath'][iObs,jCCD]
                ffits = os.path.join(dpath,'%s.fits' % dd.mx['File_name'][iObs,jCCD])
                
                ccdobj = ccd.CCD(ffits)
                
                for kQ in range(nQuad):
                    Quad = Quads[kQ]
                    
                    for reg in ['pre', 'ove']:
                        altstats = ccdobj.get_stats(Quad,sector=reg,statkeys=['mean','std'],trimscan=[5,5],
                                ignore_pover=False,extension=-1)
                        dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = altstats[0]
                        dd.mx['std_%s' % reg][iObs,jCCD,kQ] = altstats[1]
                    
                    imgstats = ccdobj.get_stats(Quad,sector=reg,statkeys=['median','std'],trimscan=[5,5],
                                ignore_pover=True,extension=-1)
                    dd.mx['chk_flu_img'][iObs,jCCD,kQ] = imgstats[0]
                    dd.mx['std_img'][iObs,jCCD,kQ] = imgstats[1]
                    
                
    # Assess metrics are within allocated boundaries
    
    warnings.warn('NOT FINISHED')
    
    return dd, report


def extract_PTC(dd,report,inputs,log=None):
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
    
    if report is not None: report.add_Section(keyword='extract',Title='PTC Extraction',level=0)    
    
    bypass = False
    
    # HARDWIRED VALUES
    wpx = 300
    hpx = 300
    
    doMask = False
    
    if 'mask' in inputs['inCDPs']:
        Maskdata = calibration.load_CDPs(inputs['inCDPs']['Mask'],ccd.CCD)
        doMask = True
        if log is not None:
            masksstr = inputs['inCDPs']['Mask'].__str__()
            masksstr = st.replace(masksstr,',',',\n')
            log.info('Applying cosmetics mask')
            log.info(masksstr)    
    
    
    label = dd.mx['label'][:,0].copy() # labels should be the same accross CCDs. PATCH.
    ObsIDs = dd.mx['ObsID'][:].copy()
    
    indices = copy.deepcopy(dd.indices)
    
    nObs,nCCD,nQuad = indices.shape
    Quads = indices[indices.names.index('Quad')].vals
    CCDs = indices[indices.names.index('CCD')].vals
    
    emptyccdobj = ccd.CCD()    
    tile_coos = dict()
    for Quad in Quads:
        tile_coos[Quad] = emptyccdobj.get_tile_coos(Quad,wpx,hpx)
    Nsectors = tile_coos[Quads[0]]['Nsamps']    
    sectornames = np.arange(Nsectors)
    
    
    Sindices = copy.deepcopy(dd.indices)
    if 'Sector' not in Sindices.names:
        Sindices.append(core.vIndex('Sector',vals=sectornames))
    
    dd.initColumn('sec_med',Sindices,dtype='float32',valini=np.nan)
    dd.initColumn('sec_var',Sindices,dtype='float32',valini=np.nan)
    
    # Pairing ObsIDs
    
    dd.initColumn('ObsID_pair',dd.mx['ObsID'].indices,dtype='int64',valini=np.nan)
    
    ulabels = np.unique(label)
    
    for ulabel in ulabels:
        six = np.where(label == ulabel)
        nsix = len(six[0])
        ixeven = np.arange(0,nsix,2)
        ixodd = np.arange(1,nsix,2)
        
        dd.mx['ObsID_pair'][six[0][ixeven]] = ObsIDs[six[0][ixodd]]
        
        
    if not bypass:
        
        if doMask:
            estimators = dict(median=np.ma.median,std=np.ma.std)
        else:
            estimators = dict(median=np.median,std=np.std)
        
        for iObs in range(nObs):
            
            _ObsID_pair = dd.mx['ObsID_pair'][iObs]
            if np.isnan(_ObsID_pair): continue
            iObs_pair = np.where(ObsIDs == _ObsID_pair)[0][0]
            
            
            for jCCD,CCD in enumerate(CCDs):
                
                CCDkey = 'CCD%i' % CCD
                
                dpathodd = dd.mx['datapath'][iObs,jCCD]
                dpatheve = dd.mx['datapath'][iObs_pair,jCCD]
                
                ffits_odd = os.path.join(dpathodd,'%s.fits' % dd.mx['File_name'][iObs,jCCD])
                ffits_eve = os.path.join(dpatheve,'%s.fits' % dd.mx['File_name'][iObs_pair,jCCD])
                
                ccdobj_odd = ccd.CCD(ffits_odd)
                ccdobj_eve = ccd.CCD(ffits_eve)
                
                evedata = ccdobj_eve.extensions[-1].data.copy()
                
                ccdobj_odd.sub_bias(evedata,extension=-1)
                
                if doMask:
                    ccdobj_odd.get_mask(Maskdata[CCDkey].extensions[-1])
                
                for kQ in range(nQuad):
                    
                    Quad = Quads[kQ]
                    
                    _tile_coos = tile_coos[Quad]
                    
                    _tiles = ccdobj_odd.get_tiles(Quad,_tile_coos,extension=-1)
                    
                    _meds = np.array(map(estimators['median'],_tiles))
                    _vars = np.array(map(estimators['std'],_tiles))**2.
                    
                    dd.mx['sec_med'][iObs,jCCD,kQ,:] = _meds.copy()
                    dd.mx['sec_var'][iObs,jCCD,kQ,:] = _vars.copy()
                    
    
    return dd, report



def meta_analysis(dd,report,inputs,log=None):
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
                solve for blooming limit (ADU)
                    convert bloom limit to electrons, using gain
    
        plot PTC curves with best-fit f.e. CCD, Q
        report on gain estimates f. e. CCD, Q (table)
        report on blooming limits (table)
    
    """
    
    if report is not None: report.add_Section(keyword='meta',Title='PTC Analysis',level=0)
    
    dIndices = copy.deepcopy(dd.indices)
    
    CCDs = dIndices[dIndices.names.index('CCD')].vals
    Quads = dIndices[dIndices.names.index('Quad')].vals    
    
    # Initializations of output data-products
    
    gain_mx = OrderedDict()
    
    g_tmp_keys = ['a0','ea0','a1','ea1','a2','ea2','gain','egain','alpha','rn']
    
    for CCD in CCDs:
        CCDkey = 'CCD%i' % CCD
        gain_mx[CCDkey] = dict()
        for Quad in Quads:
            gain_mx[CCDkey][Quad] = dict()
            for key in g_tmp_keys:
                gain_mx[CCDkey][Quad][key] = np.nan
    
    b_tmp_keys = ['bloom_ADU','bloom_e']
    
    bloom_mx = OrderedDict()

    for CCD in CCDs:
        CCDkey = 'CCD%i' % CCD
        bloom_mx[CCDkey] = dict()
        for Quad in Quads:
            bloom_mx[CCDkey][Quad] = dict()
            for key in b_tmp_keys:
                bloom_mx[CCDkey][Quad][key] = np.nan
    
    # fitting the PTCs

    
    for iCCD, CCD in enumerate(CCDs):
        
        
        CCDkey = 'C%i' % CCD
        
        for jQ, Q in enumerate(Quads):
            
            
            raw_var = dd.mx['sec_var'][:,iCCD,jQ,:]
            raw_med = dd.mx['sec_med'][:,iCCD,jQ,:]
            
            ixnonan = np.where(~np.isnan(raw_var) & ~np.isnan(raw_med))
            var = raw_var[ixnonan]
            med = raw_med[ixnonan]
            
            #fitresults = dict(fit=p,efit=ep,gain=g,cuadterm=cuadterm,rn=rn,badresult=badresult)
            
            _fitresults = ptclib.fitPTC(med,var)
            
            for zx in range(2):
                gain_mx[CCDkey][Q]['a%i' % zx] = _fitresults['fit'][2-zx]
                gain_mx[CCDkey][Q]['ea%i' % zx] = _fitresults['efit'][2-zx]
                
            gain_mx[CCDkey][Q]['gain'] = _fitresults['gain']
            gain_mx[CCDkey][Q]['rn'] = _fitresults['rn']
            gain_mx[CCDkey][Q]['quality'] = _fitresults['quality']
            
            
            _bloom = ptclib.foo_bloom(med,var)
            
            bloom_mx[CCDkey][Q]['bloom_ADU'] = _bloom['bloom']
            bloom_mx[CCDkey][Q]['bloom_ADU'] = _bloom['bloom']
    

    dd.products['gain'] = copy.deepcopy(gain_mx)
    dd.products['bloom'] = copy.deepcopy(bloom_mx)
    
    # Build Tables
    
    # Do plots
    
    # Add rerts
    
    
    return dd, report

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
    
