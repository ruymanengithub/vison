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
from vison.pipe import lib as pilib
from vison.pipe import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
import datetime
# END IMPORT

isthere = os.path.exists

HKKeys_FOCUS00 = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

dtobj_default = datetime.datetime(1980,2,21,7,0,0) # early riser


def get_FOCUS00_structure(wavelength):
    """ """
    
    FilterPos = [key for key in pilib.FW if pilib.FW[key] == '%inm' % wavelength][0]    
    mirror_nom = polib.mirror_nom[FilterPos]

    FOCUS00_structure = dict(col1=dict(N=5,Exptime=0,Mirr_pos=mirror_nom-0.2),
                          col2=dict(N=2,Exptime=10.,Mirr_pos=mirror_nom-0.2),
                          col3=dict(N=2,Exptime=10.,Mirr_pos=mirror_nom-0.1),
                          col4=dict(N=2,Exptime=10.,Mirr_pos=mirror_nom-0.0),
                          col5=dict(N=2,Exptime=10.,Mirr_pos=mirror_nom+0.1),
                          col6=dict(N=2,Exptime=10.,Mirr_pos=mirror_nom+0.2),
                   Ncols=6)

    return FOCUS00_structure

FOCUS00_structure_wnom = get_FOCUS00_structure(800)

def filterexposures_FOCUS00(inwavelength,explogf,datapath,OBSID_lims,structure=FOCUS00_structure_wnom,elvis='5.7.04'):
    """Loads a list of Exposure Logs and selects exposures from test FOCUS00.
    
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
    
    selbool = (['FOCUS' in item for item in explog['TEST']]) & \
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
        
        
        Mirr_pos = DataDict[CCDkey]['Mirr_pos'].copy()
        Exptime = DataDict[CCDkey]['Exptime'].copy()
        
        label = np.zeros(Nobs,dtype='40str')
        
        uMirr_pos = np.sort(np.unique(Mirr_pos))
        
        
        for ixMP,MP in enumerate(uMirr_pos):
            
            ixsel = np.where(Mirr_pos == uMirr_pos)
            
            iexptime = Exptime[ixsel]
            if iexptime >0:
                label[ixsel] = 'focus_%i' % ixMP
            else:
                label[ixsel] = 'BGD'
        
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
    
    
    cols_to_add = ['offset_%s','std_pre_%s','std_ove_%s','med_img_%s',
                   'std_img_%s']
                   
    # Loop over CCDs
    
    for CCDindex in range(1,4):
        
        CCDkey = 'CCD%i' % CCDindex
        
        if CCDkey in DataDict:
            
            # Loop over ObsIDs
            
            ObsIDs = DataDict[CCDkey]['ObsID'].copy()
            label = DataDict[CCDkey]['label'].copy()
            #Nobs = len(ObsIDs)
            Nobs = len(ObsIDs)
            
            for col in cols_to_add:
                for Q in pilib.Quads:
                    DataDict[CCDkey][col % Q] = np.zeros(Nobs,dtype='float32')
            
            
            for iObs, ObsID in enumerate(ObsIDs):
                
                ilabel = label[iObs]
                
                # log object-id being analyzed: ObsID
                
                if log is not None: log.info('working on ObsID %i' % ObsID)
                
                # retrieve FITS file and open it
                                
                idatapath = DataDict[CCDkey]['datapath'][iObs]
                
                fitsf = os.path.join(idatapath,DataDict[CCDkey]['Files'][iObs])
                
                
                # subtract offset and add to DataDict
                # measure STD in pre and over-scan and add to DataDict
                                
                ccdobj = ccd.CCD(fitsf)
                
                
                for Q in pilib.Quads:
                    
                    Qoffset = ccdobj.sub_offset(Q,method='median',scan='pre',trimscan=[5,5])[0]
                    DataDict[CCDkey]['offset_%s' % Q][iObs] = Qoffset
                    
                    
                    ignore,std_pre = ccdobj.get_stats(Q,area='pre',detail=False,trimscan=[5,5])
                    ignore,std_ove = ccdobj.get_stats(Q,area='ove',detail=False,trimscan=[5,5])
                    med_img,std_img = ccdobj.get_stats(Q,area='img',detail=False,trimscan=[10,10])
                    
                    DataDict[CCDkey]['std_pre_%s' % Q][iObs] = std_pre
                    DataDict[CCDkey]['std_ove_%s' % Q][iObs] = std_ove
                    DataDict[CCDkey]['std_img_%s' % Q][iObs] = std_img
                    DataDict[CCDkey]['med_img_%s' % Q][iObs] = med_img
                
                # Divide by flat-field
                
                
                iFF = FFdata[CCDkey].copy()
                
                ccdobj.divide_by_flatfield(iFF)
                
                
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


def generate_Explog_FOCUS00(wavelength,struct,elvis='6.0.0',date=dtobj_default):
    """ """
    
    Nscriptcols = struct['Ncols']
    
    columns = ELtools.columnlist[elvis]
    
    defaults = {'ObsID':0,'File_name':'','CCD':0,
    'ROE':'R01','DATE':'',
    'PROGRAM':'CALFM','TEST':'FOCUS00',
    'CCD1_SN':'C01','CCD2_SN':'C02','CCD3_SN':'C03',
    'BUNIT':'ADU','Operator':'x',
    'Lab_ver':'x.x.x','Con_file':'xxx.con',
    'Exptime':0,
    'Flsh-Rdout_e_time':0.,'C.Inj-Rdout_e_time':0.,'N_P_high':'I1I2I3',
    'Chrg_inj':0,'On_cycle':0,'Off_cycl':0,
    'Rpeat_cy':0,'pls_leng':0,'pls_del':0,
    'SerRdDel':0,'Trappump':0,'TP_Ser_S':0,
    'TP_Ver_S':0,'TP_DW_V':0,'TP_DW_H':0,'TOI_flsh':143,'TOI_pump':1000,
    'TOI_read':1000,'TOI_CInj':1000,'Invflshp':500,'Invflush':1,
    'Flushes':3,'Vstart':1,'Vend':2066,'Ovrscn_H':0,'CLK_ROE':'Normal',
    'CnvStart':0,'SumWell':0,'IniSweep':1,'SPW_clk':1,'FPGA_ver':'x.x.x',
    'EGSE_ver':elvis,'M_Steps':0,'M_St_Sze':0,
    'Wavelength':wavelength,'Mirr_pos':0,
    'RPSU_SN':'RP01','ROE_SN':'RO01','CalScrpt':'FakeScriptFOCUS00',
    'R1CCD1TT':153,'R1CCD1TB':153,'R1CCD2TT':153,'R1CCD2TB':153,'R1CCD3TT':153,
    'R1CCD3TB':153,'IDL_V':13000,'IDH_V':18000,'IG1_T1_V':4000,
    'IG1_T2_V':6000,'IG1_T3_V':4000,'IG1_B1_V':4000,'IG1_B2_V':6000,'IG1_B3_V':4000,
    'IG2_T_V':6000,'IG2_B_V':6000,'OD_T1_V':26000,'OD_T2_V':26000,'OD_T3_V':26000,
    'OD_B1_V':26000,'OD_B2_V':26000,'OD_B3_V':26000,'RD_T_V':17000,'RD_B_V':17000}
    
    explog = ELtools.iniExplog(elvis)
    
    ixObsID = 1000
    
    for iscrcol in range(1,Nscriptcols+1):
        scriptcol = struct['col%i' % iscrcol]
        N = scriptcol['N']
        inputkeys = [key for key in scriptcol.keys() if key != 'N']
        
        rowdict = {}
        
        for subixrow in range(N):
        
            for ecol in columns:
                rowdict[ecol] = defaults[ecol]
        
            for key in inputkeys:
                rowdict[key] = scriptcol[key]
            
            for ixCCD in range(1,4):
                
                dmy = date.strftime('%d%m%y')
                hms = date.strftime('%H%M%S')
            
                rowdict['ObsID'] = ixObsID
                rowdict['File_name'] = 'EUC_%i_%sD_%sT_ROE1_CCD%i' % (ixObsID,dmy,hms,ixCCD)
                rowdict['DATE'] = '%sD%sT' % (dmy,hms)
                rowdict['CCD'] = 'CCD%i' % ixCCD
                
                explog.add_row(vals=[rowdict[key] for key in columns])
                
                date = date + datetime.timedelta(seconds=90)
            
            ixObsID += 1
                    
            
    return explog

def get_dtobj(DT):
    
        date = DT[0:DT.index('D')]
        y2d = int(date[4:6])
        if y2d < 20: century = 2000
        else: century = 1900
        dd,MM,yy = int(date[0:2]),int(date[2:4]),y2d+century 
        
        time = DT[DT.index('D')+1:-1]
        
        hh,mm,ss = int(time[0:2]),int(time[2:4]),int(time[4:6])
    
        dtobj = datetime.datetime(yy,MM,dd,hh,mm,ss)
        return dtobj

def generate_HK_FOCUS00(explog,datapath,elvis='6.0.0'):
    """ """
    
    HKkeys = HKtools.allHK_keys[elvis]
    
    Nobs = len(explog['ObsID'])
    
    defaults = {'TimeStamp':'','HK_OD_Top_CCD1':27.,'HK_OD_Bottom_CCD1':27.,
'HK_OD_Top_CCD2':27.,'HK_OD_Bottom_CCD2':27.,'HK_OD_Top_CCD3':27.,'HK_OD_Bottom_CCD3':27.,
'HK_IG1_Top_CCD1':0.,'HK_IG1_Bottom_CCD1':0.,'HK_IG1_Top_CCD2':0.,'HK_IG1_Bottom_CCD2':0.,
'HK_IG1_Top_CCD3':0.,'HK_IG1_Bottom_CCD3':0.,'HK_temp_top_CCD1':153.,'HK_temp_bottom_CCD1':153.,
'HK_temp_top_CCD2':153.,'HK_temp_bottom_CCD2':153.,'HK_temp_top_CCD3':153.,
'HK_temp_bottom_CCD3':153.,'HK_RD_top':17.,'HK_RD_bot':17.,'HK_IG2_top':0.,'HK_IG2_bot':0.,
'HK_IDH':18.,'HK_IDL':13.,'HK_DD_bias':20.,'HK_OG_bias':2.,'HK_1.5V_ROE':0.,'HK_VCCD_ROE':0.,
'HK_5VA_pos_ROE':0.,'HK_5V_ref_ROE':0.,'HK_10VA_ROE':0.,'HK_5.2V_neg_ROE':0.,'HK_3V_neg_ROE':0.,
'HK_VRclk_ROE':0.,'HK_VRClk_Lo_ROE':0.,'HK_3.3V_DIG_RPSU':0.,'HK_I3.3V_DIG_RPSU':0.,
'HK_1.5V_DIG_RPSU':0.,'HK_I1.5V_DIG_RPSU':0.,'HK_28V_Pri_RPSU':0.,'HK_I28V_RPSU':0.,
'HK_VAN_pos_RPSU':0.,'HK_I+VAN_RPSU':0.,'HK_VAN_neg_RPSU':0.,'HK_I-VAN_RPSU':0.,
'HK_VCLK_RPSU':0.,'HK_IVCLK_RPSU':0.,'HK_VCCD_RPSU':0.,'HK_IVCCD_RPSU':0.,'HK_Temp1_RPSU':30.,
'HK_Temp2_RPSU':30.,'HK_Video_TOP':30.,'HK_Video_BOT':30.,'HK_FPGA_TOP':30.,'HK_FPGA_BOT':30.,
'HK_ID1':0.,'HK_ID2':0.,'HK_Viclk_ROE':0.}
    
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        idate = explog['DATE'][ixobs]
        
        idtobj = get_dtobj(idate)
        
        if ixobs < Nobs-1:
            ip1dtobj = get_dtobj(explog['DATE'][ixobs+1])
            dt = (ip1dtobj-idtobj).seconds
        else:
            dt = 90
        
        HKfilef = 'HK_%s_%s_ROE1.txt' % (obsid,idate)
        
        HKfilef = os.path.join(datapath,HKfilef)
        
        HKfile = HKtools.iniHK_QFM(elvis)
        
        for sec in range(dt):
            iidtobj = idtobj + datetime.timedelta(seconds=sec)
            
            iTimeStamp = iidtobj.strftime('%H:%M:%S')
            
            rowdict = {}
            
            for HKkey in HKkeys:
                rowdict[HKkey] = defaults[HKkey]
                
            rowdict['TimeStamp'] = iTimeStamp
            
            HKfile.add_row(vals=[rowdict[key] for key in HKkeys])
        
        
        HKfile.write(HKfilef,format='ascii',overwrite=True)
        

def generate_FITS_FOCUS00(explog,datapath,elvis='6.0.0'):
    """ """
    
    NAXIS1,NAXIS2 = 4238,4132
    
    maxexptime = explog['Exptime'].max()
    flatlevel = 100.
    biaslevel = 2000.
    
    waivedkeys = ['File_name','Flsh-Rdout_e_time','C.Inj-Rdout_e_time',
                  'Wavelength']
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        idate = explog['DATE'][ixobs]
        iCCD = explog['CCD'][ixobs]
        iexptime = explog['Exptime'][ixobs]
        
        idtobj = get_dtobj(idate)
        
        dmy = idtobj.strftime('%d%m%y')
        HMS = idtobj.strftime('%H%M%S')

        FITSf = 'EUC_%s_%sD_%sT_ROE1_%s.fits' % \
            (obsid,dmy,HMS,iCCD)
        
        FITSf = os.path.join(datapath,FITSf)
        
        ccdobj = ccd.CCD()
        
        img = np.zeros(shape=(NAXIS1,NAXIS2),dtype='float32')
        
        ccdobj.add_extension(data=None)
        ccdobj.add_extension(data=img,label='ROE1_%s' % iCCD)
        
        ccdobj.simadd_flatilum(levels=dict(E=flatlevel*iexptime*1.,
                                           F=flatlevel*iexptime*1.1,
                                           G=flatlevel*iexptime*1.2,
                                           H=flatlevel*iexptime*1.3))
        ccdobj.simadd_poisson()
        
        ccdobj.simadd_bias(levels=dict(E=biaslevel*1.,
                                           F=biaslevel*1.1,
                                           G=biaslevel*1.2,
                                           H=biaslevel*1.3))
        ccdobj.simadd_ron()
        
        ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')
        

        ccdobj.extensions[-1].header['WAVELENG'] = explog['Wavelength'][ixobs]
        
        for key in ELtools.columnlist[elvis]:
            if key not in waived:
                ccdobj.extensions[-1].header[key] = explog[key][ixobs]
        
        
        ccdobj.writeto(FITSf,clobber=True,unsigned16bit=True)
        
        stop()
    


def generate_Fake_FOCUS00(wavelength,date=dtobj_default,rootpath=''):
    """Generates Fake FOCUS00 data"""
    
    # Generate Exposure Log
    # Generate HK files
    # Generate FITS files
    
    FOCUS00_structure = get_FOCUS00_structure(wavelength)
    
    dmmmy = date.strftime('%d_%b_%y')
    dmy = date.strftime('%d%m%y')
    
    datapath = os.path.join(rootpath,dmmmy)
    
    if not isthere(datapath):
        os.system('mkdir %s' % datapath)
    
    explog = generate_Explog_FOCUS00(wavelength,FOCUS00_structure,elvis='6.0.0',
                                     date=date)
                                     
    explogf = os.path.join(datapath,'EXP_LOG_%s.txt' % dmy)
    
    explog.write(explogf,format='ascii',overwrite=True)
    
    #generate_HK_FOCUS00(explog,datapath,elvis='6.0.0')
    
    
    generate_FITS_FOCUS00(explog,datapath,elvis='6.0.0')
    
    
    
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
        structure = get_FOCUS00_structure(wavelength)
        
    try: reportroot = inputs['reportroot']
    except KeyError: reportroot = 'FOCUS00_%inm_report' % wavelength
    
    try: cleanafter = inputs['cleanafter']
    except KeyError: cleanafter = False
    
    if 'todo_flags' in inputs: todo_flags.update(inputs['todo_flags'])
        
    
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
        DataDict = pilib.addHK(DataDict,HKKeys_FOCUS00,elvis=elvis)
        save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
        
    else:
        
        DataDict, reportobj = recover_progress(DataDictFile,reportobjFile)
        
    # DATA-WORK
    
    # Prepare Data for further analysis (subtract offsets, divide by FFF, trim snapshots). 
    # Check Data has enough quality:
    #     median levels in pre-scan, image-area, overscan
    #     fluences and spot-sizes (coarse measure) match expectations for all spots
    
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


if __name__ == '__main__':
    
    pass