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
from vison.point.spot import Spot
from vison.point import display as pdspl
import datetime
# END IMPORT

isthere = os.path.exists

HKKeys_FOCUS00 = ['HK_temp_top_CCD1','HK_temp_bottom_CCD1','HK_temp_top_CCD2',
'HK_temp_bottom_CCD2','HK_temp_top_CCD3','HK_temp_bottom_CCD3'] # TESTS

dtobj_default = datetime.datetime(1980,2,21,7,0,0) # early riser

stampw = 25

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
        
    
    #for key in pilib.FW.keys():
    #    if pilib.FW[key] == '%inm' % inwavelength: Filter = key
    
    #Filter = '%i' % inwavelength # TESTS
    
    selbool = (['FOCUS00' in item for item in explog['TEST']]) & \
        (explog['ObsID'] >= OBSID_lims[0]) & \
        (explog['ObsID'] <= OBSID_lims[1]) & \
        (explog['Wavelength'] == inwavelength) # TESTS

    
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
            DataDict[CCDkey][key] = explog[key][CCDselbool].data.copy()
        
        DataDict[CCDkey]['time'] = np.array(map(pilib.get_dtobj,DataDict[CCDkey]['DATE'])).copy()
        
        
        Mirr_pos = DataDict[CCDkey]['Mirr_pos'].copy()
        Exptime = DataDict[CCDkey]['Exptime'].copy()
        
        label = np.zeros(Nobs,dtype='40str')
        
        uMirr_pos = np.sort(np.unique(Mirr_pos))
        

        for ixMP,iMP in enumerate(uMirr_pos):
            
            ixselMP = Mirr_pos == iMP
            
            label[ixselMP & (Exptime > 0)] = 'focus_%i' % ixMP
            label[ixselMP & (Exptime ==0)] = 'BGD'
            
        
        DataDict[CCDkey]['label'] = label.copy()

        
        rootFile_name = DataDict[CCDkey]['File_name'].copy()
        
        File_name  = ['%s.fits' % item for item in rootFile_name]
        
        DataDict[CCDkey]['Files'] = np.array(File_name).copy()
        
    
    return DataDict, isconsistent


def get_basic_spot_FOCUS00(stamp,x0,y0,log=None,debug=False):
    """ 
    # TODO:
    #   get basic statistics, measure and subtract background
    #   update centroid
    #   do aperture photometry
    #   pack-up results and return
    """
    
    rap = 5
    rin = 10
    rout = max(stamp.shape)
    gain = 3.5 # used to compute photom. error

    if debug:
        return dict(bgd=0.,peak=1.,fluence=1.,efluence=0.,
               x=1.,y=1.,fwhmx=1.,fwhmy=1.,stamp=stamp.copy())

    spot = Spot(stamp,log)
    spot.xcen = x0
    spot.ycen = y0
    
    bgd,sbgd = spot.measure_bgd(rin,rout)
    spot.sub_bgd(rin,rout)
    peak = spot.data.max()
    xcen,ycen,fwhmx,fwhmy = spot.get_centroid(full=True)
    centre = spot.xcen,spot.ycen

    flu,eflu = spot.doap_photom(centre,rap,rin,rout,gain=gain,doErrors=True,
                    subbgd=True)
    
    x,y = spot.xcen,spot.ycen
    res = dict(bgd=bgd,peak=peak,fluence=flu,efluence=eflu,
               x=x,y=y,fwhmx=fwhmx,fwhmy=fwhmy,stamp=spot.data.copy())
         
    return res


def get_shape_spot_FOCUS00(stamp,x0,y0,method='G',log=None,debug=False):
    """ 
    # TODO:
    #   get basic statistics, measure and subtract background
    #   update centroid
    #   do aperture photometry
    #   pack-up results and return
    """
    
    if debug and method == 'G':
        return dict()
    
    if debug and method == 'M':
        return dict()

    spot = Spot(stamp,log)
    spot.xcen = x0
    spot.ycen = y0
    
    if method == 'G':
        # gauss(): i00, xcen, self.ycen, x_stddev=0.5, y_stddev=0.5
        # G_keys = ['x','ey','y','ey','i0','ei0','sigma_x','esigma_x',
        # 'sigma_y','esigma_y','fwhm_x','efwhm_x','fwhm_y','efwhm_y',
        # 'fluence','efluence']
    
        Gpars, eGpars = spot.fit_Gauss()
        
        res = dict(i0=Gpars[0],ei0=eGpars[0],
                  x=Gpars[1],ex=eGpars[1],y=Gpars[2],ey=eGpars[2],
                  sigma_x=Gpars[3],esigma_x=eGpars[3],
                  sigma_y=Gpars[4],esigma_y=eGpars[4])
        res['fwhm_x'] = res['sigma_x']*2.355
        res['efwhm_x'] = res['esigma_x']*2.355
        res['fwhm_y'] = res['sigma_y']*2.355
        res['efwhm_y'] = res['esigma_y']*2.355
        res['fluence'] = None  # ??
        res['efluence'] = None # ??
        
        
    if method == 'M':
        #res = dict(centreX=quad['centreX']+1, centreY=quad['centreY']+1,
        #           e1=quad['e1'], e2=quad['e2'],
        #           ellipticity=quad['ellipticity'],
        #           R2=R2,
        #           R2arcsec=R2arcsec,
        #           GaussianWeighted=GaussianWeighted,
        #           a=quad['a'], b=quad['b'])
        #M_keys = ['x','y','ellip','e1','e2','R2','a','b']
        rawres = spot.measureRefinedEllipticity()
        res = dict(x=rawres['centreX']-1.,y=rawres['centreY']-1.,
                   ellip=rawres['ellipticity'],e1=rawres['e1'],
                   e2=rawres['e2'],R2=rawres['R2'],
                   a=rawres['a'],b=rawres['b'])

         
    return res    

    
def prep_data_FOCUS00(DataDict,Report,inputs,log=None,debug=False):
    """Takes Raw Data and prepares it for further analysis. 
    Also checks that data quality is enough."""
    
    # Inputs un-packing
    
    FFs = inputs['FFs'] # flat-fields for each CCD (1,2,3)
    resultspath = inputs['resultspath']
    
    # Load Flat-Field data for each CCD
    
    FFdata = dict()
    for CCDindex in range(1,4):
        CCDkey = 'CCD%i' % CCDindex
        if CCDkey in FFs.keys():
            FFdata[CCDkey] = FFing.FlatField(fitsfile=FFs[CCDkey])

    Report.add_Section(Title='Data Pre-Processing and QA',level=0)
    
    # Initialization of columns to be added to DataDict
    
    
    ccd_keys = ['spots_file']
    
    ccd_keys_Q = ['offset','med_pre','std_pre','med_ove','std_ove','med_img',
                   'std_img']
    ccd_formats = dict(spots_file='S30')
    
    spot_basic_keys = ['x','y','peak','fluence', 'efluence','bgd','fwhmy', 'fwhmx']
    
    coo_dict_CCDQ = polib.Point_CooNom['CCD1']['E']
    spotIDs = polib.Point_CooNom['names']
    
    
    for CCDindex in range(1,4):
        
        CCDkey = 'CCD%i' % CCDindex
        
        if CCDkey in DataDict:
            
            # Loop over ObsIDs
            
            Nobs = len(DataDict[CCDkey]['ObsID'].copy())
            
            for key in ccd_keys:
                try: dtype = ccd_formats[key]
                except: dtype = 'float32'
                DataDict[CCDkey][key] = (np.zeros(Nobs) + np.nan).astype(dtype)
            
            for key in ccd_keys_Q:
                try: dtype = ccd_formats[key]
                except: dtype = 'float32'
                for Q in pilib.Quads:
                    DataDict[CCDkey]['%s_%s' % (key,Q)] = (np.zeros(Nobs) + np.nan).astype(dtype)

            for key in spot_basic_keys:
                dtype = 'float32'
                for Q in pilib.Quads:
                    for spotID in spotIDs:
                        DataDict[CCDkey]['bas_%s_%s_%s' % (key,Q,spotID)] = \
                         (np.zeros(Nobs) + np.nan).astype(dtype)

    # Loop over CCDs
    
    for CCDindex in range(1,4):
        
        CCDkey = 'CCD%i' % CCDindex
        
        if CCDkey in DataDict:
            
            # Loop over ObsIDs
            
            ObsIDs = DataDict[CCDkey]['ObsID'].copy()
            label = DataDict[CCDkey]['label'].copy()
            

            for iObs, ObsID in enumerate(ObsIDs):
                
                ilabel = label[iObs]
                
                # log object-id being analyzed: ObsID
                
                if log is not None: log.info('working on ObsID %i' % ObsID)
                
                # retrieve FITS file and open it
                                
                idatapath = DataDict[CCDkey]['datapath'][iObs]
                
                fitsf = os.path.join(idatapath,DataDict[CCDkey]['Files'][iObs])
                
                
                # subtract offset and add to DataDict
                # measure basic image statistics and add to DataDict
                                
                ccdobj = ccd.CCD(fitsf)
                
                
                for Q in pilib.Quads:
                    
                    Qoffset = ccdobj.sub_offset(Q,method='median',scan='pre',trimscan=[5,5])[0]
                    DataDict[CCDkey]['offset_%s' % Q][iObs] = Qoffset
                    
                    med_pre,std_pre = ccdobj.get_stats(Q,sector='pre',statkeys=['median','std'],trimscan=[5,5])
                    med_ove,std_ove = ccdobj.get_stats(Q,sector='ove',statkeys=['median','std'],trimscan=[5,5])
                    med_img,std_img = ccdobj.get_stats(Q,sector='img',statkeys=['median','std'],trimscan=[10,10])
                    
                    DataDict[CCDkey]['med_pre_%s' % Q][iObs] = med_pre
                    DataDict[CCDkey]['std_pre_%s' % Q][iObs] = std_pre
                    DataDict[CCDkey]['med_ove_%s' % Q][iObs] = med_ove
                    DataDict[CCDkey]['std_ove_%s' % Q][iObs] = std_ove
                    DataDict[CCDkey]['med_img_%s' % Q][iObs] = med_img
                    DataDict[CCDkey]['std_img_%s' % Q][iObs] = std_img
                    
                
                if ilabel != 'BGD':
                    
                    spots_bag = dict()
                    
                    
                    # Divide by flat-field
    
                    if not debug: 
                        ccdobj.divide_by_flatfield(FFdata[CCDkey].Flat)
                        if log is not None: log.info('Divided Image by Flat-field')

                    
                    for Q in pilib.Quads:
                        
                        spots_bag[Q] = dict()
                        
                        coo_dict_CCDQ = polib.Point_CooNom[CCDkey][Q]
                        spotIDs = coo_dict_CCDQ.keys()
                        B=ccdobj.QuadBound[Q]
                    
                        for spotID in spotIDs:
                            
                            if log is not None: log.info('ObsID - CCD - spotID = %s-%s' % (ObsID,CCDkey,spotID))
                            
                            # get coordinates of spotID
                            
                            coo = coo_dict_CCDQ[spotID]
                            
                            x0Q,y0Q = tuple([int(np.round(item)) for item in coo])
                            x0Q -= B[0]
                            y0Q -= B[2]
                            
                            corners = [x0Q-stampw/2,x0Q-stampw/2+stampw,
                                       y0Q-stampw/2,y0Q-stampw/2+stampw]
                            
                            # Cut-out stamp of the spot
                            stamp = ccdobj.get_cutout(corners,Q,canonical=False)
                            
                            x0 = stampw/2
                            y0 = stampw/2

                            # do basic measurements on each spot and add to DataDict:
                            # spot_basic.keys() 
                            # ['x,','y','peak','fluence', 'efluence','bgd', 
                            # 'fwhmy', 'fwhmx']
                            
                            spot_basic = get_basic_spot_FOCUS00(stamp,x0,y0,debug=debug)
                            stamp = spot_basic['stamp'].copy()
                            
                            for key in spot_basic_keys:
                                DataDict[CCDkey]['bas_%s_%s_%s' % (key,Q,spotID)][iObs] =\
                                    spot_basic[key]
                            
                            spots_bag[Q][spotID] = spot_basic.copy()
                        
                    # save spot-data to a hard-file and add path to DataDict
                        
                    spots_bag_f = 'SPOTS_%i_%s.pick' % (ObsID,CCDkey)
                    DataDict[CCDkey]['spots_file'][iObs] = spots_bag_f
                    spots_bag_f = os.path.join(resultspath,spots_bag_f)
                    files.cPickleDumpDictionary(spots_bag,spots_bag_f)
    
    
    
    for CCDindex in range(1,4):
        if 'CCD%i' % CCDindex in DataDict:
            CCDindex_ref = CCDindex
    ObsIDs = DataDict['CCD%i' % CCDindex_ref]['ObsID'].copy()
    label = DataDict['CCD%i' % CCDindex_ref]['label'].copy()
    
    for iObs,ObsID in enumerate(ObsIDs):
        
        ilabel = label[iObs]
        
        if ilabel == 'BGD': continue
        
        spots_bag_FP = dict()
        
        for CCDindex in range(1,4):
            
            spots_bag_f = os.path.join(resultspath,DataDict[CCDkey]['spots_file'][iObs])
            spots_bag_FP['CCD%i' % CCDindex] = files.cPickleRead(spots_bag_f)
        
        spots_disp_f = os.path.join(resultspath,'STAMPS_spots_%i.png' % (ObsID,))
        spots_disp_title = 'ObsID - %i' % (ObsID,)
        pdspl.show_spots_allCCDs(spots_bag_FP,title=spots_disp_title,
                                 filename=spots_disp_f,dobar=True)
    

    # Data Quality Assessment:
    # TODO:
    #plotF00_CCDkey_vstime(DataDict,'offset')
    #plotF00_CCDkey_vstime(DataDict,'std_pre')
         
    #spot_basic_keys = ['x','y','peak','fluence', 'efluence','bgd','fwhmy', 'fwhmx']
    
    #TODO:
#    for CCDindex in range(1,4):
#        
#        CCDkey = 'CCD%i' % CCDindex
#        
#        if CCDkey not in DataDict: continue
#    
#        plotF00_spotkey_vstime(DataDict,'peak',ccdkey=CCDkey,filename='')
#        plotF00_spotkey_vstime(DataDict,'fluence',error='efluence',
#                               ccdkey=CCDkey,filename='')
#        
#        plotF00_spotkey_vstime(DataDict,'fwhmx',ccdkey=CCDkey,filename='')
#        plotF00_spotkey_vstime(DataDict,'fwhmy',ccdkey=CCDkey,filename='')
    
    
    # TODO:
    # Check all parameters are within expected ranges:
    #    offsets, STDs, peak fluence, fwhm, ellipticity
    #    save reports to RepDict
    
    
    
    return DataDict, Report


def basic_analysis_FOCUS00(DataDict,Report,inputs,log=None,debug=False):
    """Performs basic analysis on spots:
         - 2D Gaussian Model shape measurements
         - Quadrupole Moments shape measurements
    """
    
    # Inputs un-packing
    
    resultspath = inputs['resultspath']

    
    Report.add_Section(Title='Data Analysis: Basic',level=0)
    
    # Initialization of columns to be added to DataDict
    
    Quadrants = pilib.Quads
    
    G_keys = ['x','ex','y','ey','i0','ei0','sigma_x','esigma_x',
    'sigma_y','esigma_y','fwhm_x','efwhm_x','fwhm_y','efwhm_y',
    'fluence','efluence']
    
    M_keys = ['x','y','ellip','e1','e2','R2','a','b']
    
    
    spotIDs = polib.Point_CooNom['names']
    
    for CCDindex in range(1,4):
        
        CCDkey = 'CCD%i' % CCDindex
        
        if CCDkey in DataDict:
            
            # Loop over ObsIDs
            
            Nobs = len(DataDict[CCDkey]['ObsID'].copy())
            

            for key in G_keys:
                dtype = 'float32'
                for Q in pilib.Quads:
                    for spotID in spotIDs:
                        DataDict[CCDkey]['G_%s_%s_%s' % (key,Q,spotID)] = \
                         (np.zeros(Nobs) + np.nan).astype(dtype)
                         
            for key in M_keys:
                dtype = 'float32'
                for Q in pilib.Quads:
                    for spotID in spotIDs:
                        DataDict[CCDkey]['M_%s_%s_%s' % (key,Q,spotID)] = \
                         (np.zeros(Nobs) + np.nan).astype(dtype)

    # Loop over CCDs, Obsids, quadrants and spots
    
    for CCDindex in range(1,4):
        if 'CCD%i' % CCDindex in DataDict:
            CCDindex_ref = CCDindex
    ObsIDs = DataDict['CCD%i' % CCDindex_ref]['ObsID'].copy()
    label = DataDict['CCD%i' % CCDindex_ref]['label'].copy()
    
    for iObs,ObsID in enumerate(ObsIDs):
        
        ilabel = label[iObs]
        
        if ilabel == 'BGD': continue
        
        
        for CCDindex in range(1,4):
            
            spots_bag_f = os.path.join(resultspath,DataDict[CCDkey]['spots_file'][iObs])
            spots_bag = files.cPickleRead(spots_bag_f)
            
            for Q in Quadrants:
                for spotID in spotIDs:
                    data = spots_bag[Q][spotID]
                    
                    stamp = data['stamp'].copy()
                    x0,y0 = data['x'],data['y']
                    
                    spot_shape_G = get_shape_spot_FOCUS00(stamp,x0,y0,method='G',
                                                          log=log,debug=debug)
                    
                    for key in spot_shape_G:
                        DataDict[CCDkey]['G_%s_%s_%s' % (key,Q,spotID)][iObs] =\
                           spot_shape_G[key]
                    
                    spot_shape_M = get_shape_spot_FOCUS00(stamp,x0,y0,method='M',
                                                          log=log,debug=debug)
                    
                    for key in spot_shape_M:
                        DataDict[CCDkey]['M_%s_%s_%s' % (key,Q,spotID)][iObs] =\
                           spot_shape_M[key]
    
    # QA

    return DataDict, Report
             


def meta_analysis_FOCUS00(DataDict,RepDict,inputs,log=None):
    """
    
    Analyzes the relation between PSF shape and mirror position.
    
    """


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
    'SerRdDel':0,'Trappump':'x','TP_Ser_S':0,
    'TP_Ver_S':0,'TP_DW_V':0,'TP_DW_H':0,'TOI_flsh':143,'TOI_pump':1000,
    'TOI_read':1000,'TOI_CInj':1000,'Invflshp':500,'Invflush':1,
    'Flushes':3,'Vstart':1,'Vend':2066,'Ovrscn_H':0,'CLK_ROE':'Normal',
    'CnvStart':0,'SumWell':0,'IniSweep':1,'SPW_clk':1,'FPGA_ver':'x.x.x',
    'EGSE_ver':elvis,'M_Steps':0,'M_St_Sze':0,
    'Wavelength':wavelength,'Chmb_pre':1.E-6,'Mirr_pos':0,
    'RPSU_SN':'RP01','ROE_SN':'RO01','CalScrpt':'FakeScriptFOCUS00',
    'R1CCD1TT':153,'R1CCD1TB':153,'R1CCD2TT':153,'R1CCD2TB':153,'R1CCD3TT':153,
    'R1CCD3TB':153,'IDL_V':13,'IDH_V':18,'IG1_T1_V':4,
    'IG1_T2_V':6,'IG1_T3_V':4,'IG1_B1_V':4,'IG1_B2_V':6,'IG1_B3_V':4,
    'IG2_T_V':6,'IG2_B_V':6,'OD_T1_V':26,'OD_T2_V':26,'OD_T3_V':26,
    'OD_B1_V':26,'OD_B2_V':26,'OD_B3_V':26,'RD_T_V':17,'RD_B_V':17}
    
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
    
    doneObsids = []
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        if obsid in doneObsids: continue # to avoid duplications 
                                     # (each CCD has an entry in explog, so 3 entries per OBSID)
        
        idate = explog['DATE'][ixobs]
        
        idtobj = pilib.get_dtobj(idate)
        
        if ixobs < Nobs-1:
            ip1dtobj = pilib.get_dtobj(explog['DATE'][ixobs+1])
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
        
        
        HKfile.write(HKfilef,format='ascii',overwrite=True,delimiter='\t')
        
        doneObsids.append(obsid)
        

def generate_FITS_FOCUS00(wavelength,explog,datapath,elvis='6.0.0'):
    """ """
    
    NAXIS1,NAXIS2 = 4238,4132
    
    maxexptime = explog['Exptime'].max()
    flatlevel = 200.
    biaslevel = 2000.
    pflux = 1.E4 # adu
    pfwhm = 9./12. # pixels
    
    FilterID = pilib.get_FW_ID(wavelength)
    
    mirror_nom = polib.mirror_nom[FilterID]
    
    waivedkeys = ['File_name','Flsh-Rdout_e_time','C.Inj-Rdout_e_time',
                  'Wavelength']
    
    for ixobs,obsid in enumerate(explog['ObsID']):
        
        
        idate = explog['DATE'][ixobs]
        iCCD = explog['CCD'][ixobs]
        iexptime = explog['Exptime'][ixobs]
        iMirr_pos = explog['Mirr_pos'][ixobs]
        
        #if iexptime == 0.: continue # TESTS
        
        idtobj = pilib.get_dtobj(idate)
        
        dmy = idtobj.strftime('%d%m%y')
        HMS = idtobj.strftime('%H%M%S')

        FITSf = 'EUC_%s_%sD_%sT_ROE1_%s.fits' % \
            (obsid,dmy,HMS,iCCD)
        
        FITSf = os.path.join(datapath,FITSf)
        
        ccdobj = ccd.CCD()
        
        img = np.zeros(shape=(NAXIS1,NAXIS2),dtype='float32')
        
        ccdobj.add_extension(data=None)
        ccdobj.add_extension(data=img,label='ROE1_%s' % iCCD)
        
        ccdobj.simadd_flatilum(levels=dict(E=flatlevel*iexptime/maxexptime*1.,
                                           F=flatlevel*iexptime/maxexptime*1.1,
                                           G=flatlevel*iexptime/maxexptime*1.2,
                                           H=flatlevel*iexptime/maxexptime*1.3))
        
        if iexptime > 0:
            
            ipflux = pflux * iexptime/maxexptime
            ipfwhm = pfwhm * (1.+((iMirr_pos-mirror_nom)/0.2)**2.)
            
            ccdobj.simadd_points(ipflux,ipfwhm,CCDID=iCCD,dx=0,dy=0)
        
        
        ccdobj.simadd_poisson()
        
        ccdobj.simadd_bias(levels=dict(E=biaslevel*1.,
                                           F=biaslevel*1.1,
                                           G=biaslevel*1.2,
                                           H=biaslevel*1.3))
        ccdobj.simadd_ron()
        
        ccdobj.extensions[-1].data = np.round(ccdobj.extensions[-1].data).astype('int32')

        ccdobj.extensions[-1].header['WAVELENG'] = explog['Wavelength'][ixobs]
        
        for key in ELtools.columnlist[elvis]:
            if key not in waivedkeys:
                ccdobj.extensions[-1].header[key] = explog[key][ixobs]
        
        ccdobj.writeto(FITSf,clobber=True,unsigned16bit=True)
        
    


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
    
    explog.write(explogf,format='ascii',overwrite=True,delimiter='\t')
    
    generate_HK_FOCUS00(explog,datapath,elvis='6.0.0')
    
    generate_FITS_FOCUS00(wavelength,explog,datapath,elvis='6.0.0')
    
    
    
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
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
        
    else:
        
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
        
    # DATA-WORK
    
    # Prepare Data for further analysis (subtract offsets, divide by FFF, trim snapshots). 
    # Check Data has enough quality:
    #     median levels in pre-scan, image-area, overscan
    #     fluences and spot-sizes (coarse measure) match expectations for all spots
    
    debug = False
    
    if todo_flags['prep']:
        DataDict, reportobj = prep_data_FOCUS00(DataDict,reportobj,inputs,log,debug)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)        
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
    
    # Optional
    # Perform Basic Analysis : Gaussian fits and Moments shape measurements of spots
    
    if todo_flags['basic']:
        DataDict, reportobj = basic_analysis_FOCUS00(DataDict,reportobj,inputs,log)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)        
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
    
    stop()
    
    # Optional
    # Produce Summary Figures and Tables
    
    if todo_flags['meta']:
        DataDict, reportobj = meta_analysis_FOCUS00(DataDict,reportobj,inputs,log)
        pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)  
    else:
        DataDict, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
    
    # Write automatic Report of Results
    
    if todo_flags['report']:
        
        reportobj.doreport(reportroot,cleanafter)
        outfiles = reportobj.writeto(reportroot,cleanafter)
        
        for outfile in outfiles:
            os.system('mv %s %s/' % (outfile,resultspath))
    
    pilib.save_progress(DataDict,reportobj,DataDictFile,reportobjFile)
    
    if log is not None:
        log.info('Finished FOCUS00')


if __name__ == '__main__':
    
    pass