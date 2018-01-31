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
import copy
from collections import OrderedDict

#from vison.pipe import lib as pilib
from vison.support import context
from vison.ogse import ogse
from vison.flat import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
from vison.datamodel import scriptic as sc
from vison.point.spot import Spot
from vison.point import display as pdspl
from vison.support import vistime
#from vison.pipe.task import Task
import PointTask as PT
#from PointTask import PointTask, BGD_lims
import FOCUS00_lib as F00lib
from vison.image import performance
from vison.datamodel import inputs
import F00aux
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

dtobj_default = vistime.dtobj_default

stampw = 25

FOCUS00_commvalues = dict(program='CALCAMP',test='FOCUS_%i',
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,
  shuttr=1,
  vstart=0,vend=2086,
  siflsh=0,siflsh_p=500,
  mirr_on=1,
  motr_on=0,
  source='point')


FWHM_lims = OrderedDict(CCD1=OrderedDict(E=OrderedDict(ALPHA=[0.5,10.]))) # CCD-Q-Spot, pixels
for Spotname in polib.Point_CooNom['names'][1:]: FWHM_lims['CCD1']['E'][Spotname] = copy.deepcopy(FWHM_lims['CCD1']['E']['ALPHA'])
for Q in ['F','G','H']: FWHM_lims['CCD1'][Q] = copy.deepcopy(FWHM_lims['CCD1']['E'])
for iCCD in [2,3]: FWHM_lims['CCD%i' % iCCD] = copy.deepcopy(FWHM_lims['CCD1'])

Flu_lims = OrderedDict(CCD1=OrderedDict(E=OrderedDict(ALPHA=[2.*40./100.*2.**16,2.*75./100.*2.**16]))) # CCD-Q-Spot
for Spotname in polib.Point_CooNom['names'][1:]: Flu_lims['CCD1']['E'][Spotname] = copy.deepcopy(Flu_lims['CCD1']['E']['ALPHA'])
for Q in ['F','G','H']: Flu_lims['CCD1'][Q] = copy.deepcopy(Flu_lims['CCD1']['E'])
for iCCD in [2,3]: Flu_lims['CCD%i' % iCCD] = copy.deepcopy(Flu_lims['CCD1'])

class FOCUS00_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('exptime',([float],'Exposure time.')),
            ('wavelength',([int],'Wavelength')),
            ])))

class FOCUS00(PT.PointTask):
    """ """
    
    inputsclass = FOCUS00_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(FOCUS00,self).__init__(inputs,log,drill,debug)
        self.name = 'FOCUS00'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = F00aux.gt_F00figs(self.inputs['wavelength'])
        self.inputs['subpaths'] = dict(figs='figs')   
        


    def set_inpdefaults(self,**kwargs):
        
        self.inpdefaults = dict(wavelength=800,
                                exptime=60./100.*ogse.tFWC_point['nm%i' % 800])
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        self.perfdefaults['BGD_lims'] = copy.deepcopy(PT.BGD_lims)
        self.perfdefaults['FWHM_lims'] = copy.deepcopy(FWHM_lims)
        self.perfdefaults['Flu_lims'] = copy.deepcopy(Flu_lims)
         

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """Builds FOCUS00 script structure dictionary.
        
        #:param wavelength: int, [nm], wavelength.
        #:param exptime: int, [ms], exposure time.
        :param diffvalues: dict, opt, differential values.
        
        
        """
        
        wavelength = self.inputs['wavelength']
        exptime = self.inputs['exptime']
        
            
        FW_ID = ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])
        mirror_nom = polib.mirror_nom[FW_ID]
        
        #FOCUS00_sdict = dict(col1=dict(frames=5,wave=FW_IDX,exptime=0,
        #                               mirr_pos=mirror_nom-5,
        #                               comments='BGD'))
        
        FOCUS00_sdict = dict()
        
        for i,j in enumerate(range(-3,4,1)):
            FOCUS00_sdict['col%i' % (i+1,)] = dict(frames=2,
                          test='FOCUS00_%i' % wavelength,
                          exptime=exptime,
                          mirr_pos=mirror_nom+float(j)*1.,
                          wave=FW_IDX,
                          comments='F%.1f' % float(j))
        
        Ncols = len(FOCUS00_sdict.keys())    
        FOCUS00_sdict['Ncols'] = Ncols
        
        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(FOCUS00_commvalues)           
        
        if len(diffvalues)==0:
            diffvalues = self.inputs['diffvalues']
        
        FOCUS00_sdict = sc.update_structdict(FOCUS00_sdict,commvalues,diffvalues)
        
        
        return FOCUS00_sdict

    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(FOCUS00,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
                              wavedkeys=wavedkeys)
        #return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
        #                      wavedkeys=wavedkeys,elvis=elvis)
    
    
    def prep_data(self):
        """ """
        raise NotImplementedError
    
    def basic_analysis(self):
        """ """
        raise NotImplementedError
    
    def meta_analysis(self):
        """ """
        raise NotImplementedError


#def get_basic_spot_FOCUS00(stamp,x0,y0,log=None,debug=False):
#    """ 
#    # TODO:
#    #   get basic statistics, measure and subtract background
#    #   update centroid
#    #   do aperture photometry
#    #   pack-up results and return
#    """
#    
#    rap = 5
#    rin = 10
#    rout = max(stamp.shape)
#    gain = 3.5 # used to compute photom. error
#
#    if debug:
#        return dict(bgd=0.,peak=1.,fluence=1.,efluence=0.,
#               x=1.,y=1.,fwhmx=1.,fwhmy=1.,stamp=stamp.copy())
#
#    spot = Spot(stamp,log)
#    spot.xcen = x0
#    spot.ycen = y0
#    
#    bgd,sbgd = spot.measure_bgd(rin,rout)
#    spot.sub_bgd(rin,rout)
#    peak = spot.data.max()
#    xcen,ycen,fwhmx,fwhmy = spot.get_centroid(full=True)
#    centre = spot.xcen,spot.ycen
#
#    flu,eflu = spot.doap_photom(centre,rap,rin,rout,gain=gain,doErrors=True,
#                    subbgd=True)
#    
#    x,y = spot.xcen,spot.ycen
#    res = dict(bgd=bgd,peak=peak,fluence=flu,efluence=eflu,
#               x=x,y=y,fwhmx=fwhmx,fwhmy=fwhmy,stamp=spot.data.copy())
#         
#    return res
#
#
#def get_shape_spot_FOCUS00(stamp,x0,y0,method='G',log=None,debug=False):
#    """ 
#    # TODO:
#    #   get basic statistics, measure and subtract background
#    #   update centroid
#    #   do aperture photometry
#    #   pack-up results and return
#    """
#    
#    if debug and method == 'G':
#        return dict()
#    
#    if debug and method == 'M':
#        return dict()
#
#    spot = Spot(stamp,log)
#    spot.xcen = x0
#    spot.ycen = y0
#    
#    if method == 'G':
#        # gauss(): i00, xcen, self.ycen, x_stddev=0.5, y_stddev=0.5
#        # G_keys = ['x','ey','y','ey','i0','ei0','sigma_x','esigma_x',
#        # 'sigma_y','esigma_y','fwhm_x','efwhm_x','fwhm_y','efwhm_y',
#        # 'fluence','efluence']
#    
#        Gpars, eGpars = spot.fit_Gauss()
#        
#        res = dict(i0=Gpars[0],ei0=eGpars[0],
#                  x=Gpars[1],ex=eGpars[1],y=Gpars[2],ey=eGpars[2],
#                  sigma_x=Gpars[3],esigma_x=eGpars[3],
#                  sigma_y=Gpars[4],esigma_y=eGpars[4])
#        res['fwhm_x'] = res['sigma_x']*2.355
#        res['efwhm_x'] = res['esigma_x']*2.355
#        res['fwhm_y'] = res['sigma_y']*2.355
#        res['efwhm_y'] = res['esigma_y']*2.355
#        sigma_maj = np.max([res['sigma_x'],res['sigma_y']])
#        sigma_min = np.min([res['sigma_x'],res['sigma_y']])#
#        q = sigma_min/sigma_maj
#        
#        res['fluence'] = 2.*np.pi*res['i0']*q*sigma_maj**2.
#        res['efluence'] = 2.*np.pi*res['ei0']*q*sigma_maj**2.
#        
#        
#    if method == 'M':
#        #res = dict(centreX=quad['centreX']+1, centreY=quad['centreY']+1,
#        #           e1=quad['e1'], e2=quad['e2'],
#        #           ellipticity=quad['ellipticity'],
#        #           R2=R2,
#        #           R2arcsec=R2arcsec,
#        #           GaussianWeighted=GaussianWeighted,
#        #           a=quad['a'], b=quad['b'])
#        #M_keys = ['x','y','ellip','e1','e2','R2','a','b']
#        rawres = spot.measureRefinedEllipticity()
#        res = dict(x=rawres['centreX']-1.,y=rawres['centreY']-1.,
#                   ellip=rawres['ellipticity'],e1=rawres['e1'],
#                   e2=rawres['e2'],R2=rawres['R2'],
#                   a=rawres['a'],b=rawres['b'])
#
#    
#    return res    
#
#    
#def prep_data_FOCUS00(DataDict,Report,inputs,log=None,debug=False):
#    """Takes Raw Data and prepares it for further analysis. 
#    Also checks that data quality is enough."""
#    
#    # Inputs un-packing
#    
#    FFs = inputs['FFs'] # flat-fields for each CCD (1,2,3)
#    resultspath = inputs['resultspath']
#    
#    # Load Flat-Field data for each CCD
#    
#    FFdata = dict()
#    for CCDindex in range(1,4):
#        CCDkey = 'CCD%i' % CCDindex
#        if CCDkey in FFs.keys():
#            FFdata[CCDkey] = FFing.FlatField(fitsfile=FFs[CCDkey])
#    
#    ThisSectionTitle = 'Data Pre-Processing and QA'
#    try: 
#        lastsection = Report.Contents[-1]
#        if (lastsection.type == 'section') and \
#           (lastsection.Title == ThisSectionTitle):
#               _ = Report.Contents.pop(-1)
#        Report.add_Section(Title=ThisSectionTitle,level=0)
#    except:
#        pass
#    
#    Report.add_Section(Title='Data Pre-Processing and QA',level=0)
#    
#    # Initialization of columns to be added to DataDict
#    
#    
#    ccd_keys = ['spots_file']
#    
#    ccd_keys_Q = ['offset','med_pre','std_pre','med_ove','std_ove','med_img',
#                   'std_img']
#    ccd_formats = dict(spots_file='S30')
#    
#    spot_basic_keys = ['x','y','peak','fluence', 'efluence','bgd','fwhmy', 'fwhmx']
#    
#    coo_dict_CCDQ = polib.Point_CooNom['CCD1']['E']
#    spotIDs = polib.Point_CooNom['names']
#    
#    
#    for CCDindex in range(1,4):
#        
#        CCDkey = 'CCD%i' % CCDindex
#        
#        if CCDkey in DataDict:
#            
#            # Loop over ObsIDs
#            
#            Nobs = len(DataDict[CCDkey]['ObsID'].copy())
#            
#            for key in ccd_keys:
#                try: dtype = ccd_formats[key]
#                except: dtype = 'float32'
#                DataDict[CCDkey][key] = (np.zeros(Nobs) + np.nan).astype(dtype)
#            
#            for key in ccd_keys_Q:
#                try: dtype = ccd_formats[key]
#                except: dtype = 'float32'
#                for Q in pilib.Quads:
#                    DataDict[CCDkey]['%s_%s' % (key,Q)] = (np.zeros(Nobs) + np.nan).astype(dtype)
#
#            for key in spot_basic_keys:
#                dtype = 'float32'
#                for Q in pilib.Quads:
#                    for spotID in spotIDs:
#                        DataDict[CCDkey]['bas_%s_%s_%s' % (key,Q,spotID)] = \
#                         (np.zeros(Nobs) + np.nan).astype(dtype)
#
#    # Loop over CCDs
#    
#    for CCDindex in range(1,4):
#        
#        CCDkey = 'CCD%i' % CCDindex
#        
#        if CCDkey in DataDict:
#            
#            # Loop over ObsIDs
#            
#            ObsIDs = DataDict[CCDkey]['ObsID'].copy()
#            label = DataDict[CCDkey]['label'].copy()
#            
#
#            for iObs, ObsID in enumerate(ObsIDs):
#                
#                ilabel = label[iObs]
#                
#                # log object-id being analyzed: ObsID
#                
#                if log is not None: log.info('working on ObsID %i' % ObsID)
#                
#                # retrieve FITS file and open it
#                                
#                idatapath = DataDict[CCDkey]['datapath'][iObs]
#                
#                fitsf = os.path.join(idatapath,DataDict[CCDkey]['Files'][iObs])
#                
#                
#                # subtract offset and add to DataDict
#                # measure basic image statistics and add to DataDict
#                                
#                ccdobj = ccd.CCD(fitsf)
#                
#                
#                for Q in pilib.Quads:
#                    
#                    Qoffset = ccdobj.sub_offset(Q,method='median',scan='pre',trimscan=[5,5])[0]
#                    DataDict[CCDkey]['offset_%s' % Q][iObs] = Qoffset
#                    
#                    med_pre,std_pre = ccdobj.get_stats(Q,sector='pre',statkeys=['median','std'],trimscan=[5,5])
#                    med_ove,std_ove = ccdobj.get_stats(Q,sector='ove',statkeys=['median','std'],trimscan=[5,5])
#                    med_img,std_img = ccdobj.get_stats(Q,sector='img',statkeys=['median','std'],trimscan=[10,10])
#                    
#                    DataDict[CCDkey]['med_pre_%s' % Q][iObs] = med_pre
#                    DataDict[CCDkey]['std_pre_%s' % Q][iObs] = std_pre
#                    DataDict[CCDkey]['med_ove_%s' % Q][iObs] = med_ove
#                    DataDict[CCDkey]['std_ove_%s' % Q][iObs] = std_ove
#                    DataDict[CCDkey]['med_img_%s' % Q][iObs] = med_img
#                    DataDict[CCDkey]['std_img_%s' % Q][iObs] = std_img
#                    
#                
#                if ilabel != 'BGD':
#                    
#                    spots_bag = dict()
#                    
#                    
#                    # Divide by flat-field
#    
#                    if not debug: 
#                        ccdobj.divide_by_flatfield(FFdata[CCDkey].Flat)
#                        if log is not None: log.info('Divided Image by Flat-field')
#
#                    
#                    for Q in pilib.Quads:
#                        
#                        spots_bag[Q] = dict()
#                        
#                        coo_dict_CCDQ = polib.Point_CooNom[CCDkey][Q]
#                        spotIDs = coo_dict_CCDQ.keys()
#                        B=ccdobj.QuadBound[Q]
#                    
#                        for spotID in spotIDs:
#                            
#                            if log is not None: log.info('ObsID - CCD - spotID = %s-%s' % (ObsID,CCDkey,spotID))
#                            
#                            # get coordinates of spotID
#                            
#                            coo = coo_dict_CCDQ[spotID]
#                            
#                            x0Q,y0Q = tuple([int(np.round(item)) for item in coo])
#                            x0Q -= B[0]
#                            y0Q -= B[2]
#                            
#                            corners = [x0Q-stampw/2,x0Q-stampw/2+stampw,
#                                       y0Q-stampw/2,y0Q-stampw/2+stampw]
#                            
#                            # Cut-out stamp of the spot
#                            stamp = ccdobj.get_cutout(corners,Q,canonical=False)
#                            
#                            x0 = stampw/2
#                            y0 = stampw/2
#
#                            # do basic measurements on each spot and add to DataDict:
#                            # spot_basic.keys() 
#                            # ['x,','y','peak','fluence', 'efluence','bgd', 
#                            # 'fwhmy', 'fwhmx']
#                            
#                            spot_basic = get_basic_spot_FOCUS00(stamp,x0,y0,debug=debug)
#                            stamp = spot_basic['stamp'].copy()
#                            
#                            for key in spot_basic_keys:
#                                DataDict[CCDkey]['bas_%s_%s_%s' % (key,Q,spotID)][iObs] =\
#                                    spot_basic[key]
#                            
#                            spots_bag[Q][spotID] = spot_basic.copy()
#                        
#                    # save spot-data to a hard-file and add path to DataDict
#                        
#                    spots_bag_f = 'SPOTS_%i_%s.pick' % (ObsID,CCDkey)
#                    DataDict[CCDkey]['spots_file'][iObs] = spots_bag_f
#                    spots_bag_f = os.path.join(resultspath,spots_bag_f)
#                    files.cPickleDumpDictionary(spots_bag,spots_bag_f)
#    
#    
#    
#    for CCDindex in range(1,4):
#        if 'CCD%i' % CCDindex in DataDict:
#            CCDindex_ref = CCDindex
#    ObsIDs = DataDict['CCD%i' % CCDindex_ref]['ObsID'].copy()
#    label = DataDict['CCD%i' % CCDindex_ref]['label'].copy()
#    
#    for iObs,ObsID in enumerate(ObsIDs):
#        
#        ilabel = label[iObs]
#        
#        if ilabel == 'BGD': continue
#        
#        spots_bag_FP = dict()
#        
#        for CCDindex in range(1,4):
#            
#            spots_bag_f = os.path.join(resultspath,DataDict[CCDkey]['spots_file'][iObs])
#            spots_bag_FP['CCD%i' % CCDindex] = files.cPickleRead(spots_bag_f)
#        
#        spots_disp_f = os.path.join(resultspath,'STAMPS_spots_%i.png' % (ObsID,))
#        spots_disp_title = 'ObsID - %i' % (ObsID,)
#        pdspl.show_spots_allCCDs(spots_bag_FP,title=spots_disp_title,
#                                 filename=spots_disp_f,dobar=True)
#    
#
#    # Data Quality Assessment:
#    # TODO:
#    #plotF00_CCDkey_vstime(DataDict,'offset')
#    #plotF00_CCDkey_vstime(DataDict,'std_pre')
#         
#    #spot_basic_keys = ['x','y','peak','fluence', 'efluence','bgd','fwhmy', 'fwhmx']
#    
#    #TODO:
##    for CCDindex in range(1,4):
##        
##        CCDkey = 'CCD%i' % CCDindex
##        
##        if CCDkey not in DataDict: continue
##    
##        plotF00_spotkey_vstime(DataDict,'peak',ccdkey=CCDkey,filename='')
##        plotF00_spotkey_vstime(DataDict,'fluence',error='efluence',
##                               ccdkey=CCDkey,filename='')
##        
##        plotF00_spotkey_vstime(DataDict,'fwhmx',ccdkey=CCDkey,filename='')
##        plotF00_spotkey_vstime(DataDict,'fwhmy',ccdkey=CCDkey,filename='')
#    
#    
#    # TODO:
#    # Check all parameters are within expected ranges:
#    #    offsets, STDs, peak fluence, fwhm, ellipticity
#    #    save reports to RepDict
#    
#    
#    
#    return DataDict, Report
#
#
#def basic_analysis_FOCUS00(DataDict,Report,inputs,log=None,debug=False):
#    """Performs basic analysis on spots:
#         - 2D Gaussian Model shape measurements
#         - Quadrupole Moments shape measurements
#    """
#    
#    # Inputs un-packing
#    
#    resultspath = inputs['resultspath']
#    
#    
#    ThisSectionTitle = 'Data Analysis: Basic'
#    try: 
#        lastsection = Report.Contents[-1]
#        if (lastsection.type == 'section') and \
#           (lastsection.Title == ThisSectionTitle):
#               _ = Report.Contents.pop(-1)
#        Report.add_Section(Title=ThisSectionTitle,level=0)
#    except:
#        pass
#    
#    # Initialization of columns to be added to DataDict
#    
#    Quadrants = pilib.Quads
#    
#    G_keys = ['x','ex','y','ey','i0','ei0','sigma_x','esigma_x',
#    'sigma_y','esigma_y','fwhm_x','efwhm_x','fwhm_y','efwhm_y',
#    'fluence','efluence']
#    
#    M_keys = ['x','y','ellip','e1','e2','R2','a','b']
#    
#    
#    spotIDs = polib.Point_CooNom['names']
#    
#    for CCDindex in range(1,4):
#        
#        CCDkey = 'CCD%i' % CCDindex
#        
#        if CCDkey in DataDict:
#            
#            # Loop over ObsIDs
#            
#            Nobs = len(DataDict[CCDkey]['ObsID'].copy())
#            
#
#            for key in G_keys:
#                dtype = 'float32'
#                for Q in pilib.Quads:
#                    for spotID in spotIDs:
#                        DataDict[CCDkey]['G_%s_%s_%s' % (key,Q,spotID)] = \
#                         (np.zeros(Nobs) + np.nan).astype(dtype)
#                         
#            for key in M_keys:
#                dtype = 'float32'
#                for Q in pilib.Quads:
#                    for spotID in spotIDs:
#                        DataDict[CCDkey]['M_%s_%s_%s' % (key,Q,spotID)] = \
#                         (np.zeros(Nobs) + np.nan).astype(dtype)
#
#    # Loop over CCDs, Obsids, quadrants and spots
#    
#    for CCDindex in range(1,4):
#        if 'CCD%i' % CCDindex in DataDict:
#            CCDindex_ref = CCDindex
#    ObsIDs = DataDict['CCD%i' % CCDindex_ref]['ObsID'].copy()
#    label = DataDict['CCD%i' % CCDindex_ref]['label'].copy()
#    
#    for iObs,ObsID in enumerate(ObsIDs):
#        ilabel = label[iObs]
#        
#        if ilabel == 'BGD': continue
#        
#        for CCDindex in range(1,4):
#            
#            CCDkey = 'CCD%i' % CCDindex
#            
#            spots_bag_f = os.path.join(resultspath,DataDict[CCDkey]['spots_file'][iObs])
#            spots_bag = files.cPickleRead(spots_bag_f)
#            
#            for Q in Quadrants:
#                for spotID in spotIDs:
#                    data = spots_bag[Q][spotID]
#                    
#                    stamp = data['stamp'].copy()
#                    x0,y0 = data['x'],data['y']
#                    
#                    spot_shape_G = get_shape_spot_FOCUS00(stamp,x0,y0,method='G',
#                                                          log=log,debug=debug)
#                    
#                    for key in G_keys:
#                        DataDict[CCDkey]['G_%s_%s_%s' % (key,Q,spotID)][iObs] =\
#                           spot_shape_G[key]
#                    
#                    spot_shape_M = get_shape_spot_FOCUS00(stamp,x0,y0,method='M',
#                                                          log=log,debug=debug)
#                    
#                    for key in M_keys:
#                        DataDict[CCDkey]['M_%s_%s_%s' % (key,Q,spotID)][iObs] =\
#                           spot_shape_M[key]
#    
#    # QA
#    return DataDict, Report
#             
#
#
#def meta_analysis_FOCUS00(DataDict,Report,inputs,log=None):
#    """
#    
#    Analyzes the relation between PSF shape and mirror position.
#    
#    :param DataDict: Dictionary with input data
#    :param Report: Report Objects
#    :param inputs: Dictionary with inputs
#    
#    """
#    
#    ThisSectionTitle = 'Data Analysis: Meta'
#    try: 
#        lastsection = Report.Contents[-1]
#        if (lastsection.type == 'section') and \
#           (lastsection.Title == ThisSectionTitle):
#               _ = Report.Contents.pop(-1)
#        Report.add_Section(Title=ThisSectionTitle,level=0)
#    except:
#        pass
#    
#    Quadrants = pilib.Quads
#    spotIDs = polib.Point_CooNom['names']
#    
#    G_keys = ['x','y','i0','fwhm_x','fwhm_y','fluence']
#    keys_to_fit = ['fwhm_x','fwhm_y']
#    
#    # Loop over CCDs, Obsids, quadrants and spots
#    
#    for CCDindex in range(1,4):
#        if 'CCD%i' % CCDindex in DataDict:
#            CCDindex_ref = CCDindex
#    
#    #ObsIDs = DataDict['CCD%i' % CCDindex_ref]['ObsID'].copy()
#    label = DataDict['CCD%i' % CCDindex_ref]['label'].copy()
#    selix = label != 'BGD'
#    
#    shape_seq = dict()    
#    
#    for CCDindex in range(1,4):
#        CCDkey = 'CCD%i' % CCDindex
#        shape_seq[CCDkey] = dict()
#        shape_seq[CCDkey]['Mirr_pos'] = DataDict[CCDkey]['Mirr_pos'][selix].copy()
#        
#        for Q in Quadrants:
#            
#            for spotID in spotIDs:
#                
#                for key in G_keys:
#                    spkey = 'G_%s_%s_%s' % (key,Q,spotID)
#                    espkey = 'G_e%s_%s_%s' % (key,Q,spotID)
#                    shape_seq[CCDkey][spkey] = dict()
#                    val = DataDict[CCDkey][spkey][selix].copy()
#                    erval = DataDict[CCDkey][espkey][selix].copy()
#                    shape_seq[CCDkey][spkey]['val'] = val
#                    shape_seq[CCDkey][spkey]['eval'] = erval
#                    
#                    if key in keys_to_fit:
#                        
#                        pfit = F00lib.fit_focus_single(shape_seq[CCDkey]['Mirr_pos'],
#                                                     val,yerror=erval,
#                                                     degree=2,doplot=False)
#                        shape_seq[CCDkey][spkey]['coeffs'] = pfit['coeffs']
#                        shape_seq[CCDkey][spkey]['ecoeffs'] = pfit['ecoeffs']
#                        shape_seq[CCDkey][spkey]['focus'] = pfit['focus']
#        
#            
#    F00lib.fit_focus_all(shape_seq,doplot=True)
#    
#    F00lib.inspect_focus_all(shape_seq,doplot=True)
#    
#
#
#    return DataDict, Report
#
#
#
#    
#    
#    
