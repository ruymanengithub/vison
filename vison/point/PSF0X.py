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
import string as st
import warnings
import copy
from collections import OrderedDict

from vison.datamodel import core
from vison.datamodel import ccd
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.ogse import ogse
#from vison.datamodel import HKtools
from vison.datamodel import scriptic as sc
from vison.flat import FlatFielding as FFing
#from vison.point import lib as polib
#from vison.support.report import Report
from vison.support import files
from vison.image import calibration
#from vison.pipe.task import Task
from PointTask import  PointTask
import PSF0Xaux
from vison.image import performance
from vison.datamodel import inputs
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

testdefaults = dict(waves=[590,640,800,890],
                               exptimes=dict(),
                               frames=[20,15,10,4,3])

for w in testdefaults['waves']:
    testdefaults['exptimes']['nm%i' % w] = (np.array([5.,25.,50.,75.,90.])/100.*ogse.tFWC_point['nm%i' % w]).tolist()

stampw = polib.stampw

satur_fluence = 2.*2.E5/performance.gains['CCD1']
Flu_lims = OrderedDict(
        CCD1=OrderedDict(
                E=OrderedDict(
                        ALPHA = OrderedDict(
                                col1 = np.array([0.0225,0.0275])*satur_fluence, # +/-10%
                                col2 = np.array([0.1125,0.1375])*satur_fluence,
                                col3 = np.array([0.225,0.275])*satur_fluence,
                                col4 = np.array([0.3375,0.4125])*satur_fluence,
                                col5 = np.array([0.405,0.495])*satur_fluence))))
for Spot in ['BRAVO','CHARLIE','DELTA','ECHO']: Flu_lims['CCD1']['E'][Spot] = Flu_lims['CCD1']['E']['ALPHA']
for Q in ['F','G','H']: Flu_lims['CCD1'][Q] = copy.deepcopy(Flu_lims['CCD1']['E'])
for CCD in [2,3]: Flu_lims['CCD%i' % CCD] = copy.deepcopy(Flu_lims['CCD1'])


FWHM_lims = OrderedDict(CCD1=
                        OrderedDict(
                                E=OrderedDict(
                                        ALPHA=[1.2,2.5])))
for Spot in ['BRAVO','CHARLIE','DELTA','ECHO']: FWHM_lims['CCD1']['E'][Spot] = FWHM_lims['CCD1']['E']['ALPHA']
for Q in ['F','G','H']: FWHM_lims['CCD1'][Q] = copy.deepcopy(FWHM_lims['CCD1']['E'])
for CCD in [2,3]: FWHM_lims['CCD%i' % CCD] = copy.deepcopy(FWHM_lims['CCD1'])

BGD_lims = OrderedDict(CCD1=OrderedDict(
                    E=[-1.,5.]))
for Q in ['F','G','H']: BGD_lims['CCD1'][Q] = copy.deepcopy(BGD_lims['CCD1']['E'])
for CCD in [2,3]: BGD_lims['CCD%i' % CCD] = copy.deepcopy(BGD_lims['CCD1'])

class PSF0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs
    manifesto.update(OrderedDict(sorted([
            ('exptimes',([list],'Exposure times for each fluence.')),
            ('frames',([list],'Number of Frames for each fluence.')),
            ('wavelength',([int],'Wavelength')),
            ])))


class PSF0X(PointTask):
    
    inputsclass = PSF0X_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(PSF0X,self).__init__(inputs,log,drill,debug)
        self.name = 'PSF0X'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                         ('basic',self.basic_analysis),('bayes',self.bayes_analysis),
                         ('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = PSF0Xaux.PSF0Xfigs
        self.inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')
        
    def set_inpdefaults(self,**kwargs):
        
        testkey = kwargs['test']
        
        if 'PSF01' in testkey: 
            wavelength = kwargs['wavelength']
        elif 'PSF02' in testkey: 
            wavelength = 800

        exptimes = testdefaults['exptimes']['nm%i' % wavelength]
        frames = testdefaults['frames']

        self.inpdefaults = dict(wavelength=wavelength,
                                frames=frames,
                                exptimes=exptimes)
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)        
        
        self.perfdefaults['BGD_lims'] = BGD_lims # ADUs
        self.perfdefaults['Flu_lims']  = Flu_lims # ADUs
        self.perfdefaults['FWHM_lims'] = FWHM_lims # Pixels
    
    def build_scriptdict(self,diffvalues=dict(),elvis=pilib.elvis):
        """ 
        
        Builds PSF0X script structure dictionary.
        
        #:param exptimes: list of ints, [ms], exposure times.
        #:param frames: list of frame numbers. Same length as exptimes.
        #:param wavelength: int, [nm], wavelength.
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.
        
        """
        exptimes = self.inputs['exptimes']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        test = self.inputs['test']
        
        assert len(exptimes) == len(frames)
        
        FW_ID = ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])
        
        PSF0X_commvalues['wave'] = FW_IDX
        PSF0X_commvalues['mirr_pos'] = polib.mirror_nom['F%i' % FW_IDX]
        
        ncols = len(exptimes)
        
        PSF0X_sdict = dict()
        
        for ic in range(ncols):
            colid = 'col%i' % (ic+1,)
            PSF0X_sdict[colid]=dict(frames=frames[ic],exptime=exptimes[ic],
                       test=test)
    
        Ncols = len(PSF0X_sdict.keys())    
        PSF0X_sdict['Ncols'] = Ncols
        
        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(PSF0X_commvalues)               
        
        PSF0X_sdict = sc.update_structdict(PSF0X_sdict,commvalues,diffvalues)
        
        return PSF0X_sdict
    
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims,elvis=pilib.elvis):
        """
        """
        wavedkeys = []
        return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
                              wavedkeys=wavedkeys,elvis=elvis)    
            
    
    #def prep_data(dd,report,inputs,log=None):
    def prep_data(self):
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
        
        
        
        if self.report is not None: self.report.add_Section(keyword='prep_data',Title='Data Pre-Processing',level=0)
        
        # Inputs un-packing
        
        doMask = False
        doBias = False
        doOffset = True
        doFlats = False
        
        if 'mask' in self.inputs['inCDPs']:
            Maskdata = calibration.load_CDPs(self.inputs['inCDPs']['Mask'],ccd.CCD)
            doMask = True
            if self.log is not None:
                masksstr = self.inputs['inCDPs']['Mask'].__str__()
                masksstr = st.replace(masksstr,',',',\n')
                self.log.info('Applying cosmetics mask')
                self.log.info(masksstr)
        
        if 'bias' in self.inputs['inCDPs']:
            Biasdata = calibration.load_CDPs(self.inputs['inCDPs']['Bias'],ccd.CCD)
            doBias = True
            if self.log is not None:
                biasstr = self.inputs['inCDPs']['Bias'].__str__()
                biasstr = st.replace(biasstr,',',',\n')
                self.log.info('Subtracting Bias')
                self.log.info(biasstr)
        
        
        if 'FF' in self.inputs['inCDPs']:
            FFdata = calibration.load_CDPs(self.inputs['inCDPs']['FF'],FFing.FlatField)
            doFlats = True
            if self.log is not None:
                FFstr = self.inputs['inCDPs']['FF'].__str__()
                FFstr = st.replace(FFstr,',',',\n')
                self.log.info('Dividing by Flat-Fields')
                self.log.info(FFstr)
        
        # INDEX HANDLING
        
        dIndices = copy.deepcopy(self.dd.indices)
        
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
        
        self.dd.initColumn('ccdobj_name',CIndices,dtype='S100',valini='None')
        self.dd.initColumn('spots_name',CIndices,dtype='S100',valini='None')
        
        if not self.drill:
            
            rpath = self.inputs['resultspath']
            picklespath = self.inputs['subpaths']['pickles']
            spotspath = self.inputs['subpaths']['spots']
            
            for iObs in range(nObs):
                
                for jCCD,CCD in enumerate(CCDs):
                    
                    CCDkey = 'CCD%i' % CCD
                    
                    self.dd.mx['ccdobj_name'][iObs,jCCD] = '%s_proc' % self.dd.mx['File_name'][iObs,jCCD]
                    
                    dpath = self.dd.mx['datapath'][iObs,jCCD]
                    infits = os.path.join(dpath,'%s.fits' % self.dd.mx['File_name'][iObs,jCCD])
                    
                    ccdobj = ccd.CCD(infits)
                    
                    fullccdobj_name = os.path.join(picklespath,'%s.pick' % self.dd.mx['ccdobj_name'][iObs,jCCD]) 
                    
                    if doMask:
                        ccdobj.get_mask(Maskdata[CCDkey].extensions[-1])
                    
                    if doOffset:
                        for Quad in Quads:
                            ccdobj.sub_offset(Quad,method='median',scan='pre',trimscan=[5,5],
                                              ignore_pover=False)
                    if doBias:
                        ccdobj.sub_bias(Biasdata[CCDkey].extensions[-1],extension=-1)
                    
                    if doFlats:
                        ccdobj.divide_by_flatfield(FFdata[CCDkey].extensions[1].data,extension=-1)
                    
                    ccdobj.writeto(fullccdobj_name,clobber=True)
                    
                    
                    # Cut-out "spots"
                    
                    spots_array = np.zeros((nQuads,nSpots),dtype=object)
                    
                    for kQ,Quad in enumerate(Quads):
                        
                        for lS,SpotName in enumerate(SpotNames):
                        
                            coo = polib.Point_CooNom[CCDkey][Quad][SpotName]
                            lSpot = polib.extract_spot(ccdobj,coo,Quad,stampw=stampw)
                            
                            spots_array[kQ,lS] = copy.deepcopy(lSpot)
                    
                    
                    # save "spots" to a separate file and keep name in dd
                    
                    self.dd.mx['spots_name'][iObs,jCCD] = '%s_spots' % self.dd.mx['File_name'][iObs,jCCD]
                    
                    fullspots_name = os.path.join(spotspath,'%s.pick' % self.dd.mx['spots_name'][iObs,jCCD])
                    
                    spotsdict = dict(spots=spots_array)
                    
                    files.cPickleDumpDictionary(spotsdict,fullspots_name)
        
        if self.log is not None: self.log.info('Saved spot "bag" files to %s' % rpath)            
    
    
    #def basic_analysis(dd,report,inputs,log=None):
    def basic_analysis(self):
        """Performs basic analysis on spots:
             - shape from moments
             - Gaussian fitting: peak intensity, position, width_x, width_y
             
        """
        
        raise NotImplementedError
        
    
        
    #def bayes_analysis(dd,report,inputs,log=None):
    def bayes_analysis(self):
        """ 
        Performs bayesian decomposition of the spot images:
            - optomechanic PSF and detector PSF.
        Also measures R2, fwhm and ellipticity of "extracted" detector PSF.
        """
        raise NotImplementedError
    
    #def meta_analysis(dd,report,inputs,log=None):
    def meta_analysis(self):
        """
        
        Analyzes the relation between detector PSF and fluence.
        
        """
        raise NotImplementedError
        

