# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BF01

Brighter-Fatter Analysis
   Using data from test PTC01 or PTC02

Created on Wed Mar 7 10:57:00 2018

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

from vison.support import context
#from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import HKtools
from vison.datamodel import core, ccd
from vison.image import calibration
import ptc as ptclib
from vison.image import performance
#from FlatTask import FlatTask
from PTC0X import PTC0X
from vison.datamodel import inputs
#import BF01aux
from vison.analsysis import Guyonnet15 as G15

from vison.support.files import cPickleRead,cPickleDumpDictionary
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

#BF01_commvalues = dict(program='CALCAMP',
#  flushes=7,exptime=0.,vstart=0,vend=2086,
#  shuttr=1,
#  siflsh=1,siflsh_p=500,
#  wave = 4,
#  source='flat',
#  comments='')


#PTC01_relfluences = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])
#
#PTC01_exptimes = (PTC01_relfluences / 100.*ogse.tFWC_flat['nm800']).tolist() # ms
#PTC02waves = [590,640,730,800,880,0]
#
#PTC02_relfluences = np.array([10.,30.,50.,70.,80.,90.])
#
#PTC02TEMP_exptimes = (PTC02_relfluences / 100.*ogse.tFWC_flat['nm800']).tolist()
#
#testdefaults = dict(PTC01=dict(exptimes=PTC01_exptimes,
#                         frames=[10,10,10,10,10,10,10,10,4,4,4],
#                         wavelength=800),
#                    PTC02WAVE=dict(waves=PTC02waves,
#                                   frames=[4,4,4,4,4,4],
#                                   exptimes=dict()),
#                    PTC02TEMP=dict(frames=[4,4,4,4,4,4],
#                                   exptimes=PTC02TEMP_exptimes,
#                                   wavelength=800))
#
#for w in testdefaults['PTC02WAVE']['waves']:
#    testdefaults['PTC02WAVE']['exptimes']['nm%i' % w] = (PTC02_relfluences/100.*ogse.tFWC_flat['nm%i' % w]).tolist()
#
#
#plusminus10pcent = 1.+np.array([-0.10,0.10])
#
#FLU_lims_PTC01 = dict(CCD1= dict())
#for iflu,rflu in enumerate(PTC01_relfluences):
#    _cenval = min(rflu / 100.,1.) * 2.**16
#    _lims = _cenval * plusminus10pcent
#    FLU_lims_PTC01['CCD1']['col%i' % (iflu+1)] = _lims
#
#for i in [2,3]: FLU_lims_PTC01['CCD%i' % i] = copy.deepcopy(FLU_lims_PTC01['CCD1'])
#
#FLU_lims_PTC02 = dict(CCD1= dict())
#for iflu,rflu in enumerate(PTC02_relfluences):
#    _cenval = min(rflu / 100.,1.) * 2.**16
#    _lims = _cenval * plusminus10pcent
#    FLU_lims_PTC02['CCD1']['col%i' % (iflu+1)] = _lims
#
#for i in [2,3]: FLU_lims_PTC02['CCD%i' % i] = copy.deepcopy(FLU_lims_PTC02['CCD1'])


class BF01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('exptimes',([dict,list],'Exposure times for each fluence.')),
            ('frames',([list],'Number of Frames for each fluence.')),
            ('wavelength',([int],'Wavelength')),
            ])))


class BF01(PTC0X):
    """ """
    
    inputsclass = BF01_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(BF01,self).__init__(inputs,log,drill,debug)
        self.name = 'BF01'
        #self.type = 'Simple'
        self.subtasks = [('check',self.check_data),
                         ('prep',self.prepare_images),
                         ('extract',self.extract_PTC),
                    ('meta',self.meta_analysis)]
        #self.HKKeys = HKKeys
        self.figdict = dict() # BF01aux.gt_BF01figs(self.inputs['surrogate'])
        self.inputs['subpaths'] = dict(figs='figs',ccdpickles='ccdpickles')
        

    def set_inpdefaults(self,**kwargs):
        """ """
        
        # maskerading as PTC0X here... 
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']        
        super(BF01,self).set_inpdefaults(**_kwargs)        
        self.inpdefaults['test'] = kwargs['test']
        

    def set_perfdefaults(self,**kwargs):
        
        # maskerading as PTC0X here... 
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']        
        super(BF01,self).set_perfdefaults(**_kwargs)        
        
        

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """Builds PTC0X script structure dictionary.
        
        #:param exptimes: list of ints [ms], exposure times.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 800 nm.
        :param diffvalues: dict, opt, differential values.   
            
        """
        
        raise NotImplementedError("%s: This Task does not build a script, it uses data from another test" % self.name)
    
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """
        
        """    
        wavedkeys = ['motr_siz']
        return super(BF01,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
                              wavedkeys=wavedkeys,surrogate=self.inputs['surrogate'])
        
    def prepare_images(self):
        super(BF01,self).prepare_images(doExtract=True,doMask=True,
             doOffset=True,doBias=False,doFF=False)
    
    
    def extract_COV(self):
        """
        
        Performs basic analysis of images:
            - extracts COVARIANCE matrix for each fluence
        
             
        """
        
        if self.report is not None: self.report.add_Section(keyword='extractCOV',Title='Covariance-Matrices Extraction',level=0)    

        UP TO HERE        
                        
        label = self.dd.mx['label'][:,0].copy() # labels should be the same accross CCDs. PATCH.
        ObsIDs = self.dd.mx['ObsID'][:].copy()
        
        indices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuad = indices.shape
        
        Quads = indices[indices.names.index('Quad')].vals
        CCDs = indices[indices.names.index('CCD')].vals
        
#        emptyccdobj = ccd.CCD()    
#        tile_coos = dict()
#        for Quad in Quads:
#            tile_coos[Quad] = emptyccdobj.get_tile_coos(Quad,wpx,hpx)
#        Nsectors = tile_coos[Quads[0]]['Nsamps']    
#        sectornames = np.arange(Nsectors)
#        
#        Sindices = copy.deepcopy(self.dd.indices)
#        if 'Sector' not in Sindices.names:
#            Sindices.append(core.vIndex('Sector',vals=sectornames))
            
        # Initializing new columns
        
        valini = 0.
        
        self.dd.initColumn('sec_med',Sindices,dtype='float32',valini=valini)
        self.dd.initColumn('sec_var',Sindices,dtype='float32',valini=valini)
        
        # Pairing ObsIDs
        
        self.dd.initColumn('ObsID_pair',self.dd.mx['ObsID'].indices,dtype='int64',valini=0)
        
        ulabels = np.unique(label)
        
        for ulabel in ulabels:
            six = np.where(label == ulabel)
            nsix = len(six[0])
            ixeven = np.arange(0,nsix,2)
            ixodd = np.arange(1,nsix,2)
            
            self.dd.mx['ObsID_pair'][six[0][ixeven]] = ObsIDs[six[0][ixodd]]
            
        
        if not self.drill:
            
            #if self.proc_histo['Masked']:
            #    estimators = dict(median=np.ma.median,std=np.ma.std)
            #else:
            #    estimators = dict(median=np.median,std=np.std)
            
            dpath = self.inputs['subpaths']['ccdpickles']
            
            for iObs in range(nObs):
                
                _ObsID_pair = self.dd.mx['ObsID_pair'][iObs]
                if np.isnan(_ObsID_pair): continue
                iObs_pair = np.where(ObsIDs == _ObsID_pair)[0][0]
                
                for jCCD,CCD in enumerate(CCDs):
                                        
                    ccdobj_odd_f = os.path.join(dpath,self.dd.mx['ccdobj_name'][iObs,jCCD])
                    ccdobj_eve_f = os.path.join(dpath,self.dd.mx['ccdobj_name'][iObs_pair,jCCD])
                    
                    ccdobj_odd = copy.deepcopy(cPickleRead(ccdobj_odd_f)) # ['ccdobj'])
                    ccdobj_eve = copy.deepcopy(cPickleRead(ccdobj_eve_f)) # ['ccdobj'])
                    
                    evedata = ccdobj_eve.extensions[-1].data.copy()
                    
                    ccdobj_odd.sub_bias(evedata,extension=-1) # easy way to subtract one image from the other
                                        
                    for kQ in range(nQuad):
                        
                        Quad = Quads[kQ]
                        
                        _tile_coos = tile_coos[Quad]
                        
                        _meds = ccdobj_odd.get_tile_stats(Quad,_tile_coos,'median',extension=-1)
                        _vars = ccdobj_odd.get_tile_stats(Quad,_tile_coos,'std',extension=-1)**2.
                        
                        self.dd.mx['sec_med'][iObs,jCCD,kQ,:] = _meds.copy()
                        self.dd.mx['sec_var'][iObs,jCCD,kQ,:] = _vars.copy()
                        
        return None

    
    def extract_BF(self):
        """ """
    
    
    def meta_analysis(self):
        """
    
        Analyzes the BF results across fluences.
    
        
        """
        
        raise NotImplementedError
        
        if self.report is not None: self.report.add_Section(keyword='meta',Title='PTC Analysis',level=0)
        
        dIndices = copy.deepcopy(self.dd.indices)
        
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
                
                raw_var = self.dd.mx['sec_var'][:,iCCD,jQ,:]
                raw_med = self.dd.mx['sec_med'][:,iCCD,jQ,:]
                
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
        
    
        self.dd.products['gain'] = copy.deepcopy(gain_mx)
        self.dd.products['bloom'] = copy.deepcopy(bloom_mx)
        
        # Build Tables
        
        # Do plots
        
        # Add reports
        
        
