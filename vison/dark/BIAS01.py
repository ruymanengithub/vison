#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BIAS01

Bias-structure/RON analysis script

Created on Tue Aug 29 16:53:40 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
import copy
import string as st
from collections import OrderedDict

from vison.support import context
from vison.pipe import lib as pilib
from vison.datamodel import  scriptic as sc
from vison.datamodel import ccd
from vison.image import calibration
from vison.datamodel import core
import B01aux
#from vison.pipe.task import Task
from DarkTask import DarkTask
from vison.image import performance
from vison.datamodel import inputs
#from vison.support.report import Report


# END IMPORT

isthere = os.path.exists


HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']


BIAS01_commvalues = dict(program='CALCAMP',test='BIAS01',
  flushes=7,exptime=0.,shuttr=0,
  e_shuttr=0,vstart=0,vend=2086,
  siflush=0,#sinvflushp=500,
  chinj=0,
  s_tpump=0,
  v_tpump=0,
  motr_on=0,
  toi_fl=143.,toi_tp=1000.,toi_ro=1000.,toi_ch=1000.,
  wave=4,
  comments='BIAS')

class BIAS01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('N',([int],'Number of Frame Acquisitions.')),
            ])))

class BIAS01(DarkTask):
    """ """
    
    inputsclass = BIAS01_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(BIAS01,self).__init__(inputs,log,drill,debug)
        self.name = 'BIAS01'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = B01aux.B01figs
        self.inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')
        
        
    def set_inpdefaults(self,**kwargs):
        self.inpdefaults = self.inputsclass(N=25)
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        
    
    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """Builds BIAS01 script structure dictionary.
        
        ###:param N: integer, number of frames to acquire.
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.
        
        """
        
        N = self.inputs['N']
        BIAS01_sdict = dict(col1=dict(frames=N,exptime=0))
    
        Ncols = len(BIAS01_sdict.keys())
        BIAS01_sdict['Ncols'] = Ncols
                    
        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(BIAS01_commvalues)
        
        if len(diffvalues)==0:
            diffvalues = self.inputs['diffvalues']
        
        BIAS01_sdict = sc.update_structdict(BIAS01_sdict,commvalues,diffvalues)
    
        return BIAS01_sdict
    
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = []
        return super(BIAS01,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys)
    
    
    def prep_data(self):
        """
        
        BIAS01: Preparation of data for further analysis.
        applies a mask
        
        **METACODE**
        
        ::
            f.e. ObsID:
                f.e.CCD:                
                    apply cosmetic mask, if available                
                    f.e.Q:    
                        subtract offset
                    save file as a datamodel.ccd.CCD object.
        
        """
        super(self).prepare_images(doExtract=True,doMask=True,doOffset=True)
        
#==============================================================================
#         if self.report is not None: 
#             self.report.add_Section(keyword='prep_data',Title='Data Pre-Processing',level=0)
#         
#                 
#         doMask = True
#         doOffset = True
#         
#         if 'mask' in self.inputs['inCDPs']:
#             Maskdata = calibration.load_CDPs(self.inputs['inCDPs']['Mask'],ccd.CCD)
#             if self.log is not None:
#                 masksstr = self.inputs['inCDPs']['Mask'].__str__()
#                 masksstr = st.replace(masksstr,',',',\n')
#                 self.log.info('Applying cosmetics mask')
#                 self.log.info(masksstr)    
#         
#         # Initialize new columns
#     
#         Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
#         
#         
#         self.dd.initColumn('ccdobj_name',Cindices,dtype='S100',valini='None')
#         
#         DDindices = copy.deepcopy(self.dd.indices)
#         
#         nObs,nCCD,nQuad = DDindices.shape
#         Quads = DDindices[2].vals
#         CCDs = DDindices[DDindices.names.index('CCD')].vals
#         
#         if not self.drill:
#             
#             picklespath = self.inputs['subpaths']['pickles']
#             
#             
#             for iObs in range(nObs):
#                 
#                 for jCCD,CCD in enumerate(CCDs):
#                     
#                     CCDkey = 'CCD%i' % CCD
#                     
#                     ccdobj_name = '%s_proc' % self.dd.mx['File_name'][iObs,jCCD]
#                     
#                     dpath = self.dd.mx['datapath'][iObs,jCCD]
#                     infits = os.path.join(dpath,'%s.fits' % self.dd.mx['File_name'][iObs,jCCD])
#                     
#                     ccdobj = ccd.CCD(infits)
#                     
#                     fullccdobj_name = os.path.join(picklespath,'%s.pick' % self.dd.mx['ccdobj_name'][iObs,jCCD]) 
#                     
#                     if doMask:
#                         ccdobj.get_mask(Maskdata[CCDkey].extensions[-1])
#                     
#                     if doOffset:
#                         
#                         for Quad in Quads:
#                             ccdobj.sub_offset(Quad,method='median',scan='pre',trimscan=[5,5],
#                                               ignore_pover=False)
#                     
#                     ccdobj.writeto(fullccdobj_name,clobber=True)
#                     
#                     self.dd.mx['ccdobj_name'][iObs,jCCD] = ccdobj_name
#     
#         
#         return None
#         
#==============================================================================
    def basic_analysis(self):
        """ 
        
        BIAS01: Basic analysis of data.
        
        **METACODE**
        
        :: 
    
            f. e. ObsID:
               f.e.CCD:
                   
                   load ccdobj of ObsID, CCD
                   
                   with ccdobj, f.e.Q:
                       produce a 2D poly model of bias, save coefficients
                       produce average profile along rows
                       produce average profile along cols
                       save 2D model and profiles in a pick file for each OBSID-CCD
                       measure and save RON after subtracting large scale structure
            plot RON vs. time f. each CCD and Q
            plot average profiles f. each CCD and Q (color coded by time)
     
        
        """
        
        raise NotImplementedError
       
        
        
    def meta_analysis(self):
        """
        
        **METACODE**
        
        ::
        
            f. each CCD:
               f. e. Q:
                   stack all ObsIDs to produce Master Bias
                   measure average profile along rows
                   measure average profile along cols
            plot average profiles of Master Bias f. each Q
            produce table with summary of results, include in report
            show Master Bias (image), include in report
            save name of MasterBias to DataDict, report
            
        """
        
        raise NotImplementedError
