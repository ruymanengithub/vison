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

from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
from vison.image import calibration
from vison.datamodel import core
import B01aux
from vison.support import files
from vison.pipe.task import Task

#from vison.support.report import Report


# END IMPORT

isthere = os.path.exists


HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']


BIAS01_commvalues = dict(program='CALCAMP',test='BIAS01',
  flushes=7,exptime=0.,shuttr=1,
  e_shuttr=0,vstart=1,vend=2086,
  siflush=0,#sinvflushp=500,
  chinj=0,
  s_tpump=0,
  v_tpump=0,
  motr_on=0,
  toi_fl=143.,toi_tp=1000.,toi_ro=1000.,toi_ch=1000.,
  wave=4,
  comments='BIAS')
  
class BIAS01(Task):
    
    def __init__(self,inputs,log=None):
        """ """
        super(BIAS01,self).__init__(inputs,log)
        self.name = 'BIAS01'
        self.HKKeys = HKKeys
        self.figdict = B01aux.B01figs
        

    def build_scriptdict(self,N,diffvalues=dict(),elvis='6.3.0'):
        """Builds BIAS01 script structure dictionary.
        
        :param N: integer, number of frames to acquire.
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.
        
        """
        
        BIAS01_sdict = dict(col1=dict(frames=N,exptime=0))
    
        Ncols = len(BIAS01_sdict.keys())    
        BIAS01_sdict['Ncols'] = Ncols
    
                    
        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(BIAS01_commvalues)
        
        BIAS01_sdict = sc.update_structdict(BIAS01_sdict,commvalues,diffvalues)
    
        return BIAS01_sdict
    
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims,elvis='6.3.0'):
        """ """
        wavedkeys = []
        return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys,elvis=elvis)
    
    
    def check_data(self):
        """ 
        
        BIAS01: Checks quality of ingested data.
        
        **METACODE**
        
        
        **TODO**: consider to raise an exception that
              would halt execution of task if 
              processing data could be just a waste of time.
              
        **TODO**: consider add a common binary "flags" variable as 
              input/output. It could go in DataDict, and reported in 
              log and report.
        
        ::
          
          check common HK values are within safe / nominal margins
          check voltages in HK match commanded voltages, within margins
        
          f.e.ObsID:
              f.e.CCD: 
                  f.e.Q.:
                      measure offsets in pre-, img-, over-
                      measure std in pre-, img-, over-
          assess std in pre- is within allocated margins
          assess offsets in pre-, img-, over- are equal, within allocated  margins
          assess offsets are within allocated margins
        
          plot offsets vs. time
          plot std vs. time
        
          issue any warnings to log
          issue update to report
          update flags as needed
        
        """
    
        if self.report is not None: 
            self.report.add_Section(keyword='check_data',Title='Data Validation',level=0)
    
        bypass = True # TESTS
        
        # ALIASES
        dd = self.dd
        report = self.report
        inputs = self.inputs
        
        # CHECK AND CROSS-CHECK HK
        
        #print 'HK-perf' # TESTS
        report_HK_perf = HKtools.check_HK_vs_command(HKKeys,dd,limits='P',elvis=self.elvis)
        #print 'HK-safe' # TESTS
        report_HK_safe = HKtools.check_HK_abs(HKKeys,dd,limits='S',elvis=self.elvis)
        
    
        # Initialize new columns
    
        Xindices = copy.deepcopy(dd.indices)
        
        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad',vals=pilib.Quads))
        
        
        newcolnames_off = ['offset_pre','offset_img','offset_ove']
        for newcolname_off in newcolnames_off:
            dd.initColumn(newcolname_off,Xindices,dtype='float32',valini=np.nan)
        
        newcolnames_std = ['std_pre','std_img','std_ove']
        for newcolname_std in newcolnames_std:
            dd.initColumn(newcolname_std,Xindices,dtype='float32',valini=np.nan)
        
        
        nObs,nCCD,nQuad = Xindices.shape
        CCDs = Xindices[1].vals
        Quads = Xindices[2].vals
        
        # Get statistics in different regions
        
        if not bypass:
        
            for iObs in range(nObs):
                
                for jCCD in range(nCCD):
                    
                    dpath = dd.mx['datapath'][iObs,jCCD]
                    ffits = os.path.join(dpath,dd.mx['Files'][iObs,jCCD])
                    
                    ccdobj = ccd.CCD(ffits)
                    
                    for kQ in range(nQuad):
                        Quad = Quads[kQ]
                        
                        for reg in ['pre','img', 'ove']:
                            stats = ccdobj.get_stats(Quad,sector=reg,statkeys=['mean','std'],trimscan=[5,5],
                                    ignore_pover=True,extension=-1)
                            dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = stats[0]
                            dd.mx['std_%s' % reg][iObs,jCCD,kQ] = stats[1]
                    
        # Assess metrics are within allocated boundaries
    
        offset_lims = [1000.,3500.] # TESTS, should be in common limits file
        offset_diffs = dict(img=[-1.,+5.],ove=[-1.,+6.]) # TESTS, should be in common limits file
        
        std_lims = [0.5,2.] # TESTS, should be in common limits file
        
        stop()
        
        # absolute value of offsets
        
        compliance_offsets = dict()
        for reg in ['pre','img','ove']:
            
            test = ((dd.mx['offset_%s' % reg][:] <= offset_lims[0]) |\
                          (dd.mx['offset_%s' % reg][:] >= offset_lims[1]))
            compliance_offsets[reg] = np.any(test,axis=(1,2)).sum()
        
        
        # cross-check of offsets
            
        xcheck_offsets = dict()
        for reg in ['img','ove']:
            
            test = dd.mx['offset_%s' % reg][:]-dd.mx['offset_pre'][:]
            
            testBool = (test <= offset_diffs[reg][0]) | \
                   (test >= offset_diffs[reg][1])
            xcheck_offsets[reg] = np.any(testBool,axis=(1,2)).sum()
        
        # absolute value of std
            
        compliance_std = dict()
        for reg in ['pre']:
            
            test = ((dd.mx['std_%s' % reg][:] <= std_lims[0]) |\
                          (dd.mx['std_%s' % reg][:] >= std_lims[1]))
            compliance_std[reg] = np.any(test,axis=(1,2)).sum()
        
        stop()
        
        # Plot offsets vs. time
        
        for CCD in CCDs:
            pmeta = dict(CCD=CCD,path = inputs['subpaths']['figs'])
            try:
                self.doPlot('B01offsets_CCD%i' % CCD,**pmeta)
            except: pass
            self.addFigure2Report('B01offsets_CCD%i' % CCD)
        
        
        # std vs. time
        
        for CCD in CCDs:
            pmeta = dict(CCD=CCD,path = inputs['subpaths']['figs'])
            try:
                self.doPlot('B01stds_CCD%i' % CCD,**pmeta)
            except: pass
            self.addFigure2Report('B01stds_CCD%i' % CCD)
        
        
        # Update Report, raise flags, fill-in log
    
        if self.log is not None:
            self.log.info('Reporting and Flagging MISSING in BIAS01.check_data')    
        
        
        
    
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
        
        if self.report is not None: 
            self.report.add_Section(keyword='prep_data',Title='Data Pre-Processing',level=0)
        
        
        bypass = True # TESTS
        
        doMask = True
        doOffset = True
        
        if 'mask' in self.inputs['inCDPs']:
            Maskdata = calibration.load_CDPs(self.inputs['inCDPs']['Mask'],ccd.CCD)
            doMask = True
            if self.log is not None:
                masksstr = self.inputs['inCDPs']['Mask'].__str__()
                masksstr = st.replace(masksstr,',',',\n')
                self.log.info('Applying cosmetics mask')
                self.log.info(masksstr)    
        
        # Initialize new columns
    
        Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
        
        
        self.dd.initColumn('ccdobj_name',Cindices,dtype='S100',valini='None')
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs,nCCD,nQuad = DDindices.shape
        Quads = DDindices[2].vals
        CCDs = DDindices[DDindices.names.index('CCD')].vals
        
        if not bypass:
            
            picklespath = self.inputs['subpaths']['pickles']
            
            
            for iObs in range(nObs):
                
                for jCCD,CCD in enumerate(CCDs):
                    
                    CCDkey = 'CCD%i' % CCD
                    
                    ccdobj_name = '%s_proc' % self.dd.mx['File_name'][iObs,jCCD]
                    
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
                    
                    ccdobj.writeto(fullccdobj_name,clobber=True)
                    
                    self.dd.mx['ccdobj_name'][iObs,jCCD] = ccdobj_name
    
        
        return None
        
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
        
    
    
    def feeder(self,inputs,elvis='6.3.0'):
        """ """
        
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        
        N = inputs['N']
        if 'elvis' in inputs:
            self.elvis = inputs['elvis']
        if 'diffvalues' in inputs:
            diffvalues = inputs['diffvalues']
        else:
            diffvalues = {}
        
        
        scriptdict = self.build_scriptdict(N,diffvalues,elvis=elvis)
        
        inputs['structure'] = scriptdict        
        inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')
        
        return inputs




