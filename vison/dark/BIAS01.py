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

from vison.pipe import lib as pilib
from vison.datamodel import  scriptic as sc
from vison.datamodel import ccd
from vison.image import calibration
from vison.datamodel import core
import B01aux
from vison.pipe.task import Task
from vison.image import performance

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
    
    
    def __init__(self,inputs,log=None,drill=False):
        """ """
        super(BIAS01,self).__init__(inputs,log,drill)
        self.name = 'BIAS01'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = B01aux.B01figs
        self.inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')
    
    def set_inpdefaults(self,**kwargs):
        self.inpdefaults = dict(N=25)
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        
    
    def build_scriptdict(self,diffvalues=dict(),elvis='6.3.0'):
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
        
        #raise RunTimeError # TEST        
    
        if self.report is not None: 
            self.report.add_Section(keyword='check_data',Title='Data Validation',level=0)
        
        
        # CHECK AND CROSS-CHECK HK
        
        if self.report is not None: 
            self.report.add_Section(keyword='check_HK',Title='HK',level=1)
        
        report_HK_perf = self.check_HK(HKKeys,reference='command',limits='P',tag='Performance',
                      doReport=self.report is not None,
                                 doLog=self.log is not None)
                
        HK_perf_ok = np.all([value for key,value in report_HK_perf.iteritems()])
        
        report_HK_safe = self.check_HK(HKKeys,reference='abs',limits='S',tag='Safe',
                      doReport = self.report is not None,
                          doLog = self.log is not None)
        
        HK_safe_ok = np.all([value for ke,value in report_HK_safe.iteritems()])
        
        if (not HK_perf_ok) or (not HK_safe_ok): self.dd.flags.add('HK_OOL')
        
        # Initialize new columns
    
        Xindices = copy.deepcopy(self.dd.indices)
        
        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad',vals=pilib.Quads))
        
        
        newcolnames_off = ['offset_pre','offset_img','offset_ove']
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off,Xindices,dtype='float32',valini=np.nan)
        
        newcolnames_std = ['std_pre','std_img','std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std,Xindices,dtype='float32',valini=np.nan)
        
        
        nObs,nCCD,nQuad = Xindices.shape
        CCDs = Xindices[1].vals
        Quads = Xindices[2].vals
        
        # Get statistics in different regions
        
        if not self.drill:
            
            for iObs in range(nObs):
                for jCCD in range(nCCD):
                    dpath = self.dd.mx['datapath'][iObs,jCCD]
                    ffits = os.path.join(dpath,'%s.fits' % \
                                         self.dd.mx['File_name'][iObs,jCCD])                    
                    ccdobj = ccd.CCD(ffits)
                    
                    for kQ in range(nQuad):
                        Quad = Quads[kQ]
                        
                        for reg in ['pre','img', 'ove']:
                            stats = ccdobj.get_stats(Quad,sector=reg,statkeys=['mean','std'],trimscan=[5,5],
                                    ignore_pover=True,extension=-1)
                            self.dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = stats[0]
                            self.dd.mx['std_%s' % reg][iObs,jCCD,kQ] = stats[1]
        
        #  METRICS ASSESSMENT
        
        # Assess metrics are within allocated boundaries
        
               
        if self.report is not None: 
            self.report.add_Section(keyword='check_ronoffset',Title='Offsets and RON',level=1)
        
        # absolute value of offsets
        
        offsets_lims = self.perflimits['offsets_lims']
        for reg in ['pre','img','ove']:
            arr = self.dd.mx['offset_%s' % reg]
            _compliance_offsets = self.check_stat_perCCD(arr,offsets_lims,CCDs)
            
            if not self.IsComplianceMatrixOK(_compliance_offsets): self.dd.flags.add('POORQUALDATA')
            if self.log is not None: self.addComplianceMatrix2Log(_compliance_offsets,label='COMPLIANCE OFFSETS [%s]:' % reg)        
            if self.report is not None: self.addComplianceMatrix2Report(_compliance_offsets,label='COMPLIANCE OFFSETS [%s]:' % reg)

        # cross-check of offsets: referred to pre-scan

        offsets_gradients = self.perflimits['offsets_gradients']
        for ireg,reg in enumerate(['img','ove']):            
            _lims = dict()
            for CCD in CCDs: _lims['CCD%i'%CCD] = offsets_gradients['CCD%i'%CCD][ireg+1]
            arr = self.dd.mx['offset_%s' % reg][:]-self.dd.mx['offset_pre'][:]
            _xcheck_offsets = self.check_stat_perCCD(arr,_lims,CCDs)
            
            if not self.IsComplianceMatrixOK(_xcheck_offsets): self.dd.flags.add('POORQUALDATA')
            if self.log is not None: self.addComplianceMatrix2Log(_xcheck_offsets,label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)        
            if self.report is not None: self.addComplianceMatrix2Report(_xcheck_offsets,label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)        
        
        # absolute value of std
        
        RONs_lims = self.perflimits['RONs_lims']
        for reg in ['pre','img','ove']:
            _compliance_std = self.check_stat_perCCD(self.dd.mx['std_%s' % reg],RONs_lims,CCDs)
        
            if not self.IsComplianceMatrixOK(_compliance_std): 
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('RON_OOL')
            if self.log is not None: self.addComplianceMatrix2Log(_compliance_std,label='COMPLIANCE RON [%s]:' % reg)        
            if self.report is not None: self.addComplianceMatrix2Report(_compliance_std,label='COMPLIANCE RON [%s]:' % reg)        
                
        #### PLOTS
        
        if self.report is not None: self.report.add_Section(keyword='check_plots',Title='Plots',level=1)
        
        # Plot offsets vs. time
                
        try:
            pmeta = dict(path = self.inputs['subpaths']['figs'],
                     stat='offset')
            self.doPlot('B01checks_offsets',**pmeta)
            self.addFigure2Report('B01checks_offsets')
        except:
            self.skipMissingPlot('BS_checkoffsets',ref='B01checks_offsets')

        # std vs. time
        
        try:
            pmeta = dict(path = self.inputs['subpaths']['figs'],
                     stat='std')
            self.doPlot('B01checks_stds',**pmeta)
            self.addFigure2Report('B01checks_stds')
        except:
            self.skipMissingPlot('BS_checkstds',ref='B01checks_stds')
            
        
        # Update Report, raise flags, fill-in
        
    
        if self.log is not None:
            self.addFlagsToLog()
        
        if self.report is not None:
            self.addFlagsToReport()
        
        
    
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
        
        if not self.drill:
            
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
        
    
    



