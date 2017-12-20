#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:00:10 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
import os

from vison.pipe.task import Task
from vison.datamodel import core
from vison.datamodel import ccd
from vison.pipe import lib as pilib
# END IMPORT

class FlatTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(FlatTask,self).__init__(*args,**kwargs)
    
    
    def check_data(self):
        """ """
        test = self.inputs['test']
        if 'FLAT01'in test: # AD-HOC modification of test label
            kwargs = dict()
        elif 'FLAT02' in test:
            kwargs = dict()
        elif test == 'PTC01':
            kwargs = dict()
        elif 'PTC02' in test:
            kwargs = dict()
        elif test == 'NL01':
            kwargs = dict()
            
        Task.check_data(self,**kwargs)
    
    
    def get_checkstats_ST(self,**kwargs):
        """ """
        
        #test = self.inputs['test']
        
        # Initialize new columns
    
        Xindices = copy.deepcopy(self.dd.indices)
        
        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad',vals=pilib.Quads))
        
        
        newcolnames_off = ['offset_pre','offset_ove']
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off,Xindices,dtype='float32',valini=np.nan)
        
        self.dd.initColumn('flu_med_img',Xindices,dtype='float32',valini=np.nan)
        self.dd.initColumn('flu_var_img',Xindices,dtype='float32',valini=np.nan)
        
        newcolnames_std = ['std_pre','std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std,Xindices,dtype='float32',valini=np.nan)
        
        nObs,_,_ = Xindices.shape
        CCDs = Xindices[Xindices.names.index('CCD')].vals
        Quads = Xindices[Xindices.names.index('Quad')].vals
        
        # Get statistics in different regions
        
        
        if not self.drill:
            
            for iObs in range(nObs):
                for jCCD, CCD in enumerate(CCDs):
                    dpath = self.dd.mx['datapath'][iObs,jCCD]
                    ffits = os.path.join(dpath,'%s.fits' % \
                                         self.dd.mx['File_name'][iObs,jCCD])                    
                    ccdobj = ccd.CCD(ffits)
                    
                    for kQ,Quad in enumerate(Quads):
                        
                        for reg in ['pre','ove']:
                            stats_bias = ccdobj.get_stats(Quad,sector=reg,statkeys=['median','std'],trimscan=[5,5],
                                    ignore_pover=True,extension=-1)
                            self.dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = stats_bias[0]
                            self.dd.mx['std_%s' % reg][iObs,jCCD,kQ] = stats_bias[1]
                            
                        stats_img = ccdobj.get_stats(Quad,sector='img',statkeys=['median','std'],trimscan=[5,5],
                                ignore_pover=True,extension=-1)
                        self.dd.mx['flu_med_img'][iObs,jCCD,kQ] = stats_img[0]
                        self.dd.mx['flu_var_img'][iObs,jCCD,kQ] = stats_img[1]**2.

    
    def check_metrics_ST(self,**kwargs):
        """ 
        
        TODO:
            - offset levels (pre and over-scan), abs. and relative
            - RON in pre and overscan
            - fluence in image area [script-column-dependent]
            - variance in image area [script-column-dependent]
        
        """
        
        #test = self.inputs['test']
        
        Xindices = self.dd.indices
        CCDs = Xindices[Xindices.names.index('CCD')].vals
        
        if self.report is not None: 
            self.report.add_Section(keyword='check_ronoffset',Title='Offsets and RON',level=1)
        
        
        # absolute value of offsets
        
        offsets_lims = self.perflimits['offsets_lims']
        for reg in ['pre','ove']:
            arr = self.dd.mx['offset_%s' % reg]
            _compliance_offsets = self.check_stat_perCCD(arr,offsets_lims,CCDs)
            
            if not self.IsComplianceMatrixOK(_compliance_offsets): self.dd.flags.add('POORQUALDATA')
            if self.log is not None: self.addComplianceMatrix2Log(_compliance_offsets,label='COMPLIANCE OFFSETS [%s]:' % reg)        
            if self.report is not None: self.addComplianceMatrix2Report(_compliance_offsets,label='COMPLIANCE OFFSETS [%s]:' % reg)

        # cross-check of offsets: referred to pre-scan

        offsets_gradients = self.perflimits['offsets_gradients']
        for ireg,reg in enumerate(['ove']):            
            _lims = dict()
            for CCD in CCDs: _lims['CCD%i'%CCD] = offsets_gradients['CCD%i'%CCD][ireg+1]
            arr = self.dd.mx['offset_%s' % reg][:]-self.dd.mx['offset_pre'][:]
            _xcheck_offsets = self.check_stat_perCCD(arr,_lims,CCDs)
            
            if not self.IsComplianceMatrixOK(_xcheck_offsets): self.dd.flags.add('POORQUALDATA')
            if self.log is not None: self.addComplianceMatrix2Log(_xcheck_offsets,label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)        
            if self.report is not None: self.addComplianceMatrix2Report(_xcheck_offsets,label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)        
        
        # absolute value of std
        
        regs_std = ['pre','ove']
        RONs_lims = self.perflimits['RONs_lims']
        
        for reg in regs_std:
            _compliance_std = self.check_stat_perCCD(self.dd.mx['std_%s' % reg],RONs_lims,CCDs)
        
            if not self.IsComplianceMatrixOK(_compliance_std): 
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('RON_OOL')
            if self.log is not None: self.addComplianceMatrix2Log(_compliance_std,label='COMPLIANCE RON [%s]:' % reg)        
            if self.report is not None: self.addComplianceMatrix2Report(_compliance_std,label='COMPLIANCE RON [%s]:' % reg)
            
        
        # IMG FLUENCES        
        FLU_lims = self.perflimits['FLU_lims'] # dict
        
        _compliance_flu = self.check_stat_perCCDandCol(self.dd.mx['flu_med_img'],FLU_lims,CCDs)
        
        if not self.IsComplianceMatrixOK(_compliance_flu): 
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCE_OOL')
        if self.log is not None: self.addComplianceMatrix2Log(_compliance_flu,label='COMPLIANCE FLUENCE:')        
        if self.report is not None: self.addComplianceMatrix2Report(_compliance_flu,label='COMPLIANCE FLUENCE:')
            
