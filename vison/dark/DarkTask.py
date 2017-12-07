#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:54:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
import os

from vison.datamodel import core,ccd
from vison.pipe import lib as pilib
from vison.pipe.task import Task
# END IMPORT

class DarkTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(DarkTask,self).__init__(*args,**kwargs)
        
    def get_checkstats_ST(self,**kwargs):
        """ """
        
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
        #CCDs = Xindices[1].vals
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
    
    def check_data(self):
        """ """
        if self.name == 'BIAS01':
            kwargs = dict(figkeys=['B01checks_offsets','B01checks_stds'])
        Task.check_data(self,**kwargs)

    def check_metrics_ST(self,*args,**kwargs):
        """ """
        
        Xindices = self.dd.indices
        CCDs = Xindices[Xindices.names.index('CCD')].vals
        
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
            