#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:56:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import  numpy as np
from pdb import  set_trace as stop
import copy
import os
from collections import OrderedDict

#from vison.support import context
from vison.inject import extract_injection_lines
from vison.inject import lib as ilib
from vison.pipe.task import Task
from vison.datamodel import core,ccd
#from vison.pipe import lib as pilib
from vison.support import context
#from vison.pipe.task import Task
# END IMPORT

lineoffsets = ilib.lineoffsets


class InjTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(InjTask,self).__init__(*args,**kwargs)
        
    
    def check_data(self,**kwargs):
        """ """
        
        test = self.inputs['test']
        if test == 'CHINJ01':
            kwargs = dict().update(kwargs)
        elif test == 'CHINJ02':
            kwargs = dict().update(kwargs)
        Task.check_data(self,**kwargs)
        
    def predict_expected_injlevels(self,teststruct):
        """ """        
        CCDs = [1,2,3]
        Quads = ['E','F','G','H']
        
        SecTags = dict(E='B',F='B',G='T',H='T')
        
        Ncols = teststruct['Ncols']
        
        expectation = OrderedDict()
        
        for CCD in CCDs:
            CCDkey = 'CCD%i' % CCD
            expectation[CCDkey] = OrderedDict()
            
            for Q in Quads:
                expectation[CCDkey][Q] = OrderedDict()                
                sectag = SecTags[Q]
                
                for icol in range(1,Ncols+1):
                    coldict = teststruct['col%i' % icol]
                    chinj = coldict['chinj']
                    IG2 = coldict['IG2_%s' % sectag]
                    IG1 = coldict['IG1_%i_%s' % (CCD,sectag)]
                    IDL = coldict['IDL']
                    IDH = coldict['IDH']
                    id_wid = coldict['id_wid']
                    id_dly = coldict['id_dly']
                    toi_ch = coldict['toi_ch']
                    
                    if chinj==0:
                        _inj = 0.
                    else:
                        _inj = ilib.predict_inj_level(IDL,IDH,IG1,IG2,id_wid,id_dly,toi_ch,sectag)
                    
                    expectation[CCDkey][Q]['col%i' % icol] = _inj
                    
        
        return expectation
    
    def get_FluenceAndGradient_limits(self):
        """ """
        
        tmpstructure = self.build_scriptdict(diffvalues={},elvis=self.elvis)
        inj_exp = self.predict_expected_injlevels(tmpstructure)
        
        Flu_lims = OrderedDict()
        FluGrad_lims = OrderedDict()
        
        for CCDkey in inj_exp.keys():
            Flu_lims[CCDkey] = OrderedDict()
            FluGrad_lims[CCDkey] = OrderedDict()
            
            for Q in inj_exp[CCDkey].keys():
                Flu_lims[CCDkey][Q] = OrderedDict()
                FluGrad_lims[CCDkey][Q] = OrderedDict()
                
                for colkey in inj_exp[CCDkey][Q].keys():
                    _inj = inj_exp[CCDkey][Q][colkey]
                    
                    if np.isnan(_inj):
                        Flu_lims[CCDkey][Q][colkey] = [-10.,1.01*2.**16]
                        FluGrad_lims[CCDkey][Q][colkey] = [10.,1.E4]
                    else:
                        Flu_lims[CCDkey][Q][colkey] = _inj * (1.+np.array([-0.5,0.5]))
                        FluGrad_lims[CCDkey][Q][colkey] = _inj * 0.3 * (1.+np.array([-0.5,0.5]))
                        
        return Flu_lims, FluGrad_lims

    def get_checkstats_ST(self,**kwargs):
        """ """
        
        #test = self.inputs['test']
        
        if 'pattern' in kwargs:
            pattern = kwargs['pattern']
        
        # Initialize new columns
    
        Xindices = copy.deepcopy(self.dd.indices)
        
        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad',vals=context.Quads))
        
        newcolnames_off = ['offset_pre','offset_ove']
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off,Xindices,dtype='float32',valini=np.nan)
        
        newcolnames_std = ['std_pre','std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std,Xindices,dtype='float32',valini=np.nan)
        
        for chk_inj_col in ['chk_mea_inject','chk_med_inject','chk_std_inject']:
            self.dd.initColumn(chk_inj_col,Xindices,dtype='float32',valini=np.nan)
        
        
        nObs,_,_ = Xindices.shape
        CCDs = Xindices[Xindices.names.index('CCD')].vals
        Quads = Xindices[Xindices.names.index('Quad')].vals
        
        # Get statistics in different regions
        
        if not self.drill:
            
            for iObs in range(nObs):
                
                if self.debug:
                    print 'check_data: processing Obsid %i/%i' % (iObs+1,nObs)
                
                for jCCD,CCD in enumerate(CCDs):
                    dpath = self.dd.mx['datapath'][iObs,jCCD]
                    ffits = os.path.join(dpath,'%s.fits' % \
                                         self.dd.mx['File_name'][iObs,jCCD])                    
                    ccdobj = ccd.CCD(ffits)
                    
                    vstart = self.dd.mx['vstart'][iObs][jCCD]
                    vend = self.dd.mx['vend'][iObs][jCCD]
                    dochinj = self.dd.mx['chinj'][iObs][jCCD]
                    if not 'pattern' in locals():
                        non = self.dd.mx['chinj_on'][iObs][jCCD]
                        noff = self.dd.mx['chinj_of'][iObs][jCCD]
                        nrep = (vend-vstart)/(non+noff)+1
                        pattern = (non,noff,nrep)
                    
                    for kQ,Quad in enumerate(Quads):
                        
                        
                        for reg in ['pre', 'ove']:
                            stats = ccdobj.get_stats(Quad,sector=reg,statkeys=['median','std'],trimscan=[5,5],
                                    ignore_pover=True,extension=-1)
                            self.dd.mx['offset_%s' % reg][iObs,jCCD,kQ] = stats[0]
                            self.dd.mx['std_%s' % reg][iObs,jCCD,kQ] = stats[1]
                                              
                        if dochinj:
                        
                            ccdobj.sub_offset(Quad,method='row',scan='pre',trimscan=[5,5],
                                              ignore_pover=True,extension=-1)
                            quaddata = ccdobj.get_quad(Quad,canonical=True,extension=-1)                            
                            extract_res = extract_injection_lines(quaddata,pattern,VSTART=vstart,
                                        VEND=vend,suboffmean=False,lineoffset=lineoffsets[Quad])
                            
                            self.dd.mx['chk_mea_inject'][iObs,jCCD,kQ] = extract_res['avinjection']
                            self.dd.mx['chk_med_inject'][iObs,jCCD,kQ] = extract_res['stats_injection'][0]
                            self.dd.mx['chk_std_inject'][iObs,jCCD,kQ] = extract_res['stats_injection'][1]
    
    def check_metrics_ST(self,**kwargs):
        """ 

        TODO:
            
            - offset levels (pre and over-scan), abs. and relative
            - RON in pre and overscan
            - mean fluence/signal in image area [script-column-dependent]
            - med fluence/signal in image area [script-column-dependent]
            - std in image area [script-column-dependent]
        
        
        """
        # test = self.inputs['test']
        
        Xindices = self.dd.indices
        CCDs = Xindices[Xindices.names.index('CCD')].vals
        
        if self.report is not None: 
            self.report.add_Section(keyword='check_metrics',Title='Offsets, RON, Injection levels',level=1)
                
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
            
        
        # IMG Signal Levels 
        Flu_lims = self.perflimits['Flu_lims'] # dict
        
        _compliance_flu = self.check_stat_perCCDQandCol(self.dd.mx['chk_med_inject'],Flu_lims,CCDs)
        
        
        if not self.IsComplianceMatrixOK(_compliance_flu): 
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCE_OOL')
        if self.log is not None: self.addComplianceMatrix2Log(_compliance_flu,label='COMPLIANCE FLUENCE:')        
        if self.report is not None: self.addComplianceMatrix2Report(_compliance_flu,label='COMPLIANCE FLUENCE:')
        
        # IMG std (injection-noise) levels
        
        FluGrad_lims = self.perflimits['FluGrad_lims'] # dict
        
        _compliance_flugrad = self.check_stat_perCCDQandCol(self.dd.mx['chk_std_inject'],FluGrad_lims,CCDs)
                
        if not self.IsComplianceMatrixOK(_compliance_flugrad): 
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCEGRAD_OOL')
        if self.log is not None: self.addComplianceMatrix2Log(_compliance_flugrad,label='COMPLIANCE FLUENCE GRADIENT:')        
        if self.report is not None: self.addComplianceMatrix2Report(_compliance_flugrad,label='COMPLIANCE FLUENCE GRADIENT:')
        