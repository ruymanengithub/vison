#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:55:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import copy
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.pipe.task import Task
from vison.point import lib as polib
from vison.datamodel import core, ccd
from vison.pipe import lib as pilib
# END IMPORT



class PointTask(Task):
    
    stampw = polib.stampw
    
    def __init__(self,*args,**kwargs):
        super(PointTask,self).__init__(*args,**kwargs)

    def check_data(self):
        """ """
        test = self.inputs['test']
        if 'PSF01' in test: 
            kwargs = dict()
        elif 'PSF02' in test:
            kwargs = dict()
        elif 'FOCUS00' in test:
            kwargs = dict()
        
        Task.check_data(self,**kwargs)
    
    
    def get_checkstats_ST(self,**kwargs):
        """ """
        
        # Initialize new columns
    
        Qindices = copy.deepcopy(self.dd.indices)
        
        if 'Quad' not in Qindices.names:
            Qindices.append(core.vIndex('Quad',vals=pilib.Quads))
        
        
        newcolnames_off = ['offset_pre','offset_ove']
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off,Qindices,dtype='float32',valini=np.nan)
        
        self.dd.initColumn('bgd_img',Qindices,dtype='float32',valini=np.nan)
        
        newcolnames_std = ['std_pre','std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std,Qindices,dtype='float32',valini=np.nan)
        
        Sindices = copy.deepcopy(Qindices)
        if 'Spot' not in Sindices.names:
            Sindices.append(core.vIndex('Spot',vals=polib.Point_CooNom['names']))
        
        
        self.dd.initColumn('chk_x',Sindices,dtype='float32',valini=np.nan)
        self.dd.initColumn('chk_y',Sindices,dtype='float32',valini=np.nan)
        self.dd.initColumn('chk_peak',Sindices,dtype='float32',valini=np.nan)
        self.dd.initColumn('chk_fluence',Sindices,dtype='float32',valini=np.nan)
        self.dd.initColumn('chk_fwhmx',Sindices,dtype='float32',valini=np.nan)
        self.dd.initColumn('chk_fwhmy',Sindices,dtype='float32',valini=np.nan)
         
        chkkeycorr = dict(chk_x='x',chk_y='y',chk_peak='peak',chk_fluence='fluence',
                          chk_fwhmx='fwhmx',chk_fwhmy='fwhmy')
        
        
        nObs,_,_ = Qindices.shape
        CCDs = Qindices[Qindices.names.index('CCD')].vals
        Quads = Qindices[Qindices.names.index('Quad')].vals
        Spots = Sindices[Sindices.names.index('Spot')].vals
        
        
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
                        
                        # To measure the background we mask out the sources
                         
                        alt_ccdobj = copy.deepcopy(ccdobj)
                         
                        mask_sources = polib.gen_point_mask(CCD,Quad,width=self.stampw,sources='all')
                         
                        alt_ccdobj.get_mask(mask_sources)
                        
                        alt_ccdobj.sub_offset(Quad,method='row',scan='pre',trimscan=[5,5],
                                              ignore_pover=True,extension=-1)
                         
                        imgstats = alt_ccdobj.get_stats(Quad,sector='img',statkeys=['median'],trimscan=[5,5],
                                    ignore_pover=True,extension=-1)
                         
                        self.dd.mx['bgd_img'][iObs,jCCD,kQ] = imgstats[0]
                         
                        alt_ccdobj = None
                         
                        for xSpot,SpotName in enumerate(Spots):
                                                        
                            coo = polib.Point_CooNom['CCD%i' % CCD][Quad][SpotName]                            
    
                            spot = polib.extract_spot(ccdobj,coo, Quad,log=self.log,
                                                      stampw=self.stampw)
                            
                            res_bas = spot.measure_basic(rin=10,rap=10,rout=-1)
                            
                            for chkkey in chkkeycorr:
                                self.dd.mx[chkkey][iObs,jCCD,kQ,xSpot] = res_bas[chkkeycorr[chkkey]]


    def check_stat_perCCDQSpot(self,arr,lims,CCDs=[1,2,3]):
        """ """
        #Qs = lims['CCD%i' % CCDs[0]].keys()
        Spots = polib.Point_CooNom['names']
        Qs = ccd.Quads
        
        compliance = OrderedDict()
        
        for iCCD,CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            compliance[CCDkey] = OrderedDict()
            
            for jQ,Q in enumerate(Qs):
                compliance[CCDkey][Q] = OrderedDict()                
                for kSpot,Spot in enumerate(Spots):
                    compliance[CCDkey][Q][Spot] = OrderedDict()
                    
                    _lims  = lims[CCDkey][Q][Spot]
                    
                    if isinstance(_lims,dict) or isinstance(_lims,OrderedDict):
                        colnames = _lims.keys()
                        compliance[CCDkey][Q][Spot] = OrderedDict()
                        for kcol,colname in enumerate(colnames):
                            compliance
                            _lims_col = lims[CCDkey][Q][Spot][colname]
                            ixsel = np.where(self.dd.mx['label'] == colname)
                            test = (np.isnan(arr[ixsel,iCCD,jQ,kSpot]) |\
                             (arr[ixsel,iCCD,jQ,kSpot] <= _lims_col[0]) | (arr[ixsel,iCCD,jQ,kSpot] >= _lims_col[1]))
                            compliance[CCDkey][Q][Spot][colname] = not np.any(test,axis=(0,1)).sum()
                    else:
                        
                        test = (np.isnan(arr[:,iCCD,jQ,kSpot]) |\
                             (arr[:,iCCD,jQ,kSpot] <= _lims[0]) | (arr[:,iCCD,jQ,kSpot] >= _lims[1]))
                        compliance[CCDkey][Q][Spot] = not np.any(test).sum()
                
        return compliance
    
    def check_metrics_ST(self,**kwargs):
        """         
        TO-CHECK:
            - offset levels (pre and over-scan), abs. and relative
            - RON in pre and overscan
            - background level in image area
            - spot fluences
            - spot sizes
        
        """
        
        #test = self.inputs['test']
        
        Xindices = self.dd.indices
        CCDs = Xindices[Xindices.names.index('CCD')].vals
        
        if self.report is not None: 
            self.report.add_Section(keyword='check_metrics',Title='Offsets, RON, Spots metrics',level=1)
        
        
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
            
        
        # Background Level
        
        BGD_lims = self.perflimits['BGD_lims'] # dict
        _compliance_bgd = self.check_stat_perCCDandQ(self.dd.mx['bgd_img'],BGD_lims,CCDs)
        
        if not self.IsComplianceMatrixOK(_compliance_bgd): 
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('BGD_OOL')
        if self.log is not None: self.addComplianceMatrix2Log(_compliance_bgd,label='COMPLIANCE BGD:')        
        if self.report is not None: self.addComplianceMatrix2Report(_compliance_bgd,label='COMPLIANCE BGD:')
        
        # Spot FWHM-(x**2+y**2_)**0.5
        
        FWHM_lims = self.perflimits['FWHM_lims'] # dict
        chk_fwhm = (self.dd.mx['chk_fwhmx'][:]**2.+self.dd.mx['chk_fwhmy'][:]**2.)**0.5
        _compliance_fwhm = self.check_stat_perCCDQSpot(chk_fwhm,FWHM_lims,CCDs)
        
        if not self.IsComplianceMatrixOK(_compliance_fwhm): 
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FOCUS_OOL')
        if self.log is not None: self.addComplianceMatrix2Log(_compliance_fwhm,label='COMPLIANCE FWHM(x2+y2)**0.5:')        
        if self.report is not None: self.addComplianceMatrix2Report(_compliance_fwhm,label='COMPLIANCE FWHM(x2+y2)**0.5:')        

        
        # Spot Fluence
        
        Flu_lims = self.perflimits['Flu_lims'] # dict
        _compliance_flu = self.check_stat_perCCDQSpot(self.dd.mx['chk_fluence'],Flu_lims,CCDs)

        if not self.IsComplianceMatrixOK(_compliance_flu): 
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCE_OOL')
        if self.log is not None: self.addComplianceMatrix2Log(_compliance_bgd,label='COMPLIANCE FLUENCE:')     
        if self.report is not None: self.addComplianceMatrix2Report(_compliance_bgd,label='COMPLIANCE FLUENCE:')
        
        
