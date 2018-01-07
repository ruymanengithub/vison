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

from vison.inject import extract_injection_lines
from vison.lib import lineoffsets
from vison.pipe.task import Task
from vison.datamodel import core,ccd
from vison.pipe import lib as pilib
from vison.pipe.task import Task
# END IMPORT

class InjTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(InjTask,self).__init__(*args,**kwargs)
        
    
    def check_data(self):
        """ """
        test = self.inputs['test']
        if test == 'CHINJ01':
            kwargs = dict()
        elif test == 'CHINJ02':
            kwargs = dict()
        Task.check_data(self,**kwargs)
        
    def get_checkstats_ST(self,**kwargs):
        """ """
        
        #test = self.inputs['test']
        
        if 'pattern' in kwargs:
            pattern = kwargs['pattern']
        
        # Initialize new columns
    
        Xindices = copy.deepcopy(self.dd.indices)
        
        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad',vals=pilib.Quads))
        
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
        
        """
        
        test = self.inputs['test']
        
        stop()