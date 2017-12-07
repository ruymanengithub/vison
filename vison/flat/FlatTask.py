#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:00:10 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF

from vison.pipe.task import Task
# END IMPORT

class FlatTask(Task):
    
    def __init__(self,*args,**kwargs):
        super(FlatTask,self).__init__(*args,**kwargs)
    
    
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
