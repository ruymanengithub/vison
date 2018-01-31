#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Generic Task (Test) Class.

Created on Tue Nov 14 14:20:04 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import datetime
import string as st
import copy
from collections import OrderedDict
import sys,traceback

from vison.support.report import Report
from vison.support import vistime
from vison.support import files
import lib as pilib
from vison.support import context
import task_lib as tlib
from vison.datamodel import ccd
from vison.image import calibration
from vison.ogse import ogse
from vison import __version__
# END IMPORT

isthere = os.path.exists

class Task(object):
    """ """
    
    from task_lib import check_HK,filterexposures
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        self.ID = None
        self.internals = dict()
        #self.inputs = dict()
        self.inputs = self.inputsclass()
        self.inpdefaults = dict()
        self.perfdefaults = dict()
        #self.set_defaults()
        self.elvis = context.elvis
        self.log = log
        self.name = ''
        self.type = 'Simple'
        self.HKKeys = []
        self.figdict = dict()
        self.subtasks = [()]
        self.perflimits = dict()
        self.drill = drill
        self.debug = debug
        
        self.set_inpdefaults(**inputs)
        _inputs = self.inpdefaults        
        _inputs.update(inputs)
        self.inputs.update(_inputs)
        
        if 'elvis' in self.inputs:
             self.elvis = self.inputs['elvis']
        
        self.set_perfdefaults(**inputs)
        _perfdefaults = self.perfdefaults
        self.perflimits.update(_perfdefaults)
        if 'perflimits' in self.inputs:
            self.perflimits.update(self.inputs['perflimits']) 
        
        if 'diffvalues' in inputs: diffvalues = inputs['diffvalues']
        else: diffvalues = {}
        
        self.inputs['structure'] = self.build_scriptdict(diffvalues,elvis=self.elvis)
    
    def set_inpdefaults(self,**kwargs):
        pass
    
    def set_perfdefaults(self,**kwargs):
        pass
    
    def build_scriptdict(self,diffvalues={},elvis=context.elvis):
        """ """
        return dict()
    
    def __call__(self):
        """Generic test master function."""
        
        # INPUTS
        
        subtasks = self.subtasks
        
        
        # inputs loading
        resultspath = self.inputs['resultspath']
        try: _paths = self.inputs['subpaths']
        except: _paths = dict()
        testkey = self.inputs['test']
        todo_flags = self.inputs['todo_flags']
        try: reportroot = self.inputs['reportroot']
        except KeyError: reportroot = '%s_report' % testkey
        try: cleanafter = self.inputs['cleanafter']
        except KeyError: cleanafter = False

        DataDictFile = os.path.join(resultspath,'%s_DataDict.pick' % testkey)
        reportobjFile = os.path.join(resultspath,'%s_Report.pick' % testkey)
        
        if todo_flags['init']:
            
            
            # Creating resultspath
            if not isthere(resultspath):
                os.system('mkdir %s' % resultspath)
            
            # Creating subresultspath
            for _pathkey in _paths:
                subpath = os.path.join(resultspath,_paths[_pathkey])            
                if not isthere(subpath): os.system('mkdir %s' % subpath)
                _paths[_pathkey] = subpath
            
            self.inputs['subpaths'] = _paths

            # Let's start from scratch
            
            if os.path.exists(DataDictFile): os.system('rm %s' % DataDictFile)
            if os.path.exists(reportobjFile): os.system('rm %s' % reportobjFile)
            
            os.system('rm %s/*' % self.inputs['resultspath'])
            if 'subpaths' in self.inputs:
                for _pathkey in self.inputs['subpaths'].keys():
                    subpath = self.inputs['subpaths'][_pathkey]
                    os.system('rm %s/*' % subpath)
            
            # Initialising Report Object
        
            if todo_flags['report']:
                self.report = Report(TestName=self.name)
                self.report.add_Section(keyword='init',Title='Inputs \& Data Ingestion',level=0)
                self.add_inputs_to_report()                
            else:
                self.report = None
            
            
            if self.type == 'Simple':
                self.ingest_data_SimpleTest()           
            elif self.type == 'Meta':
                self.ingest_data_MetaTest()
            
            #self.save_progress(dd,reportobj,DataDictFile,reportobjFile)
            self.save_progress(DataDictFile,reportobjFile)
            
        else:
            #dd, reportobj = pilib.recover_progress(DataDictFile,reportobjFile)
            self.recover_progres(DataDictFile,reportobjFile)
        
        
        # DATA-WORK and ANALYSIS        
        
        for subtask in subtasks:
            
            subtaskname, subtaskmethod = subtask
            
            if self.log is not None:
                self.log.info('%s: %s' % (subtaskname,subtaskmethod.__module__))
            
            if subtaskname not in todo_flags:
                todo_flags[subtaskname] = False
            
            if todo_flags[subtaskname]:
                
                tini = datetime.datetime.now()
                try: 
                    subtaskmethod()
                    tend = datetime.datetime.now()
                    dtm = ((tend-tini).seconds)/60.
                    if self.log is not None: 
                        self.log.info('%.1f minutes in running Sub-task: %s' % (dtm,subtaskname))
                    self.save_progress(DataDictFile,reportobjFile)
                except:
                    self.catchtraceback()
                    self.save_progress(DataDictFile,reportobjFile)
                    if not self.debug:
                        if self.log is not None:
                            self.log.info('SUBTASK "%s:%s" FAILED, QUITTING!' % (subtaskname,subtaskmethod.__name__))
                        break
                    else:
                        sys.exit()
                
            else:
                self.recover_progress(DataDictFile,reportobjFile)
        

        # Write automatic Report of Results

        if todo_flags['report']:
            self.report.doreport(reportroot,cleanafter)
            outfiles = self.report.writeto(reportroot,cleanafter)
            
            for outfile in outfiles:
                os.system('mv %s %s/' % (outfile,resultspath))
        

        self.save_progress(DataDictFile,reportobjFile)
        
        if self.log is not None:
            self.log.info('Finished %s' % self.name)


    def catchtraceback(self):
        """ """
        msg_trbk = traceback.format_exc()
        if self.log is not None and not self.debug:
            self.log.info(msg_trbk)
        else:
            print msg_trbk
    
    
    def addHK_2_dd(self):
        """ """
        self.dd = pilib.addHK(self.dd,self.HKKeys,elvis=self.elvis)

    def ingest_data_SimpleTest(self):
        
        testkey = self.inputs['test']
        datapath = self.inputs['datapath']
        OBSID_lims = self.inputs['OBSID_lims']
        structure = self.inputs['structure']
        explogf = self.inputs['explogf']
        #elvis = self.inputs['elvis']
        
        # META-DATA WORK
        
        explog, checkreport = self.filterexposures(structure,explogf,datapath,OBSID_lims)
        
        if self.log is not None:
            self.log.info('%s acquisition consistent with expectations: %s' % (testkey,checkreport['checksout']))
            if len(checkreport['failedcols'])>0:
                self.log.info('%s failed columns: %s' % (testkey,checkreport['failedcols']))
            if len(checkreport['failedkeys'])>0:
                self.log.info('%s failed keys: %s' % (testkey,checkreport['failedkeys']))
        
        if self.report is not None:
            ntestkey = st.replace(testkey,'_','\_')
            nchecksout = ['\\bf{%s}' % checkreport['checksout']]
            nchecksout = [st.replace(item,'False','$\\textcolor{red}{\\bf{False}}$') for item in nchecksout][0]
            self.report.add_Text('%s acquisition consistent with expectations: %s\\newline' % (ntestkey,nchecksout))
            
            if (checkreport['failedcols'])>0:          
                nfailedcols = st.replace(checkreport['failedcols'].__repr__(),'_','\_')
                self.report.add_Text('%s failed columns: %s' % (ntestkey,nfailedcols))
            if len(checkreport['failedkeys'])>0:
                nfailedkeys = st.replace(checkreport['failedkeys'].__repr__(),'_','\_')
                self.report.add_Text('%s failed keys: %s' % (ntestkey,nfailedkeys))
        
        # Adding Time Axis            
        
        explog['time'] = np.array(map(vistime.get_dtobj,explog['date'])).copy()
        
        # Building DataDict 
        
        self.dd = pilib.DataDict_builder(explog,self.inputs,structure)
        
        if not checkreport['checksout']:self.dd.flags.add('MISSDATA')
        
        # Add HK information
        self.addHK_2_dd()


    def ingest_data_MetaTest(self):
        raise NotImplementedError("Method implemented in child-class")
    
    
    def save_progress(self,DataDictFile,reportobjFile):
        """ """
        files.cPickleDumpDictionary(self.dd,DataDictFile)
        files.cPickleDump(self.report,reportobjFile)
        
        
    def recover_progress(self,DataDictFile,reportobjFile):
        """ """
        self.dd = files.cPickleRead(DataDictFile)
        self.report = files.cPickleRead(reportobjFile)
    
    
    def addFigure2Report(self,figkey):
        """ """
        figobj = self.figdict[figkey][0]
        figname = figobj.figname
        texfraction = figobj.texfraction
        caption = figobj.caption
        assert os.path.exists(figname)
        epsname = '%s.eps' % os.path.splitext(figname)[0]
        os.system('convert %s %s' % (figname,epsname))
        self.report.add_Figure(epsname,texfraction=texfraction,
                               caption=caption,label=figkey)
        
    def doPlot(self,figkey,**kwargs):
        """ """
        figobj = self.figdict[figkey][0]()
        figobj.configure(**kwargs)
        figobj.build_data(self)
        if 'meta' in kwargs: meta = kwargs['meta']
        else: meta = {}       
        figobj.plot(**meta)
        self.figdict[figkey][0] = copy.deepcopy(figobj)
        
    def addComplianceMatrix2Log(self,complidict,label=''):
        """ """
        st_compl = complidict.__str__()
        self.log.info('%s\n%s' % (label,st_compl))
        
        
    def addComplianceMatrix2Report(self,complidict,label=''):
        """ """
        nicelabel = st.replace(label,' ','\ ')
        #st_compl = complidict.__str__()
        st_compl = tlib.convert_compl_to_nesteditemlist(complidict)
        nice_st_compl = [st.replace(item,'False','$\\textcolor{red}{\\bf{False}}$') for item in st_compl]
        msgList = ['$\\bf{%s}$' % nicelabel] +\
                  ['\\begingroup'] +\
                  ['\\scriptsize'] +\
                  nice_st_compl +\
                  ['\endgroup']
#                   ['\\\\']+\
                   
        self.report.add_Text(msgList)
        
    def IsComplianceMatrixOK(self,complidict):
        """ """

        def traverse_tree(dictionary,isOK):
        
            for key,value in dictionary.items():
                #print 'Upper: %s' % key
                if isinstance(value,(dict,OrderedDict)):
                    #print key,value
                    isOK = isOK and traverse_tree(value,isOK)
                else:
                    #print key,value
                    isOK = isOK and value
            return isOK
        
        isOK = traverse_tree(complidict,True)
                
        return isOK
        
        
    def addFlagsToLog(self):
        """ """
        flagstxt = st.join(self.dd.flags.getFlagsOnList(),', ')
        self.log.info('FLAGS ON}:\n%s' % flagstxt)
    
    def addFlagsToReport(self):
        """ """
        niceflagnames = [st.replace(item,'_','\_') for item in self.dd.flags.getFlagsOnList()]
        flagstxt = st.join(niceflagnames,', ')
        msgList = ['$\\bf{FLAGS\ ON}$: ',flagstxt]
        self.report.add_Text(msgList)
        
    def skipMissingPlot(self,key,ref):
        """ """
        self.figdict[key] = copy.deepcopy(self.figdict['BlueScreen'])
        niceref = st.replace(ref,'_','\_')
        figspath = self.inputs['subpaths']['figs']
        pmeta = dict(path=figspath,
                caption = '$\\bf{MISSING}:$ %s' % niceref,
                tag=ref,
                     meta=dict(title=niceref))
        self.doPlot(key,**pmeta)
        self.addFigure2Report(key)
    
    
    def check_stat_perCCD(self,arr,CCDlims,CCDs=[1,2,3]):
        """ """
        compliance = OrderedDict()
        for iCCD,CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            #compliance[CCDkey] = OrderedDict()
            _lims = CCDlims[CCDkey]
            test = (np.isnan(arr[:,iCCD,...]) |\
                    (arr[:,iCCD,...] <= _lims[0]) | (arr[:,iCCD,...] >= _lims[1]))
            compliance[CCDkey] = not np.any(test,axis=(0,1)).sum()
        return compliance
    
    def check_stat_perCCDandQ(self,arr,CCDQlims,CCDs=[1,2,3]):
        """ """
        compliance = OrderedDict()
        for iCCD,CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            compliance[CCDkey] = OrderedDict()
            for jQ,Q in enumerate(ccd.Quads):
                _lims = CCDQlims[CCDkey][Q]
                test = (np.isnan(arr[:,iCCD,jQ,...]) |\
                    (arr[:,iCCD,jQ,...] <= _lims[0]) | (arr[:,iCCD,jQ,...] >= _lims[1]))
                
                compliance[CCDkey][Q] = not np.any(test).sum()
        return compliance
    
    def check_stat_perCCDandCol(self,arr,lims,CCDs=[1,2,3]):
        """ """
        colnames = lims['CCD%i' % CCDs[0]].keys()
        
        compliance = OrderedDict()
        for iCCD,CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            compliance[CCDkey] = OrderedDict()
            for jcol,colname in enumerate(colnames):            
                _lims = lims[CCDkey][colname]
                ixsel = np.where(self.dd.mx['label'][:] == colname)
                
                test = (np.isnan(arr[ixsel,iCCD,...]) |\
                        (arr[ixsel,iCCD,...] <= _lims[0]) | (arr[ixsel,iCCD,...] >= _lims[1]))
                
                compliance[CCDkey][colname] = not (np.any(test,axis=(0,1)).sum() | (ixsel[0].shape[0]==0))
        return compliance
    
    def check_stat_perCCDQandCol(self,arr,lims,CCDs=[1,2,3]):
        """ """
        Qs = ['E','F','G','H']
        colnames = lims['CCD%i' % CCDs[0]][Qs[0]].keys()
        
        compliance = OrderedDict()
        for iCCD,CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            compliance[CCDkey] = OrderedDict()
            for iQ,Q in enumerate(Qs):
                compliance[CCDkey][Q] = OrderedDict()
                for jcol,colname in enumerate(colnames):            
                    _lims = lims[CCDkey][Q][colname]
                    ixsel = np.where(self.dd.mx['label'][:] == colname)
                    
                    test = (np.isnan(arr[ixsel,iCCD,...]) |\
                            (arr[ixsel,iCCD,...] <= _lims[0]) | (arr[ixsel,iCCD,...] >= _lims[1]))
                    
                    compliance[CCDkey][Q][colname] = not (np.any(test,axis=(0,1)).sum() | (ixsel[0].shape[0]==0))
        return compliance
    
    def check_data(self,**kwargs):
        """Generic check_data method"""
        if self.report is not None: 
            self.report.add_Section(keyword='check_data',Title='Data Validation',level=0)
        # CHECK AND CROSS-CHECK HK        
        self.check_HK_ST()
        # OBTAIN METRICS FROM IMAGES        
        self.get_checkstats_ST(**kwargs)
        # METRICS ASSESSMENT               
        self.check_metrics_ST(**kwargs)
        # PLOTs
        if self.report is not None: self.report.add_Section(keyword='check_plots',Title='Plots',level=1)
        self.addFigures_ST(**kwargs)
        # Update Report, raise flags, fill-in
        if self.log is not None:
            self.addFlagsToLog()
        if self.report is not None:
            self.addFlagsToReport()
        stop()
        
    def check_HK_ST(self):
        """ """
        HKKeys = self.HKKeys
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
        
    
    def addFigures_ST(self,**kwargs):
        """ """
        try: figkeys = kwargs['figkeys']
        except: figkeys=[]
        
        figspath = self.inputs['subpaths']['figs']
        
        for figkey in figkeys:
            try:
                pmeta = self.figdict[figkey][1]
                pmeta['path'] = figspath
                self.doPlot(figkey,**pmeta)
                self.addFigure2Report(figkey)
            except:
                self.catchtraceback()                
                nfigkey = 'BS_%s' % figkey                
                self.skipMissingPlot(nfigkey,ref=figkey)

            
    def prepare_images(self,doExtract=True,doMask=False,doOffset=False,doBias=False,
                       doFF=False):
        """ """
        
        if self.report is not None: 
            self.report.add_Section(keyword='prep_data',Title='Images Pre-Processing',level=0)
        
        if not doExtract: 
            self.report.add_Text('Not extracting FITS files: Nothing done.')
            return
        
        def _loadCDP(cdpkey,msg):
            CDPData = calibration.load_CDPs(self.inputs['inCDPs'][cdpkey],ccd.CCD)
            if self.log is not None:
                cdpstr = self.inputs['inCDPs'][cdpkey].__str__()
                cdpstr = st.replace(cdpstr,',',',\n')
                self.log.info(msg)
                self.log.info(cdpstr)
            return CDPData
        
        if doMask and 'mask' in self.inputs['inCDPs']:
            self.inputs['inCDPs']['Mask']['CCD%i']
            MaskData = _loadCDP('Mask','Applying cosmetics mask...')
            
        offsetkwargs = self.inputs['preprocessing']['offsetkwargs']

        if doBias and 'bias' in self.inputs['inCDPs']:
            BiasData = _loadCDP('Bias','Substracting Bias Structure...')

        
        if doFF and 'FF' in self.inputs['inCDPs']:
            FFData = _loadCDP('FF','Dividing by Flat-Field Map...')
        
        # Initialize new columns
    
        Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
        
        
        self.dd.initColumn('ccdobj_name',Cindices,dtype='S100',valini='None')
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs,nCCD,nQuad = DDindices.shape
        Quads = DDindices[2].vals
        CCDs = DDindices[DDindices.names.index('CCD')].vals
        
        if not self.drill:
            
            picklespath = self.inputs['subpaths']['ccdpickles']
            
            for iObs in range(nObs):
                
                if doFF:
                    FW_ID = self.dd.mx['wavelength'][iObs]
                    wavelength = ogse.FW['F%i' % FW_ID]
                
                for jCCD,CCD in enumerate(CCDs):
                    
                    CCDkey = 'CCD%i' % CCD
                    
                    ccdobj_name = '%s_proc' % self.dd.mx['File_name'][iObs,jCCD]
                    
                    dpath = self.dd.mx['datapath'][iObs,jCCD]
                    infits = os.path.join(dpath,'%s.fits' % self.dd.mx['File_name'][iObs,jCCD])
                    
                    ccdobj = ccd.CCD(infits) # converting the FITS file into a CCD Object
                    
                    fullccdobj_name = os.path.join(picklespath,'%s.pick' % self.dd.mx['ccdobj_name'][iObs,jCCD]) 
                    
                    if doMask:
                        ccdobj.get_mask(MaskData[CCDkey].extensions[-1])
                    
                    if doOffset:
                        for Quad in Quads:
                            #ccdobj.sub_offset(Quad,method='median',scan='pre',trimscan=[5,5],
                            #                  ignore_pover=False)
                            ccdobj.sub_offsets(Quad,**offsetkwargs)
                    
                    if doBias:
                        ccdobj.sub_bias(BiasData[CCDkey],extension=-1)
                    
                    if doFF:
                        FF = FFData['nm%i' % wavelength][CCDkey]
                        ccdobj.divide_by_flatfield(FF,extension=-1)
                    
                    ccdobj.writeto(fullccdobj_name,clobber=True)                    
                    self.dd.mx['ccdobj_name'][iObs,jCCD] = ccdobj_name
    
        
        return None
    
    def add_inputs_to_report(self):
        """ """
        
        self.report.add_Text('\\textbf{Test Inputs}')
        
        caption = 'Inputs of Task %s , Test $%s$' % (self.name,self.inputs['test'])
        #ncaption = st.replace(caption,'_','\\_')
        
        names = ['Parameter','Value']
        
        keys = self.inputs.manifesto.keys()
        values = []
        for key in keys:
            _val = self.inputs[key]
            if isinstance(_val,dict):
                values.append('dict()')
            else:
                _val = _val.__repr__()
                #n_val = st.replace(_val,'_','\\_')
                n_val = st.replace(_val,'&','\\&')
                values.append(n_val)
                
        tDict = OrderedDict(Parameter=keys,Value=values)
        #formats = dict(Parameter='s',Value='char')
        
        self.report.add_Table(tDict,names=names,caption=caption)
        
