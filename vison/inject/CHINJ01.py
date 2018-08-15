#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ01

Charge injection calibration (part 1)
    Injection vs. IG1-IG2

Created on Tue Aug 29 17:36:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import pandas as pd

from vison.support import context
from vison.datamodel import scriptic as sc
from vison.support import utils
from vison.datamodel import inputs
from InjTask import InjTask
import CH01aux
from vison.support import files
from vison.inject import lib as ilib
from vison.datamodel import cdp
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

CHINJ01_commvalues = dict(program='CALCAMP', test='CHINJ01',
                          IPHI1=1, IPHI2=1, IPHI3=1, IPHI4=0,
                          rdmode='fwd_bas',
                          flushes=7, exptime=0., vstart=0, vend=2086,
                          siflsh=1, siflsh_p=500,
                          shuttr=0,
                          chinj=1, chinj_on=30,
                          chinj_of=100,
                          id_wid=60, chin_dly=1,
                          operator='who',
                          comments='')


class CHINJ01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('IDL', ([float], 'Injection Drain Low Voltage.')),
        ('IDH', ([float], 'Injection Drain High Voltage.')),
        ('IG2', ([float], 'Injection Gate 2 Voltage.')),
        ('IG1s', ([list], 'Injection Gate 1 Voltages, [min, max].')),
        ('dIG1', ([float], 'Injection Gate 1 Voltage Step.')),
        ('id_delays', ([list], 'Injection Drain Delays.')),
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
    ])))


        
    
class CHINJ01(InjTask):
    """ """

    inputsclass = CHINJ01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(CHINJ01, self).__init__(inputs, log, drill, debug)
        self.name = 'CHINJ01'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = CH01aux.CH01figs.copy()
        self.inputs['subpaths'] = dict(figs='figs',
                   ccdpickles='ccdpickles',
                   products='products')
        

    def set_inpdefaults(self, **kwargs):
        """ """
        toi_chinj = 500

        self.inpdefaults = dict(
            IDL=11.,
            IDH=18.,
            IG2=5.5,
            IG1s=[2., 6.],
            dIG1=0.25,
            id_delays=[toi_chinj*2.5, toi_chinj*1.5],
            toi_chinj=toi_chinj
        )

    def set_perfdefaults(self, **kwargs):
        super(CHINJ01, self).set_perfdefaults(**kwargs)

        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()

        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()
        
    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """
        Builds CHINJ01 script structure dictionary.
        
        #:param IDL: float, [mV], value of IDL (Inject. Drain Low).
        #:param IDH: float, [mV], Injection Drain High.
        #:param IG2: float, [mV], Injection Gate 2.
        #:param IG1s: list of 2 floats, [mV], [min,max] values of IG1.
        #:param id_delays: list of 2 floats, [mV], injection drain delays (2).
        #:param toi_chinj: int, [us], TOI-charge injection.
        :param diffvalues: dict, opt, differential values.
        
        """
        
        IDL = self.inputs['IDL']
        IDH = self.inputs['IDH']
        IG2 = self.inputs['IG2']
        IG1s = self.inputs['IG1s']
        dIG1 = self.inputs['dIG1'] # 0.25  # V
        id_delays = self.inputs['id_delays']
        toi_chinj = self.inputs['toi_chinj']
        
        CCDs = [1, 2, 3]
        halves = ['T', 'B']
        
        assert len(IG1s) == 2
        assert len(id_delays) == 2
        
        NIG1 = (IG1s[1]-IG1s[0])/dIG1+1
        IG1v = np.arange(NIG1)*dIG1+IG1s[0]
        
        CHINJ01_sdict = dict()
        
        # First Injection Drain Delay
        
        colcounter = 1
        for i, IG1 in enumerate(IG1v):
            colkey = 'col%i' % colcounter
            #print colkey
            CHINJ01_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                         IG2_T=IG2, IG2_B=IG2,
                                         id_dly=id_delays[0], toi_ch=toi_chinj)
           
            for CCD in CCDs:
                for half in halves:
                    CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD, half)] = IG1

            colcounter += 1

        # Second Injection Drain Delay

        for j, IG1 in enumerate(IG1v):
            colkey = 'col%i' % colcounter
            #print colkey
            CHINJ01_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                         IG2_T=IG2,IG2_B=IG2,
                                         id_dly=id_delays[1], toi_ch=toi_chinj)

            for CCD in CCDs:
                for half in halves:
                    CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD, half)] = IG1

            colcounter += 1

        Ncols = len(CHINJ01_sdict.keys())
        CHINJ01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(CHINJ01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        CHINJ01_sdict = sc.update_structdict(
            CHINJ01_sdict, commvalues, diffvalues)

        return CHINJ01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(CHINJ01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                    wavedkeys=wavedkeys)
    
    def prepare_images(self):
        super(CHINJ01, self).prepare_images(doExtract=True, 
             doMask=True, # ON TESTS!
             doOffset=True, 
             doBias=False, doFF=False)

    def basic_analysis(self):
        self.BROKEN_basic_analysis()

    def old_basic_analysis(self):
        """ 

        Basic analysis of data.

        **METACODE**

        ::

            f. e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        extract average 2D injection pattern (and save)
                        produce average profile along/across lines
                        measure charge-inj. non-uniformity
                        measure charge spillover into non-injection
                        measure stats of injection (mean, med, std, min/max, percentiles)

            plot average inj. profiles along lines f. each CCD, Q and IG1
                save as a rationalized set of curves
            plot average inj. profiles across lines f. each CCD, Q and IG1
                save as a rationalized set of  curves
       
            Report injection stats as a table/tables

        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='CHINJ01 Extraction', level=0)
        
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuad = DDindices.shape[0:3]
        Quads = DDindices.get_vals('Quad')
        CCDs = DDindices.get_vals('CCD')
        
        ccdpicklespath = self.inputs['subpaths']['ccdpickles']
        prodspath = self.inputs['subpaths']['products']
        
        # Initializing new columns
        

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        valini = 0.
        
        # measure charge-inj. non-uniformity
        # measure charge spillover into non-injection
        # measure stats of injection (mean, med, std, min/max, percentiles)

        self.dd.initColumn('chinj_nonuni', DDindices, dtype='float32', valini=valini)
        self.dd.initColumn('chinj_spill', DDindices, dtype='float32', valini=valini)
        statkeys = ['mean','std','min','max','p5','p25','p50','p75',
                    'p95']
        for statkey in statkeys:
             self.dd.initColumn('chinj_%s' % statkey, DDindices, dtype='float32', valini=valini)

        # EXTRACTION TABLE
        
        NP = nObs * nCCD * nQuad
        
        # OBSID CCD Q IG1 id_dly MEAN MEDIAN NONUNI
        
        CH01_dd = OrderedDict()
        CH01_dd['ObsID'] = np.zeros(NP,dtype='int32')
        CH01_dd['CCD'] = np.zeros(NP,dtype='int32')
        CH01_dd['Q'] = np.zeros(NP,dtype='int32')
        CH01_dd['IG1'] = np.zeros(NP,dtype='float32')
        CH01_dd['ID_DLY'] = np.zeros(NP,dtype='float32')
        CH01_dd['MEAN_INJ'] = np.zeros(NP,dtype='float32')
        CH01_dd['MED_INJ'] = np.zeros(NP,dtype='float32')
        CH01_dd['NU_INJ'] = np.zeros(NP,dtype='float32')
        
        # Initializing injection profiles
        
        
        prof_alrow_cdp = cdp.CDP()
        prof_alrow_cdp.header = CDP_header.copy()
        prof_alrow_cdp.path = prodspath
        prof_alrow_cdp.data = OrderedDict()
        
        prof_alcol_cdp = cdp.CDP()
        prof_alcol_cdp.header = CDP_header.copy()
        prof_alcol_cdp.path = prodspath
        prof_alcol_cdp.data = OrderedDict()
        
        xdummy = np.arange(10,dtype='float32')
        ydummy = np.zeros(10,dtype='float32')
        for jCCD, CCDk in enumerate(CCDs):
            prof_alrow_cdp.data[CCDk] = OrderedDict()
            prof_alcol_cdp.data[CCDk] = OrderedDict()
            
            for kQ, Q in enumerate(Quads):
                prof_alrow_cdp.data[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                      y=OrderedDict())
                prof_alcol_cdp.data[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                  y=OrderedDict())
                
                for iObs in range(nObs):
                    
                    IG1_key = 'IG1_%i_%s' % (jCCD+1,_get_CCDhalf(Q))
                    IG1_val = self.dd.mx[IG1_key][iObs,jCCD]
                    IG1_tag = 'IG1_%.2fV' % IG1_val
                    prof_alrow_cdp.data[CCDk][Q]['x'][IG1_tag] = xdummy.copy()
                    prof_alrow_cdp.data[CCDk][Q]['y'][IG1_tag] = ydummy.copy()
                    prof_alcol_cdp.data[CCDk][Q]['x'][IG1_tag] = xdummy.copy()
                    prof_alcol_cdp.data[CCDk][Q]['y'][IG1_tag] = ydummy.copy()
        
        
        # The hardwork
        
        if not self.drill:
            
            for iObs in range(nObs):
                
                ObsID = self.dd.mx['ObsID'][iObs]
                
                for jCCD, CCDk in enumerate(CCDs):
                    
                    ccdpickf = os.path.join(ccdpicklespath,\
                                            '%s.pick' % self.dd.mx['ccdobj_name'][iObs,jCCD])
                    
                    ccdobj = files.cPickleRead(ccdpickf)
                    
                    vstart = self.dd.mx['vstart'][iObs,jCCD]
                    vend = self.dd.mx['vend'][iObs,jCCD]
                    dochinj = self.dd.mx['chinj'][iObs,jCCD]
                    non = self.dd.mx['chinj_on'][iObs,jCCD]
                    noff = self.dd.mx['chinj_of'][iObs,jCCD]
                    nrep = (vend-vstart)/(non+noff)+1
                    pattern = (non, noff, nrep)
                    
                    for kQ, Q in enumerate(Quads):
                        
                        ix = iObs * nCCD * nQuad + jCCD * nQuad + kQ
                    
                        if dochinj:
                            
                            IG1_key = 'IG1_%i_%s' % (jCCD+1,_get_CCDhalf(Q))
                            IG1_val = self.dd.mx[IG1_key][iObs,jCCD]
                            IG1_tag = 'IG1_%.2fV' % IG1_val
                            
                            
                            id_dly = self.dd.mx['id_dly'][iObs,jCCD]
                            
                            ext_res = ilib.extract_injection_lines(ccdobj, Q, pattern, VSTART=vstart,
                                                                  VEND=vend, suboffmean=False)
                            
                            stats = ext_res['stats_injection']
                            
                            for skey in statkeys:
                                self.dd.mx['chinj_%s' % skey][iObs,jCCD,kQ] = \
                                       stats[skey]
                            
                            
                            nonuni = (stats['p95']-stats['p5'])/stats['p50']
                            
                            self.dd.mx['chinj_nonuni'][iObs,jCCD,kQ] = nonuni
                            
                            spill = ilib.get_spill(ext_res['avprof_alcol'],pattern)
                            
                            self.dd.mx['chinj_spill'][iObs,jCCD,kQ] = spill
                            
                                      
                            yalrows = ext_res['avprof_alrow'].copy()
                            yalcols = ext_res['avprof_alcol'].copy()
                            prof_alrow_cdp.data[CCDk][Q]['y'][IG1_tag] = yalrows.copy()
                            prof_alrow_cdp.data[CCDk][Q]['x'][IG1_tag] = \
                                          np.arange(len(yalrows),dtype='float32')
                            
                            prof_alcol_cdp.data[CCDk][Q]['y'][IG1_tag] = yalcols.copy()
                            prof_alcol_cdp.data[CCDk][Q]['x'][IG1_tag] = \
                                          np.arange(len(yalcols),dtype='float32')
                            
                            
                            CH01_dd['ObsID'][ix] = ObsID
                            CH01_dd['CCD'][ix] = jCCD
                            CH01_dd['Q'][ix] = kQ
                            CH01_dd['IG1'][ix] = IG1_val
                            CH01_dd['ID_DLY'][ix] = id_dly
                            CH01_dd['MEAN_INJ'][ix] = self.dd.mx['chinj_mean'][iObs,jCCD,kQ]
                            CH01_dd['MED_INJ'][ix] = self.dd.mx['chinj_p50'][iObs,jCCD,kQ]
                            CH01_dd['NU_INJ'][ix] = self.dd.mx['chinj_nonuni'][iObs,jCCD,kQ]
                            
                            
        # plot average inj. profiles along/across lines 
        # save as a rationalized set of curves
        
        maxmedinjection = np.nanmax(self.dd.mx['chinj_p50'][:])

        fdict_alrow = self.figdict['CH01_alrow'][1]
        fdict_alrow['data'] = prof_alrow_cdp.data.copy()
        fdict_alrow['meta']['ylim'] = [0.,maxmedinjection*1.5]

        
        fdict_alcol = self.figdict['CH01_alcol'][1]
        fdict_alcol['data'] = prof_alcol_cdp.data.copy()
        fdict_alcol['meta']['ylim'] = [0.,maxmedinjection*1.5]
        
        self.pack_CDP_to_dd(prof_alrow_cdp,'PROFS_ALROW')
        self.pack_CDP_to_dd(prof_alcol_cdp,'PROFS_ALCOL')
        
        if self.report is not None:
            self.addFigures_ST(figkeys=['CH01_alrow',
                                        'CH01_alcol'], 
                               dobuilddata=False)
                
        # Report injection stats as a table/tables
        
        # OBSID CCD Q IG1 id_dly MEAN MEDIAN NONUNI
        
        EXT_dddf = OrderedDict(EXTRACT=pd.DataFrame.from_dict(CH01_dd))
        EXT_cdp = CH01aux.CDP_lib['EXTRACT']
        EXT_cdp.path = prodspath
        EXT_cdp.ingest_inputs(
                data=EXT_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy())
        
        EXT_cdp.init_wb_and_fillAll(header_title='CHINJ01: EXTRACTION')
        self.save_CDP(EXT_cdp)
        self.pack_CDP_to_dd(EXT_cdp, 'EXTRACT_CDP')
        
        if self.report is not None:
            fi = lambda x: '%i' % x
            fccd = lambda x: CCDs[x-1]
            fq = lambda x: Quads[x-1]
            ff = lambda x: '%.2f' % x
            
            ext_formatters=[fi,fccd,fq,ff,ff,ff,ff,ff]
            
            caption = 'CHINJ01: EXTRACTION TABLE' 
            Etex = EXT_cdp.get_textable(sheet='EXTRACT', caption=caption,
                                               fitwidth=True,
                                               formatters=ext_formatters)
            
            Etex = ['\\tiny']+Etex+['\\normalsize']
            self.report.add_Text(Etex)
        
        
    def meta_analysis(self):
        """ 
        
        Find injection threshold: Min IG1
        Plot and model charge injection vs. IG1
        Find notch injection amount.
        
        
        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='CHINJ01 Analysis ("Meta")', level=0)
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuad = DDindices.shape[0:3]
        Quads = DDindices.get_vals('Quad')
        CCDs = DDindices.get_vals('CCD')
        
        prodspath = self.inputs['subpaths']['products']
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        
        
        toi_ch = self.dd.mx['toi_ch'][0,0]
        assert np.all(np.isclose(self.dd.mx['toi_ch'][:],toi_ch))
        
        
        # ANALYSIS TABLE
        
        NP = nCCD * nQuad
        
        MCH01_dd = OrderedDict()
        MCH01_dd['CCD'] = np.zeros(NP,dtype='int32')
        MCH01_dd['Q'] = np.zeros(NP,dtype='int32') 
        MCH01_dd['ID_DLY'] = np.zeros(NP,dtype='float32') + np.nan
        
        fitkeys = ['BGD','DROP','IG1_THRESH','IG1_NOTCH','SLOPE','NOTCH']
                
        for fitkey in fitkeys:
            MCH01_dd[fitkey] = np.zeros(NP,dtype='float32') + np.nan
        
        
        # INJECTION CURVES
        
        inj_curves_cdp = cdp.CDP()
        inj_curves_cdp.header = CDP_header.copy()
        inj_curves_cdp.path = prodspath
        inj_curves_cdp.data = OrderedDict()
        
        xdummy = np.arange(10,dtype='float32')
        ydummy = np.zeros(10,dtype='float32')
        for jCCD, CCDk in enumerate(CCDs):
            
            inj_curves_cdp.data[CCDk] = OrderedDict()
            
            for kQ, Q in enumerate(Quads):
                inj_curves_cdp.data[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                      y=OrderedDict())
                
                for tag in ['input','bestfit']:
                    inj_curves_cdp.data[CCDk][Q]['x'][tag] = xdummy.copy()
                    inj_curves_cdp.data[CCDk][Q]['y'][tag] = ydummy.copy()
        
        
        for jCCD, CCDk in enumerate(CCDs):
            
            for kQ, Q in enumerate(Quads):
                
                ix = jCCD * nQuad + kQ
                
                
                if Q in ['E','F']:
                    id_dly_opt = toi_ch * 2.5
                elif Q in['G','H']:
                    id_dly_opt = toi_ch * 1.5
                                
                selix = np.where((self.dd.mx['chinj'][:,jCCD]==1) &
                        (np.isclose(self.dd.mx['id_dly'][:,jCCD],id_dly_opt)))
                
                
                CCDhalf = _get_CCDhalf(Q)
                IG1_key = 'IG1_%i_%s' % (jCCD+1,CCDhalf)
                
                IG1 = self.dd.mx[IG1_key][selix,jCCD].flatten().copy()                
                med_inj = self.dd.mx['chinj_p50'][selix,jCCD,kQ].flatten().copy()
                
                res = ilib.fit_Inj_vs_IG1(IG1,med_inj,doPlot=False)
                didfit = res['didfit']
                
                inj_curves_cdp.data[CCDk][Q]['x']['input'] = IG1.copy()
                inj_curves_cdp.data[CCDk][Q]['y']['input'] = med_inj.copy() / 2.**16
                
                
                if didfit:
                
                    MCH01_dd['CCD'][ix] = jCCD
                    MCH01_dd['Q'][ix] = kQ
                    MCH01_dd['ID_DLY'][ix] = id_dly_opt
                    
                    for fitkey in fitkeys:                    
                        MCH01_dd[fitkey][ix] = res[fitkey]

                    xbf = res['IG1_BF'].copy()
                    ybf = res['NORMINJ_BF']                        
                    
                else:
                    
                    xbf = IG1.copy()
                    ybf = np.zeros_like(xbf)
                
                inj_curves_cdp.data[CCDk][Q]['x']['bestfit'] = xbf.copy()
                inj_curves_cdp.data[CCDk][Q]['y']['bestfit'] = ybf.copy()
                    
        
        # PLOT
        
        fdict_meta_plot = self.figdict['CH01_meta'][1]
        fdict_meta_plot['data'] = inj_curves_cdp.data.copy()
        
        if self.report is not None:
            self.addFigures_ST(figkeys=['CH01_meta'], 
                               dobuilddata=False)
        
        # REPORT RESULTS AS TABLE CDP
        
        
        MCH01_dddf = OrderedDict(ANALYSIS=pd.DataFrame.from_dict(MCH01_dd))
        MCH01_cdp = CH01aux.CDP_lib['META']
        MCH01_cdp.path = prodspath
        MCH01_cdp.ingest_inputs(
                data=MCH01_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )
        
        MCH01_cdp.init_wb_and_fillAll(header_title='CHINJ01: META-ANALYSIS')
        self.save_CDP(MCH01_cdp)
        self.pack_CDP_to_dd(MCH01_cdp,'META_CDP')
        
        if self.report is not None:
            
            fccd = lambda x: CCDs[x-1]
            fq = lambda x: Quads[x-1]
            ff = lambda x: '%.2f' % x
            
            selcolumns = ['CCD','Q','BGD','IG1_THRESH','IG1_NOTCH','NOTCH']
            
            ext_formatters=[fccd,fq]+[ff,ff,ff,ff]
            
            caption = 'CHINJ01: META-ANALYSIS TABLE'
            
            Mtex = MCH01_cdp.get_textable(sheet='ANALYSIS', 
                                          columns=selcolumns,
                                          caption=caption,
                                          fitwidth=True,
                                          formatters=ext_formatters)
            
            Mtex = ['\\tiny']+Mtex+['\\normalsize']
            self.report.add_Text(Mtex)  
        
        

