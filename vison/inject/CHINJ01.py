#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ01

Charge injection calibration (part 1)
    Injection vs. IG1-IG2

Created on Tue Aug 29 17:36:00 2017

:author: Ruyman Azzollini

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
from InjTask import InjTask, _get_CCDhalf
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
                          flushes=7, siflsh=1, siflsh_p=500,
                          inisweep=1,
                          vstart=0, vend=2086,
                          chinj=1, chinj_on=30, chinj_of=100,
                          id_wid=60,
                          exptime=0., shuttr=0, e_shuttr=0,
                          mirr_on=0,
                          wave=4,
                          motr_on=0,
                          source='flat',
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
        self.CDP_lib = CH01aux.CDP_lib
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
            IG2=6.5,
            IG1s=[2., 7.],
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
        
        #:param IDL: float, [V], value of IDL (Inject. Drain Low).
        #:param IDH: float, [V], Injection Drain High.
        #:param IG2: float, [V], Injection Gate 2.
        #:param IG1s: list of 2 floats, [V], [min,max] values of IG1.
        #:param id_delays: list of 2 floats, [us], injection drain delays.
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
            colkey = 'col%03i' % colcounter
            #print colkey
            CHINJ01_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                         IG2_T=IG2, IG2_B=IG2,
                                         id_dly=id_delays[0], 
                                         toi_ch=toi_chinj)
           
            for CCD in CCDs:
                for half in halves:
                    CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD, half)] = IG1

            colcounter += 1

        # Second Injection Drain Delay

        for j, IG1 in enumerate(IG1v):
            colkey = 'col%03i' % colcounter
            #print colkey
            CHINJ01_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                         IG2_T=IG2,IG2_B=IG2,
                                         id_dly=id_delays[1], 
                                         toi_ch=toi_chinj)

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

        
    def meta_analysis(self):
        """ 
        
        Plot and model charge injection vs. IG1
        Find injection threshold: Min IG1
        Find notch injection amount.
        
        
        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='CHINJ01 Analysis ("Meta")', level=0)
        
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
        inj_curves_cdp.data['labelkeys'] = ['input','bestfit']
        
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
                
                MCH01_dd['CCD'][ix] = jCCD+1
                MCH01_dd['Q'][ix] = kQ+1
                
                if didfit:
                
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
        MCH01_cdp = self.CDP_lib['META']
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
            
            selcolumns = ['CCD','Q','BGD','IG1_THRESH','IG1_NOTCH','SLOPE','NOTCH']
            
            ext_formatters=[fccd,fq]+[ff,ff,ff,ff,ff]
            
            caption = 'CHINJ01: META-ANALYSIS TABLE'
            
            Mtex = MCH01_cdp.get_textable(sheet='ANALYSIS', 
                                          columns=selcolumns,
                                          caption=caption,
                                          fitwidth=True,
                                          formatters=ext_formatters,
                                          index=False)
            
            Mtex = ['\\tiny']+Mtex+['\\normalsize']
            self.report.add_Text(Mtex)  
        
        

