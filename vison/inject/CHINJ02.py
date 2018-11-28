#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ02

Charge injection calibration (part 2)
    Injection vs. IDL (injection threshold)

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

from vison.pipe.task import HKKeys
from vison.support import context, utils
from vison.datamodel import scriptic as sc
from vison.datamodel import inputs
from InjTask import InjTask
from vison.image import performance
import CH02aux
# END IMPORT

isthere = os.path.exists

#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

IG1comm = 6.
IG2comm = 4.

CHINJ02_commvalues = dict(program='CALCAMP', test='CHINJ02',                          
                          flushes=7, siflsh=1,siflsh_p=500,
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


class CHINJ02_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('IDLs', ([list], 'Injection Drain Low Voltages List: [min, max].')),
        ('dIDL', ([float], 'Injection Drain Voltage Step.')),
        ('IDH', ([float], 'Injection Drain High Voltage.')),
        ('IG1', ([float], 'Injection Gate 1 Voltage.')),
        ('IG2', ([float], 'Injection Gate 2 Voltage.')),
        ('id_delays', ([list], 'Injection Drain Delays.')),
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
    ])))


class CHINJ02(InjTask):
    """ """

    inputsclass = CHINJ02_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(CHINJ02, self).__init__(inputs, log, drill, debug)
        self.name = 'CHINJ02'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.CDP_lib = CH02aux.get_CDP_lib()
        self.figdict = CH02aux.get_CH02figs()
        self.inputs['subpaths'] = dict(figs='figs',
                   ccdpickles='ccdpickles',
                   products='products')
        

    def set_inpdefaults(self, **kwargs):
        """ """
        toi_chinj = 500

        self.inpdefaults = dict(
            IDLs=[10., 13.],
            dIDL=0.25,
            IG1=IG1comm,
            IG2=IG2comm,
            IDH=18.,
            id_delays=[toi_chinj*2.5, toi_chinj*1.5],
            toi_chinj=toi_chinj
        )

    def set_perfdefaults(self, **kwargs):
        super(CHINJ02, self).set_perfdefaults(**kwargs)

        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()

        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ 
        Builds CHINJ02 script structure dictionary.

        #:param IDLs: list of 2 ints, [V], [min,max] values of IDL (Inject. Drain Low).
        #:param IDH: int, [V], Injection Drain High.
        #:param id_delays: list of 2 ints, [us], injection drain delays.
        #:param toi_chinj: int, [us], TOI-charge injection.
        :param diffvalues: dict, opt, differential values.

        """

        IDLs = self.inputs['IDLs']
        dIDL = self.inputs['dIDL']
        IDH = self.inputs['IDH']
        IG1 = self.inputs['IG1']
        IG2 = self.inputs['IG2']
        id_delays = self.inputs['id_delays']
        toi_chinj = self.inputs['toi_chinj']

        assert len(IDLs) == 2
        assert len(id_delays) == 2

        NIDL = (IDLs[1]-IDLs[0])/dIDL+1
        IDLv = np.arange(NIDL)*dIDL+IDLs[0]

        CHINJ02_sdict = dict()

        # First Injection Drain Delay

        colcounter = 1
        for i, IDL in enumerate(IDLv):
            colkey = 'col%03i' % (i+1,)
            CHINJ02_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                         IG1_1_T=IG1, IG1_2_T=IG1, IG1_3_T=IG1,
                                         IG1_1_B=IG1, IG1_2_B=IG1, IG1_3_B=IG1,
                                         IG2_T=IG2, IG2_B=IG2,
                                         id_dly=id_delays[0], toi_ch=toi_chinj)
            colcounter += 1

        # Second Injection Drain Delay

        colstart = colcounter

        for j, IDL in enumerate(IDLv):
            colkey = 'col%03i' % (colstart+j,)
            CHINJ02_sdict[colkey] = dict(frames=1, IDL=IDL, IDH=IDH,
                                         IG1_1_T=IG1, IG1_2_T=IG1, IG1_3_T=IG1,
                                         IG1_1_B=IG1, IG1_2_B=IG1, IG1_3_B=IG1,
                                         IG2_T=IG2, IG2_B=IG2,
                                         id_dly=id_delays[1], toi_ch=toi_chinj)

        Ncols = len(CHINJ02_sdict.keys())
        CHINJ02_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(CHINJ02_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        CHINJ02_sdict = sc.update_structdict(
            CHINJ02_sdict, commvalues, diffvalues)

        return CHINJ02_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(CHINJ02, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                    wavedkeys=wavedkeys)


    def meta_analysis(self):
        """ 

        Finds the Injection Threshold for each CCD half.

        **METACODE**

        ::

            f.e.CCD:
                f.e.Q:
                    load injection vs. IDL cuve
                    find&save injection threshold on curve

            report injection threshold as a table

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='CHINJ02 Analysis ("Meta")', level=0)
        
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
        
        CCDhalves = ['top','bottom']
        
        NP = nObs
        
        MCH02_dd = OrderedDict()
        MCH02_dd['meta'] = OrderedDict()
        
        _Quads_dict = dict(top = ['E','F'],
                           bottom = ['G','H'])
        
        for CCDhalf in CCDhalves:
            MCH02_dd[CCDhalf] = OrderedDict()
            MCH02_dd['meta'][CCDhalf] = OrderedDict()
            #MCH02_dd[CCDhalf]['CCD'] = np.zeros(NP,dtype='int32')
            #MCH02_dd[CCDhalf]['Q'] = np.zeros(NP,dtype='int32')
            
            MCH02_dd[CCDhalf]['IDL'] = np.zeros(NP,dtype='float32') + np.nan
            
            for jCCD,CCDk in enumerate(CCDs):
                for Q in _Quads_dict[CCDhalf]:                    
                    MCH02_dd[CCDhalf]['INJ_%s%s' % (jCCD+1,Q)] = np.zeros(NP,dtype='float32') + np.nan
            
        
        for CCDhalf in CCDhalves:
            
            _Quads = _Quads_dict[CCDhalf]
            
            if CCDhalf == 'top':
                id_dly_opt = toi_ch * 2.5                
            elif CCDhalf == 'bottom':
                id_dly_opt = toi_ch * 1.5
                
            for jCCD, CCDk in enumerate(CCDs):
            
                selix = np.where((self.dd.mx['chinj'][:,jCCD]==1) &
                            (np.isclose(self.dd.mx['id_dly'][:,jCCD],id_dly_opt)))
                                  
                IDL = self.dd.mx['IDL'][selix,jCCD].flatten().copy()
                 
                for Q in _Quads:
                      
                    kQ = Quads.index(Q)
                     
                    med_inj = self.dd.mx['chinj_p50'][selix,jCCD,kQ].flatten().copy()
             
                    MCH02_dd[CCDhalf]['IDL'] = IDL.copy()
                    MCH02_dd[CCDhalf]['INJ_%s%s' % (jCCD+1, Q)] = med_inj.copy()
                    
                
                MCH02_dd['meta'][CCDhalf]['toi_ch'] = toi_ch
                MCH02_dd['meta'][CCDhalf]['id_dly'] = id_dly_opt
                MCH02_dd['meta'][CCDhalf]['IDH'] = self.dd.mx['IDH'][selix,jCCD].flatten()[0]
                IG1key = 'IG1_1_%s' % CCDhalf[0].upper()
                MCH02_dd['meta'][CCDhalf]['IG1'] = self.dd.mx[IG1key][selix,jCCD].flatten()[0]
                IG2key = 'IG2_%s' % CCDhalf[0].upper()
                MCH02_dd['meta'][CCDhalf]['IG2'] = self.dd.mx[IG2key][selix,jCCD].flatten()[0]
                
        
        # REPORT RESULTS AS TABLE CDPs 
        
        for CCDhalf in CCDhalves:
            
            MCH02half_dddf = OrderedDict(ANALYSIS=pd.DataFrame.from_dict(MCH02_dd[CCDhalf]))
            MCH02half_cdp = self.CDP_lib['META']
            MCH02half_cdp.path = prodspath
            MCH02half_cdp.rootname += '_%s' % CCDhalf
            MCH02half_cdp.ingest_inputs(
                    data=MCH02half_dddf.copy(),
                    meta=MCH02_dd['meta'][CCDhalf].copy(),
                    header=CDP_header.copy()
                    )
            
            MCH02half_cdp.init_wb_and_fillAll(header_title='CHINJ02: META-ANALYSIS [%s]' % CCDhalf)
            self.save_CDP(MCH02half_cdp)
            self.pack_CDP_to_dd(MCH02half_cdp,'META_CDP_%s' % CCDhalf)
            
            
            if self.report is not None:
                
                ff = lambda x: '%.2f' % x
                
                selcolumns = ['IDL']
                
                for jCCD,CCDk in enumerate(CCDs):
                    for Q in _Quads_dict[CCDhalf]:                    
                        selcolumns.append('INJ_%s%s' % (jCCD+1,Q)) 
                
                ext_formatters=[ff]*len(selcolumns)
                
                caption = 'CHINJ02: META-ANALYSIS TABLE. CCD-half = %s. '+\
                 'IDH = %.2f V,  IG1 = %.2f V, IG2 = %.2f V, toi\_chj = % us, '+\
                 'id\_dly=%.1f us.'
                caption = caption % (CCDhalf, 
                                     MCH02half_cdp.meta['IDH'], 
                                     MCH02half_cdp.meta['IG1'], 
                                     MCH02half_cdp.meta['IG2'], 
                                     MCH02half_cdp.meta['toi_ch'], 
                                     MCH02half_cdp.meta['id_dly'])
                
                Mtex = MCH02half_cdp.get_textable(sheet='ANALYSIS', 
                                              columns=selcolumns,
                                              caption=caption,
                                              fitwidth=True,
                                              tiny=True,
                                              formatters=ext_formatters,
                                              index=False)
                
                self.report.add_Text(Mtex)  
        
