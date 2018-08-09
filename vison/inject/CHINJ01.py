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

from vison.support import context
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#import datetime
from vison.datamodel import inputs
from InjTask import InjTask
from vison.image import performance
import CH01aux
from vison.support import files
from vison.inject import lib as ilib
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

CHINJ01_commvalues = dict(program='CALCAMP', test='CHINJ01',
                          IG2_T=5.5, IG2_B=5.5,
                          IPHI1=1, IPHI2=1, IPHI3=1, IPHI4=0,
                          rdmode='fwd_bas',
                          flushes=7, exptime=0., vstart=0, vend=2086,
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
        ('IG1s', ([list], 'Injection Gate 1 Voltages.')),
        ('dIG1', ([float], 'Injection Gate 1 Voltage Step.')),
        ('id_delays', ([list], 'Injection Drain Delays.')),
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
    ])))


def _get_CCDhalf(Q):
    if Q in ['E','F']:
        return 'B'
    elif Q in ['G','H']:
        return 'T'
        
    
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

        #:param IDL: int, [mV], value of IDL (Inject. Drain Low).
        #:param IDH: int, [mV], Injection Drain High.
        #:param IG1s: list of 2 ints, [mV], [min,max] values of IG1.
        #:param id_delays: list of 2 ints, [mV], injection drain delays (2).
        #:param toi_chinj: int, [us], TOI-charge injection.
        :param diffvalues: dict, opt, differential values.

        """

        IDL = self.inputs['IDL']
        IDH = self.inputs['IDH']
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
        
        # Initializing new columns

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
        
        profiles_alrow = OrderedDict()
        profiles_alcol = OrderedDict()
        
        
        # Initializing injection profiles
        
        xdummy = np.arange(10,dtype='float32')
        ydummy = np.zeros(10,dtype='float32')
        for jCCD, CCDk in enumerate(CCDs):
            profiles_alrow[CCDk] = OrderedDict()
            profiles_alcol[CCDk] = OrderedDict()
            
            for kQ, Q in enumerate(Quads):
                profiles_alrow[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                      y=OrderedDict())
                profiles_alcol[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                  y=OrderedDict())
                
                for iObs in range(nObs):
                    
                    IG1_key = 'IG1_%i_%s' % (jCCD+1,_get_CCDhalf(Q))
                    IG1_val = self.dd.mx[IG1_key][iObs,jCCD]
                    IG1_tag = 'IG1_%.2fV' % IG1_val
                    profiles_alrow[CCDk][Q]['x'][IG1_tag] = xdummy.copy()
                    profiles_alrow[CCDk][Q]['y'][IG1_tag] = ydummy.copy()
                    profiles_alcol[CCDk][Q]['x'][IG1_tag] = xdummy.copy()
                    profiles_alcol[CCDk][Q]['y'][IG1_tag] = ydummy.copy()
        
        
        ccdpicklespath = self.inputs['subpaths']['ccdpickles']
        prodspath = self.inputs['subpaths']['products']
        
        # The hardwork
        
        if not self.drill:
            
            
            for iObs in range(nObs):
                
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
                    
                    if dochinj:
                    
                        for kQ, Q in enumerate(Quads):
                            
                            IG1_key = 'IG1_%i_%s' % (jCCD+1,_get_CCDhalf(Q))
                            IG1_val = self.dd.mx[IG1_key][iObs,jCCD]
                            IG1_tag = 'IG1_%.2fV' % IG1_val
                            
                            
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
                            profiles_alrow[CCDk][Q]['y'][IG1_tag] = yalrows.copy()
                            profiles_alrow[CCDk][Q]['x'][IG1_tag] = \
                                          np.arange(len(yalrows),dtype='float32')
                            
                            profiles_alcol[CCDk][Q]['x'][IG1_tag] = yalcols.copy()
                            profiles_alcol[CCDk][Q]['y'][IG1_tag] = \
                                          np.arange(len(yalrows),dtype='float32')
                            
                            
        # plot average inj. profiles along lines 
        # save as a rationalized set of curves
        
        # PENDING
        
        # plot average inj. profiles across lines
        # save as a rationalized set of  curves
        
        # PENDING
        
        # Report injection stats as a table/tables
        
        # PENDING
                    
        
    def meta_analysis(self):
        """ 
        
        find injection threshold: Min IG1
        plot and model charge injection vs. IG1
        find notch injection amount
        
        
        """
        return

