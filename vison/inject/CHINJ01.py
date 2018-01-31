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
from copy import deepcopy
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
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

CHINJ01_commvalues = dict(program='CALCAMP',test='CHINJ01',
  IG2_T=5.5,IG2_B=5.5,
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,exptime=0.,vstart=0,vend=2086,
  shuttr=0,
  chinj=1,chinj_on=30,
  chinj_of=100,
  id_wid=60,chin_dly=1,
  operator='who',
  comments='')

class CHINJ01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('IDL',([float],'Injection Drain Low Voltage.')),
            ('IDH',([float],'Injection Drain High Voltage.')),
            ('IG1s',([list],'Injection Gate 1 Voltages.')),
            ('id_delays',([list],'Injection Drain Delays.')),
            ('toi_chinj',([int],'TOI Charge Injection.')),
            ])))


class CHINJ01(InjTask):
    """ """
    
    inputsclass = CHINJ01_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(CHINJ01,self).__init__(inputs,log,drill,debug)
        self.name = 'CHINJ01'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('extract',self.extract_data),
                         ('basic',self.basic_analysis)]
        self.HKKeys = HKKeys
        self.figdict = CH01aux.CH01figs.copy()
        self.inputs['subpaths'] = dict(figs='figs')


    def set_inpdefaults(self,**kwargs):
        """ """
        toi_chinj = 500
        
        self.inpdefaults = dict(
                IDL = 11.,
                IDH = 18.,
                IG1s = [2.,6.],
                id_delays = [toi_chinj*3.,toi_chinj*2.],
                toi_chinj = toi_chinj
                )
        
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        
        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()
        
        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()
        
        

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
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
        id_delays = self.inputs['id_delays']
        toi_chinj = self.inputs['toi_chinj']
        
        CCDs = [1,2,3]
        halves = ['T','B']
        
        assert len(IG1s) == 2
        assert len(id_delays) == 2
        
        
        dIG1 = 0.25  # V
        NIG1 = (IG1s[1]-IG1s[0])/dIG1+1
        IG1v = np.arange(NIG1)*dIG1+IG1s[0]
        
        CHINJ01_sdict = dict()
        
        # First Injection Drain Delay
        
        colcounter = 1
        for i,IG1 in enumerate(IG1v):
            colkey = 'col%i' % colcounter
            #print colkey
            CHINJ01_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                         id_dly=id_delays[0],toi_ch=toi_chinj)
            
            for CCD in CCDs:
                for half in halves:
                    CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD,half)] = IG1
            
            
            colcounter += 1
        
        # Second Injection Drain Delay
    
        for j,IG1 in enumerate(IG1v):
            colkey = 'col%i' % colcounter
            #print colkey
            CHINJ01_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                         id_dly=id_delays[1],toi_ch=toi_chinj)
            
            for CCD in CCDs:
                for half in halves:
                    CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD,half)] = IG1
    
            colcounter += 1
                    
                    
        Ncols = len(CHINJ01_sdict.keys())    
        CHINJ01_sdict['Ncols'] = Ncols
        
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(CHINJ01_commvalues)
        
        if len(diffvalues)==0:
            diffvalues = self.inputs['diffvalues']
        
        CHINJ01_sdict = sc.update_structdict(CHINJ01_sdict,commvalues,diffvalues)
        
        return CHINJ01_sdict
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(CHINJ01,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys)
    
    
    
    
    def extract_data(self):
        """
        
        **NEEDED?** Could be merged with basic_analysis
        
        **METACODE**
        
        ::
        
            Preparation of data for further analysis:
    
            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        subtract offset
                        extract average 2D injection pattern and save
    
        
        """
        
        raise NotImplementedError

        
    def basic_analysis(self):
        """ 
    
        Basic analysis of data.
    
        **METACODE**
        
        ::
    
            f. e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        load average 2D injection pattern
                        produce average profile along lines
                        measure charge-inj. non-uniformity
                        produce average profile across lines
                        measure charge spillover into non-injection
                        measure stats of injection (mean, med, std, min/max, percentiles)
                        
            plot average inj. profiles along lines f. each CCD, Q and IG1
                save as a rationalized set of curves
            plot average inj. profiles across lines f. each CCD, Q and IG1
                save as a rationalized set of  curves
            
            plot charge injection vs. IG1
            report injection stats as a table
        
        """
        
        raise NotImplementedError
        

