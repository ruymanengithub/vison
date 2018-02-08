#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: CHINJ02

Charge injection calibration (part 2)
    Injection vs. IDL (injection threshold)

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
#from vison.pipe import lib as pilib
#from vison.point import lib as polib
#from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#import datetime
#from vison.pipe.task import Task
from vison.datamodel import inputs
from InjTask import InjTask
from vison.image import performance
import CH02aux
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

IG1comm = 6.
IG2comm = 4.

CHINJ02_commvalues = dict(program='CALCAMP',test='CHINJ02',
  IG1_1_T=IG1comm,IG1_2_T=IG1comm,IG1_3_T=IG1comm,
  IG1_1_B=IG1comm,IG1_2_B=IG1comm,IG1_3_B=IG1comm,
  IG2_T=IG2comm,IG2_B=IG2comm,
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,vstart=0,vend=2086,
  exptime=0.,shuttr=0,
  siflsh=0,
  chinj=1,chinj_on=30,
  chinj_of=100,
  id_wid=60,
  chin_dly=1,
  comments='')

class CHINJ02_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('IDLs',([list],'Injection Drain Low Voltages List.')),
            ('IDH',([float],'Injection Drain High Voltage.')), 
            ('id_delays',([list],'Injection Drain Delays.')),
            ('toi_chinj',([int],'TOI Charge Injection.')),
            ])))

class CHINJ02(InjTask):
    """ """
    
    inputsclass = CHINJ02_inputs

    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(CHINJ02,self).__init__(inputs,log,drill,debug)
        self.name = 'CHINJ02'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('extract',self.extract_data),
                         ('basic',self.basic_analysis),('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = CH02aux.CH02figs.copy()
        self.inputs['subpaths'] = dict(figs='figs')
        
    
    def set_inpdefaults(self,**kwargs):
        """ """
        toi_chinj = 500
        
        self.inpdefaults = dict(
        IDLs = [10.,13.],
        IDH = 18.,
        id_delays = [toi_chinj*2.5,toi_chinj*1.5],
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
        Builds CHINJ02 script structure dictionary.
        
        #:param IDLs: list of 2 ints, [mV], [min,max] values of IDL (Inject. Drain Low).
        #:param IDH: int, [mV], Injection Drain High.
        #:param id_delays: list of 2 ints, [mV], injection drain delays (2).
        #:param toi_chinj: int, [us], TOI-charge injection.
        :param diffvalues: dict, opt, differential values.
        
        """
        
        IDLs = self.inputs['IDLs']
        IDH = self.inputs['IDH']
        id_delays = self.inputs['id_delays']
        toi_chinj = self.inputs['toi_chinj']
        
        assert len(IDLs) == 2
        assert len(id_delays) == 2
        
        dIDL = 0.25 # V
        NIDL = (IDLs[1]-IDLs[0])/dIDL+1
        IDLv = np.arange(NIDL)*dIDL+IDLs[0]
        
        CHINJ02_sdict = dict()
        
        
        # First Injection Drain Delay
        
        colcounter = 1
        for i,IDL in enumerate(IDLv):
            colkey = 'col%i' % (i+1,)
            CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                         id_dly=id_delays[0],toi_ch=toi_chinj)
            colcounter += 1
        
        # Second Injection Drain Delay
        
        colstart  = colcounter
    
        for j,IDL in enumerate(IDLv):
            colkey = 'col%i' % (colstart+j,)
            CHINJ02_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                         id_dly=id_delays[1],toi_ch=toi_chinj)
        
        Ncols = len(CHINJ02_sdict.keys())    
        CHINJ02_sdict['Ncols'] = Ncols
        
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(CHINJ02_commvalues)
        
        if len(diffvalues)==0:
            try: diffvalues = self.inputs['diffvalues']
            except: diffvalues = diffvalues = dict()
        
        CHINJ02_sdict = sc.update_structdict(CHINJ02_sdict,commvalues,diffvalues)
        
        return CHINJ02_sdict

    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(CHINJ02,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
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
        AS IT IS, REPEATS WHAT'S DONE IN THE CHECK_DATA. CONSIDER MERGING/SKIPPING
    
        **METACODE**
        
        ::
    
            f. e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        load average 2D injection pattern
                        produce average profile along lines
                        [measure charge-inj. non-uniformity]
                        [produce average profile across lines]
                        [measure charge spillover into non-injection]
                        measure stats of injection (mean, med, std, min/max, percentiles)
                        
            [plot average inj. profiles along lines f. each CCD, Q and IG1]
            [    save as a rationalized set of curves]
            [plot average inj. profiles across lines f. each CCD, Q and IG1]
            [    save as a rationalized set of  curves]
            
            save&plot charge injection vs. IDL
            report injection stats as a table
        
        """
        
       
        
        raise NotImplementedError
    
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
        
        raise NotImplementedError
    
