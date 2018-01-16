#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (vertical)

Created on Tue Aug 29 17:37:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from copy import deepcopy

from vison.pipe import lib as pilib
from vison.support import context
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
from vison.pipe.task import Task
from PumpTask import PumpTask
from vison.image import performance
# END IMPORT 

isthere = os.path.exists

HKKeys = []

IG1=6
IG2=5

TP01_commvalues = dict(program='CALCAMP',test='TP01',
  flushes=7,exptime=0.,shuttr=0,
  e_shuttr=0,vstart=1,vend=2066,
  siflsh=1,siflsh_p=500,
  IDL=11,IDH=18,
  IG1_1_T=IG1,IG1_2_T=IG1,IG1_3_T=IG1,
  IG1_1_B=IG1,IG1_2_B=IG1,IG1_3_B=IG1,
  IG2_T=IG2,IG2_B=IG2,
  chinj=1,chinj_on=2066,chinj_of=0,
  chin_dly=0,
  v_tpump=1,s_tpump=0,
  v_tp_cnt=1000,
  dwell_v=0,dwell_s=0,
  motr_on=0,
  comments='')
  

class TP01(PumpTask):
    """ """

    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(TP01,self).__init__(inputs,log,drill,debug)
        self.name = 'TP01'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = dict()
        self.inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')


    def set_inpdefaults(self,**kwargs):
        """ """
        
        toi_chinjTP01 = 250
        self.inpdefaults = dict(toi_chinj=toi_chinjTP01,
                         Nshuffles_V=5000,
                         id_delays=np.array([3.,2.]) * toi_chinjTP01,
                         toi_tpv=[200,1000,2000,4000,8000],
                         vpumpmodes = [123,234,341,412])
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """ """
        
        Nshuffles_V = self.inputs['Nshuffles_V']
        toi_tpv = self.inputs['toi_tpv']
        id_delays = self.inputs['id_delays']
        vpumpmodes = self.inputs['vpumpmodes']
        
        assert len(id_delays) == 2
        
        TP01_sdict = dict()
        
        TP01_commvalues['ver_shuffles'] = Nshuffles_V
        
        # First Injection Drain Delay
        
        TP01_sdict['col1'] = dict(frames=1,v_tpump=0,comments='BGD',
                  id_delay=id_delays[0])
        
        colcounter = 2
        for i,toi_tp in enumerate(toi_tpv):
            
            for k,vpumpmode in enumerate(vpumpmodes):
                colkey = 'col%i' % colcounter
                TP01_sdict[colkey] = dict(frames=1,toi_tp=toi_tp,
                         id_dly=id_delays[0],v_tpmod=vpumpmode)
                
                colcounter += 1
        
        # Second Injection Drain Delay
        
        
        TP01_sdict['col%i' % colcounter] = dict(frames=1,v_tpump=0,comments='BGD',
                  id_delay=id_delays[1])
        colcounter += 1
    
        for j,toi_tp in enumerate(toi_tpv):
            
            for k,vpumpmode in enumerate(vpumpmodes):
            
                colkey = 'col%i' % colcounter
                #print colkey
                TP01_sdict[colkey] = dict(frames=1,toi_tp=toi_tp,
                             id_dly=id_delays[1],v_tpmod=vpumpmode)    
                
                colcounter += 1
        
        
        Ncols = len(TP01_sdict.keys())    
        TP01_sdict['Ncols'] = Ncols
    
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP01_commvalues)    
       
        TP01_sdict = sc.update_structdict(TP01_sdict,commvalues,diffvalues)
        
        return TP01_sdict

    def filterexposures(self,structure,explogf,datapath,OBSID_lims,elvis=context.elvis):
        """ """
        wavedkeys = []
        return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys,elvis=elvis)
        
    
    def check_data(self):
        """ 
    
        TP01: Checks quality of ingested data.
        
    
        **METACODE**
        
        ::
    
            check common HK values are within safe / nominal margins
            check voltages in HK match commanded voltages, within margins
        
            f.e.ObsID:
                f.e.CCD:
                    f.e.Q.:
                        measure offsets in pre-, over-
                        measure std in pre-, over-
                        measure mean in img-
            
            assess std in pre- (~RON) is within allocated margins
            assess offsets in pre-, and over- are equal, within allocated margins
            assess offsets are within allocated margins
            assess injection level is within expected margins
        
            plot histogram of injected levels for each Q
            [plot std vs. time]
        
            issue any warnings to log
            issue update to report          
    
        
        """
        raise NotImplementedError
    
    def prep_data(self):
        """
        
        **METACODE**
        
        ::
        
            Preparation of data for further analysis:
    
            f.e. ObsID [images with TPing only]:
                f.e.CCD:
                    f.e.Q:
                        subtract offset
                        divide by reference image wo TPing
                        save "map of relative pumping"
    
        
        """
        
        raise NotImplementedError
        
    
    def basic_analysis(self):
        """
        
        Basic analysis of data.
    
        **METACODE**
        
        ::
    
            f. e. ObsID [there are different TOI_TP and TP-patterns]:
                f.e.CCD:
                    f.e.Q:
                        load "map of relative pumping"
                        find_dipoles:
                            x, y, rel-amplitude, orientation
    
            produce & report:  
                map location of dipoles
                PDF of dipole amplitudes (for N and S)
                Counts of dipoles (and N vs. S)
        
        """
        raise NotImplementedError
        
    def meta_analysis(self):
        """
        
        Meta-analysis of data:
            
            Try to identify tau and pixel-phase location for each trap.
            Need to associate dipoles across TOI_TPs and TP-patterns
    
    
        **METACODE**
        
        ::
            
            across TOI_TP, patterns:
                build catalog of traps: x,y,I-phase, Amp
                from Amp(TOI) -> tau, Pc
            
            Report on :
                Histogram of Taus
                Histogram of Pc (capture probability)
                Histogram of I-phases (larger phases should have more traps, 
                                  statistically) -> check
    
                Total Count of Traps
    
        """
        raise NotImplementedError
        
    
