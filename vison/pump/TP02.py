#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (serial)

Created on Tue Aug 29 17:38:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from copy import deepcopy
from collections import OrderedDict

from vison.pipe import lib as pilib
from vison.support import context
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.pipe.task import Task
from PumpTask import PumpTask
from vison.image import performance
from vison.datamodel import inputs
import TP02aux
# END IMPORT 

isthere = os.path.exists

HKKeys = []

IG1=6
IG2=5

TP02_commvalues = dict(program='CALCAMP',test='TP02',
  exptime=0.,shuttr=0,
  vstart=0,vend=100,
  siflsh=1,siflsh_p=500,
  IDL=11,IDH=18,
  IG1_1_T=IG1,IG1_2_T=IG1,IG1_3_T=IG1,
  IG1_1_B=IG1,IG1_2_B=IG1,IG1_3_B=IG1,
  IG2_T=IG2,IG2_B=IG2,
  chinj=1,chinj_on=2066,chinj_of=0,
  chin_dly=0,
  s_tpump=1,s_tp_cnt=5000,
  v_tp_cnt=0,dwell_v=0,dwell_s=0,
  comments='')
  

class TP02_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('toi_chinj',([int],'TOI Charge Injection.')),            
            ('Nshuffles_H',([int],'Number of Shuffles, Horizontal/Serial Pumping.')),
            ('dwell_sv',([list],'Dwell Times list [serial].')),
            ('id_delays',([list],'Injection Drain Delays [2, one per CCDs section].')),
            ('spumpmodes',([list],'Horizontal/Serial Pumping Starting points.'))
            ])))

class TP02(PumpTask):
    """ """
    
    inputsclass = TP02_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(TP02,self).__init__(inputs,log,drill,debug)
        self.name = 'TP02'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = TP02aux.TP02figs.copy()
        self.inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')
        
    def set_inpdefaults(self,**kwargs):
        """ """
        toi_chinj=250
        self.inpdefaults = dict(toi_chinj=toi_chinj,
                         Nshuffles_H=5000,
                         dwell_sv=[0.,4.75,14.3,28.6],
                         id_delays=np.array([3.,2.])*toi_chinj,
                         spumpmodes=[23,31])
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)
        
        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()
        
        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()


    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """ """
        
        Nshuffles_H = self.inputs['Nshuffles_H']
        dwell_sv = self.inputs['dwell_sv']
        id_delays = self.inputs['id_delays']
        spumpmodes = self.inputs['spumpmodes']
        
        assert len(id_delays) == 2
        
        TP02_sdict = dict()
        
        TP02_commvalues['ser_shuffles'] = Nshuffles_H
        
        # First Injection Drain Delay
        
        TP02_sdict['col1'] = dict(frames=1,v_tpump=0,s_tpump=0,
                  comments='BGD',id_dly=id_delays[0])
        
        colcounter = 2
        for i,dwell_s in enumerate(dwell_sv):
            
            for k,sermode in enumerate(spumpmodes):
                colkey = 'col%i' % colcounter
                TP02_sdict[colkey] = dict(frames=1,dwell_s=dwell_s,
                         id_dly=id_delays[0],s_tpmod=sermode)
                
                colcounter += 1
        
        # Second Injection Drain Delay
        
        TP02_sdict['col%i' % colcounter] = dict(frames=1,v_tpump=0,s_tpump=0,
                  comments='BGD',id_dly=id_delays[1])
        colcounter += 1 
        
        for j,dwell_s in enumerate(dwell_sv):
            
            for k,sermode in enumerate(spumpmodes):
            
                colkey = 'col%i' % colcounter
                #print colkey
                TP02_sdict[colkey] = dict(frames=1,dwell_s=dwell_s,
                             id_dly=id_delays[1],s_tpmod=sermode)    
                
                colcounter += 1
                
        
        Ncols = len(TP02_sdict.keys())    
        TP02_sdict['Ncols'] = Ncols
                  
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP02_commvalues)
        
        if len(diffvalues)==0:
            diffvalues = self.inputs['diffvalues']
        
        TP02_sdict = sc.update_structdict(TP02_sdict,commvalues,diffvalues)
        
        return TP02_sdict
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(TP02,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys)

    
    
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
                        average across readout lines (iterations)
                        save raw 1D map of relative pumping
    
        
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
                        load raw 1D map of relative pumping (from extract_data)
                        identify dipoles:
                            x, rel-amplitude, orientation (E or W)
    
            produce & report:  
                map location of dipoles
                PDF of dipole amplitudes (for E and W)
                Counts of dipoles (and E vs. W)
        
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
                build catalog of traps: x,y,R-phase, amp(dwell)
                from Amp(dwell) -> tau, Pc
                
            Report on :
               Histogram of Taus
               Histogram of Pc (capture probability)
               Histogram of R-phases
    
               Total Count of Traps
    
        
    
        """
        raise NotImplementedError
        
        
    
    
