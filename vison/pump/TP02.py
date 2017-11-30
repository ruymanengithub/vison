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

from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
from vison.pipe.task import Task
# END IMPORT 

isthere = os.path.exists


IG1=6
IG2=5

TP02_commvalues = dict(program='CALCAMP',test='TP02',
  exptime=0.,shuttr=0,
  vstart=1,vend=100,
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
  

class TP02(Task):
    """ """

def build_TP02_scriptdict(Nshuffles_H,dwell_sv,id_delays,
                diffvalues=dict(),elvis='6.3.0'):
    """MISSING: different starting points (not implemented in ELVIS yet)."""
    
    assert len(id_delays) == 2
    
    sermodes = [23,31]
              
    TP02_sdict = dict()
    
    TP02_commvalues['ser_shuffles'] = Nshuffles_H
    
    # First Injection Drain Delay
    
    TP02_sdict['col1'] = dict(frames=1,v_tpump=0,s_tpump=0,
              comments='BGD',id_dly=id_delays[0])
    
    colcounter = 2
    for i,dwell_s in enumerate(dwell_sv):
        
        for k,sermode in enumerate(sermodes):
            colkey = 'col%i' % colcounter
            TP02_sdict[colkey] = dict(frames=1,dwell_s=dwell_s,
                     id_dly=id_delays[0],s_tpmod=sermode)
            
            colcounter += 1
    
    
    # Second Injection Drain Delay
    
    TP02_sdict['col%i' % colcounter] = dict(frames=1,v_tpump=0,s_tpump=0,
              comments='BGD',id_dly=id_delays[1])
    colcounter += 1 
    

    for j,dwell_s in enumerate(dwell_sv):
        
        for k,sermode in enumerate(sermodes):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP02_sdict[colkey] = dict(frames=1,dwell_s=dwell_s,
                         id_dly=id_delays[1],s_tpmod=sermode)    
            
            colcounter += 1
            
    
    Ncols = len(TP02_sdict.keys())    
    TP02_sdict['Ncols'] = Ncols
              
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(TP02_commvalues)
    
    TP02_sdict = sc.update_structdict(TP02_sdict,commvalues,diffvalues)
    
    return TP02_sdict


def check_data(DataDict,report,inputs,log=None):
    """ 

    TP02: Checks quality of ingested data.
    

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

def prep_data(DataDict,report,inputs,log=None):
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
    
    return DataDict,report

def basic_analysis():
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
    
def meta_analysis():
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
    
