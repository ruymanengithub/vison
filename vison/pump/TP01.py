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
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc

from copy import deepcopy
# END IMPORT 

isthere = os.path.exists

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
  


def build_TP01_scriptdict(Nshuffles_V,TOI_TPv,id_delays,
                    diffvalues=dict(),elvis='6.3.0'):
    """ """
    
    assert len(id_delays) == 2
    
    vpumpmodes = [123,234,341,412]
              
    TP01_sdict = dict()
    
    
    TP01_commvalues['ver_shuffles'] = Nshuffles_V
    
    # First Injection Drain Delay
    
    TP01_sdict['col1'] = dict(frames=1,v_tpump=0,comments='BGD',
              id_delay=id_delays[0])
    
    colcounter = 2
    for i,TOI_TP in enumerate(TOI_TPv):
        
        for k,vpumpmode in enumerate(vpumpmodes):
            colkey = 'col%i' % colcounter
            TP01_sdict[colkey] = dict(frames=1,toi_tp=TOI_TP,
                     id_dly=id_delays[0],v_tpmod=vpumpmode)
            
            colcounter += 1
    
    # Second Injection Drain Delay
    
    
    TP01_sdict['col%i' % colcounter] = dict(frames=1,v_tpump=0,comments='BGD',
              id_delay=id_delays[1])
    colcounter += 1

    for j,TOI_TP in enumerate(TOI_TPv):
        
        for k,vpumpmode in enumerate(vpumpmodes):
        
            colkey = 'col%i' % colcounter
            #print colkey
            TP01_sdict[colkey] = dict(frames=1,toi_tp=TOI_TP,
                         id_dly=id_delays[1],v_tpmod=vpumpmode)    
            
            colcounter += 1
    
    
    Ncols = len(TP01_sdict.keys())    
    TP01_sdict['Ncols'] = Ncols

    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(TP01_commvalues)    
   
    TP01_sdict = sc.update_structdict(TP01_sdict,commvalues,diffvalues)
    
    return TP01_sdict



def check_data(DataDict,report,inputs,log=None):
    """ 

    TP01: Checks quality of ingested data.
    

    **METACODE**
    
    ::

        check common HK values are within safe / nominal margins
        check voltages in HK match commanded voltages, within margins
    
        f.e.ObsID:
            f.e.CCD:
                f.e.Q.:
                    measure offsets in pre-, img-, over-
                    measure std in pre-, img-, over-
                    measure offsets/means in pre-, img-, over-
        
        assess std in pre- is within allocated margins
        assess offsets in pre-, and over- are equal, within allocated  margins
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
                    [subtract /? divide by reference image wo TPing]
                    save raw map of dipoles

    
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
                    load raw map of dipoles (from extract_data)
                    identify dipoles:
                        x, y, rel-amplitude, orientation

        produce & report:  
            map location of dipoles
            PDF of dipole amplitudes (for N and S)
            Counts of dipoles (and N vs. S)
    
    """
    
def meta_analysis():
    """
    
    Meta-analysis of data:
        
        Try to identify tau and pixel-phase location for each trap.
        Need to associate dipoles across TOI_TPs and TP-patterns
        
        
        NEED ADVICE FROM O.U. - ask Jesper, Ben, Julia


    **METACODE**
    
    ::
        
        across TOI_TP, patterns:
            build catalog of unique dipoles / traps
            get amplitude vs. TOI -> P(tau)
            location of trap within pixel I-PHASE
            
        DON'T KNOW HOW TO DEAL WITH DIFFERENT PATTERNS, NOR HOW TO 
        PRODUCE A USEFUL MAP OF TRAP LOCATIONS (X,Y,I-PHASE), P(TAU)
    

    """
    
