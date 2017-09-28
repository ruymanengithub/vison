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
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.datamodel import EXPLOGtools as ELtools
#from vison.datamodel import HKtools
#from vison.datamodel import ccd
#from vison.datamodel import generator
#import datetime
from copy import deepcopy
# END IMPORT

isthere = os.path.exists



CHINJ01_commvalues = dict(program='CALCAMP',test='CHINJ01',
  IG2_T=5500,IG2_B=5500,
  iphi1=1,iphi2=1,iphi3=1,iphi4=0,
  readmode_1='Normal',readmode_2='Normal',
  vertical_clk = 'Tri-level',serial_clk='Even mode',
  flushes=7,exptime=0.,shutter='Thorlabs SC10',
  electroshutter=0,vstart=1,vend=2066,
  sinvflush=0,chinj=1,chinj_rows_on=30,
  chinj_rows_off=100,chinj_repeat=15,id_width=60,chinj_ser_wait=1,
#  id_delay=100,
  tpump=0,ser_shuffles=1,
  ver_shuffles=1,dwell_v=0,dwell_h=0,motor=0,
  matrix_size=2,step_size=100,add_h_overscan=0,
  add_v_overscan=0,toi_flush=143.,toi_tpump=1000.,
  toi_rdout=1000.,toi_chinj=1000.,
  wavelength='Filter 4',pos_cal_mirror=polib.mirror_nom['Filter4'],
  operator='who',sn_ccd1='x',sn_ccd2='y',sn_ccd3='z',
  sn_roe='rr',sn_rpsu='pp',
  comments='')
  


def build_CHINJ01_scriptdict(IDL,IDH,IG1s,id_delays,toi_chinj,diffvalues=dict(),
                             elvis='6.0.0'):
    """
    Builds CHINJ01 script structure dictionary.
    
    :param IDL: int, [mV], value of IDL (Inject. Drain Low).
    :param IDH: int, [mV], Injection Drain High.
    :param IG1s: list of 2 ints, [mV], [min,max] values of IG1.
    :param id_delays: list of 2 ints, [mV], injection drain delays (2).
    :param toi_chinj: int, [us], TOI-charge injection.
    :param diffvalues: dict, opt, differential values.
    
    """
    
    CCDs = [1,2,3]
    halves = ['T','B']
    
    assert len(IG1s) == 2
    assert len(id_delays) == 2
    
    
    dIG1 = 0.25 * 1.E3 # Vx1E3
    NIG1 = (IG1s[1]-IG1s[0])/dIG1+1
    IG1v = np.arange(NIG1)*dIG1+IG1s[0]
    
    CHINJ01_sdict = dict()
    
    # First Injection Drain Delay
    
    colcounter = 1
    for i,IG1 in enumerate(IG1v):
        colkey = 'col%i' % (i+1,)
        #print colkey
        CHINJ01_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[0],toi_chinj=toi_chinj)
        
        for CCD in CCDs:
            for half in halves:
                CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD,half)] = IG1
        
        
        colcounter += 1
    
    # Second Injection Drain Delay
    
    colstart = colcounter

    for j,IG1 in enumerate(IG1v):
        colkey = 'col%i' % (colstart+j,)
        #print colkey
        CHINJ01_sdict[colkey] = dict(frames=1,IDL=IDL,IDH=IDH,
                     id_delay=id_delays[1],toi_chinj=toi_chinj)
        
        for CCD in CCDs:
            for half in halves:
                CHINJ01_sdict[colkey]['IG1_%i_%s' % (CCD,half)] = IG1

    
    Ncols = len(CHINJ01_sdict.keys())    
    CHINJ01_sdict['Ncols'] = Ncols
    
    commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
    commvalues.update(CHINJ01_commvalues)
    
    CHINJ01_sdict = sc.update_structdict(CHINJ01_sdict,commvalues,diffvalues)
    
    return CHINJ01_sdict


def check_data(DataDict,report,inputs,log=None):
    """ 

    CHINJ01: Checks quality of ingested data.
    

    **METACODE**
    
    ::

        check common HK values are within safe / nominal margins
        check voltages in HK match commanded voltages, within margins
    
        f.e.ObsID:
            f.e.CCD:
                f.e.Q.:
                    measure offsets in pre-, img-, over-
                    measure std in pre-, img-, over-
                    extract 2D chinj-pattern:
                        measure average level of injection
                        measure average level of non-injection
        
        assess std in pre- is within allocated margins
        assess offsets in pre-, and over- are equal, within allocated  margins
        assess offsets are within allocated margins
        assess non-injection level is within expected margins
        assess injection level is within expected margins
    
        [plot offsets vs. time]
        [plot std vs. time]
        plot injected level vs. IG1
    
    
        issue any warnings to log
        issue update to report          

    
    """
    

    
    return DataDict, report
    

def prep_data(DataDict,report,inputs,log=None):
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
    
    return DataDict,report
    
def basic_analysis(DataDict,report,inputs,log=None):
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
    
   
    
    return DataDict,report
    
#def meta_analysis(DataDict,report,inputs,log=None):
#    """ 
#    
#    CURRENTLY NOT NEEDED
#    
#    """
#    
#    return DataDict,report

