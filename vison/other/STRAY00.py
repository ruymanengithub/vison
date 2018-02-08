#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: STRAY00 - used to investigate STRAY-LIGHT sources in OGSE.
  NOT intended for performance evaluation.
  COMMISSIONING.


Created on Thu Feb 08 14:07:00 2018

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
#from InjTask import InjTask
from vison.image import performance
from vison.datamodel import inputs
from vison.dark.DarkTask import DarkTask
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

STRAY00_commvalues = dict(program='CALCAMP',test='STRAY00',
  flushes=7,exptime=0.,shuttr=0,
  e_shuttr=0,vstart=0,vend=2086,
  siflush=0,#sinvflushp=500,
  chinj=0,
  s_tpump=0,
  v_tpump=0,
  motr_on=0,
  toi_fl=143.,toi_tp=1000.,toi_ro=1000.,toi_ch=1000.,
  wave=4,
  comments='')
  

class STRAY00_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    #manifesto.update(OrderedDict(sorted([
    #        ('IDH',([float],'Injection Drain High Voltage.')),
    #        ('toi_chinj',([int],'TOI Charge Injection.')),
    #        ('chinj_on',([int],'Number of lines injected per cycle.')),
    #        ('chinj_of',([int],'Number of lines NON injected per cycle.'))
    #        ])))

class STRAY00(DarkTask):
    """ """
    
    inputsclass = STRAY00_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(STRAY00,self).__init__(inputs,log,drill,debug)
        self.name = 'STRAY00'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data)]
        self.HKKeys = HKKeys
        self.figdict = dict()
        self.inputs['subpaths'] = dict(figs='figs')

    def set_inpdefaults(self,**kwargs):
        """ """        
        self.inpdefaults = dict()
        
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """
        Builds STRAY00 script structure dictionary.        
        :param diffvalues: dict, opt, differential values.
        
        """
        

        STRAY00_sdict = dict()
        
        # start with all lights on in LAB
        STRAY00_sdict['col1'] = dict(frames=1,exptime=0,shuttr=0,wave=1,
                     source='flat',
                     comment = 'LABLIT')
        STRAY00_sdict['col2'] = dict(frames=1,exptime=100,shuttr=0,wave=1,
                     source='flat',
                     comment='LABLIT')
        # switch off all lights in LAB
        STRAY00_sdict['col3'] = dict(frames=1,exptime=0,shuttr=0,wave=1,
                     source='flat',
                     comment='LABOUTW1')
        STRAY00_sdict['col4'] = dict(frames=1,exptime=100,shuttr=0,wave=1,
                     source='flat',
                     comment='LABOUTW1')
        # change light source wavelength
        STRAY00_sdict['col5'] = dict(frames=1,exptime=0,shuttr=0,wave=4,
                     source='flat',
                     comment='LABOUTW4')
        STRAY00_sdict['col6'] = dict(frames=1,exptime=100,shuttr=0,wave=4,
                     source='flat',
                     comment='LABOUTW4')                
        # switch off Light Source
        STRAY00_sdict['col7'] = dict(frames=1,exptime=0,shuttr=0,wave=1,
                     comment='TUNGSOFF')
        STRAY00_sdict['col8'] = dict(frames=1,exptime=100,shuttr=0,wave=1,
                     source='flat',
                     comment='TUNGSOFF')
        # switch off Pressure Gauge - Light Source still Off
        STRAY00_sdict['col9'] = dict(frames=1,exptime=0,shuttr=0,wave=1,
                     comment='GAUGEOFF')
        STRAY00_sdict['col10'] = dict(frames=1,exptime=100,shuttr=0,wave=1,
                     source='flat',
                     comment='GAUGEOFF')
        # Switch on Pressure Gauge and Light Source
        STRAY00_sdict['col11'] = dict(frames=1,exptime=0,shuttr=0,wave=1,
                     source='flat',
                     comment='TUNGGAUON')
        STRAY00_sdict['col12'] = dict(frames=1,exptime=100,shuttr=0,wave=1,
                     source='flat',
                     comment='TUNGGAUON')
        
        # NON-supervised, Lights Off, Gauge On, Tungsten On
        
        # DARK FF-source wave = 1, 4

        STRAY00_sdict['col13'] = dict(frames=1,exptime=300,shuttr=0,wave=1,
                     source='flat',
                     comment='DARK-FF')
        
        STRAY00_sdict['col14'] = dict(frames=1,exptime=300,shuttr=0,wave=4,
                     source='flat',
                     comment='DARK-FF')        
        
        # DARK PSF-source wave = 1, mirror at near end

        STRAY00_sdict['col15'] = dict(frames=1,exptime=300,shuttr=0,wave=1,
                     source='point',
                     mirr_on=1,mirr_pos=1.,
                     comment='DARK-PSF')
        
        # DARK PSF-source wave = 1, mirror at far end
        
        STRAY00_sdict['col16'] = dict(frames=1,exptime=300,shuttr=0,wave=1,
                     source='point',
                     mirr_on=1,mirr_pos=99.,
                     comment='DARK-PSF')        
        
        # DARK PSF-source wave = 4, mirror at near end
        
        STRAY00_sdict['col17'] = dict(frames=1,exptime=300,shuttr=0,wave=4,
                     source='point',
                     mirr_on=1,mirr_pos=1.,
                     comment='DARK-PSF')
        
        # DARK PSF-source wave = 4, mirror at far end
        
        STRAY00_sdict['col18'] = dict(frames=1,exptime=300,shuttr=0,wave=4,
                     source='point',
                     mirr_on=1,mirr_pos=99.,
                     comment='DARK-PSF')        
        
        
        # PSF - Short, wave = 1
        
        STRAY00_sdict['col19'] = dict(frames=1,exptime=5.,shuttr=1,wave=1,
                     source='point',
                     mirr_on=1,mirr_pos=50.,
                     comment='PSF-S')
        
        STRAY00_sdict['col20'] = dict(frames=1,exptime=50.,shuttr=1,wave=1,
                     source='point',
                     mirr_on=1,mirr_pos=50.,
                     comment='PSF-L')
        
        # PSF - Long, wave = 4
        
        STRAY00_sdict['col21'] = dict(frames=1,exptime=5.,shuttr=1,wave=4,
                     source='point',
                     mirr_on=1,mirr_pos=50.,
                     comment='PSF-S')
        
        STRAY00_sdict['col22'] = dict(frames=1,exptime=50.,shuttr=1,wave=4,
                     source='point',
                     mirr_on=1,mirr_pos=50.,
                     comment='PSF-L')        
        
        
        
        Ncols = len(STRAY00_sdict.keys())    
        STRAY00_sdict['Ncols'] = Ncols
        
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(STRAY00_commvalues)
        

        if len(diffvalues)==0:
            try: diffvalues = self.inputs['diffvalues']
            except: diffvalues = diffvalues = dict()
        
        STRAY00_sdict = sc.update_structdict(STRAY00_sdict,commvalues,diffvalues)
        
        return STRAY00_sdict
    
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = []
        return super(STRAY00,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys)
    
    
    
