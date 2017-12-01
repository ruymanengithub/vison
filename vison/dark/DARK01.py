#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: DARK01

"Dark Current" analysis script

Created on Tue Aug 29 17:21:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
from copy import deepcopy

from vison.pipe import lib as pilib
from vison.point import  lib as polib
from vison.datamodel import scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
from vison.pipe.task import Task
from vison.image import performance
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']


DARK01_commvalues = dict(program='CALCAMP',test='DARK01',
  flushes=7,shuttr=1,
  siflsh=1,siflsh_p=500,
  comments='DARK')
  

class DARK01(Task):
    """ """

    def __init__(self,inputs,log=None,drill=False):
        """ """
        super(DARK01,self).__init__(inputs,log,drill)
        self.name = 'DARK01'
        self.HKKeys = HKKeys
        self.figdict = dict()
        
        self.perflimits.update(performance.perf_rdout)
    
    

    def build_scriptdict(self,N,exptime,diffvalues=dict(),elvis='6.3.0'):
        """Builds DARK01 script structure dictionary.
        
        :param N: integer, number of frames to acquire.
        :param exptime: integer, ms, exposure time.
        :param diffvalues: dict, opt, differential values.
        
        """
        
        DARK01_sdict = dict(col1=dict(frames=N,exptime=exptime))
        
        Ncols = len(DARK01_sdict.keys())    
        DARK01_sdict['Ncols'] = Ncols
        
                    
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(DARK01_commvalues)
        
        DARK01_sdict = sc.update_structdict(DARK01_sdict,commvalues,diffvalues)    
        
        return DARK01_sdict
    
    
    def check_data(self):
        """ 
        DARK0: Checks quality of ingested data.
        
        **METACODE**
        
        ::
        
            check common HK values are within safe / nominal margins
            check voltages in HK match commanded voltages, within margins
        
            f.e.ObsID:
                f.e.CCD:
                    f.e.Q.:
                        measure offsets/means in pre-, img-, over-
                        measure std in pre-, img-, over-
            assess std in pre- is within allocated margins
            assess offsets/means in pre-, img-, over- are equal, within allocated  margins
            assess offsets/means are within allocated margins
        
            plot offsets/means vs. time
            plot std vs. time
        
            issue any warnings to log
            issue update to report
        
        """
        
        
        raise NotImplementedError
        
    
    def prep_data(self):
        """
        
        DARK01: Preparation of data for further analysis.
        
        **METACODE**
        
        ::
            
            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        subtract offset: save to pick, update filename
    
        
        """
        
        
        raise NotImplementedError
        
        
    def basic_analysis(self):
        """ 
    
        DARK01: Basic analysis of data.
        
        **METACODE**
        
        ::
    
            f. e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        produce mask of hot pixels
                        count hot pixels / columns
                        produce a 2D poly model of masked-image, save coefficients
                        produce average profile along rows
                        produce average profile along cols
                        measure and save RON after subtracting large scale structure
                    save 2D model and profiles in a pick file for each OBSID-CCD
        
            plot average profiles f. each CCD and Q (color coded by time)
        
        """
        
        
        raise NotImplementedError
        
        
    def meta_analysis(self):
        """ 
        
        **METACODE**
        
        ::
    
            f. each CCD:
                f. e. Q:
                    stack all ObsIDs to produce Master Dark
                    produce mask of hot pixels / columns
                    count hot pixels / columns
                    measure average profile along rows
                    measure average profile along cols
            
            plot average profiles of Master Bias f. each CCD,Q
            show Master Dark (images), include in report
            report stats of defects, include in report
            save name of MasterDark to DataDict, report
            save name of Defects in Darkness Mask to DD, report
        
        
        """
        
        
        raise NotImplementedError
    
    
    def feeder(self,inputs,elvis='6.3.0'):
        """ """
        
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('basic',self.basic_analysis),
                    ('meta',self.meta_analysis)]
        
        N = inputs['N']
        exptime = inputs['exptime']
        if 'elvis' in inputs:
            elvis = inputs['elvis']
        if 'diffvalues' in inputs:
            diffvalues = inputs['diffvalues']
        else:
            diffvalues = {}
        
        scriptdict = self.build_scriptdict(N,exptime,diffvalues=diffvalues,elvis=elvis)
    
        inputs['structure'] = scriptdict
        inputs['subpaths'] = dict()
        
        if 'perflimits' in inputs:
            self.perflimits.update(inputs['perflimits'])
        
        
        return inputs
