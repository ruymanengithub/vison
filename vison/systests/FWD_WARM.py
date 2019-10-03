#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: FPA-FWD (WARM)

Created on Thu Sep 26 16:11:14 2019

@author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import sys

from vison.image import bits
from vison.other import MOT_FFaux
from vison.systests import fpatask
from vison.datamodel import inputs

# END IMPORT

class FWD_WARM_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()

class FWD_WARM(fpatask.FpaTask):
    
    inputsclass = FWD_WARM_inputs
    
    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        
        self.subtasks = [('check', self.check_data), 
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        
        super(FWD_WARM, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                        cleanafter=cleanafter)
        
        self.name = 'FWD_WARM'
        self.type = 'Simple'
        
        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')
    

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(preprocessing=dict(
                                            offsetkwargs=dict(
                                            ignore_pover= True, 
                                            trimscan = [25, 5], 
                                            method = 'median',
                                            extension= -1, 
                                            scan = 'pre' 
                                                )))    
    
    def basic_analysis(self):
        """        
        To-Do:
            Extract average profiles along columns (image area)
            Extract HER profiles and metrics
            Extract bits histograms
        """
        
        if self.report is not None:
            self.report.add_Section(keyword='extract', Title='Extraction', level=0)
        
        return # TESTS
        
        iObs = 0
        
        LE1fits = self.dd.mx['File_name'][iObs]

        fullLE1fits = os.path.join(self.datapath, LE1fits)
        
        LE1 = self.load_LE1(fullLE1fits)
        vstart = 0
        vend = 2086
                
        for jY in range(1,LE1.fpamodel.NSLICES+1):
            for iX in range(1,LE1.fpamodel.NSLICES+1):
                CCDID = 'C_%i%i' % (jY,iX)
                
                kccdobj = LE1.get_ccdobj(CCDID)
                
                for Q in LE1.Quads:
                    
                    # vertical profile: RAMP
                    
                    vQ = kccdobj.get_1Dprofile(Q=Q, orient='ver', area='img', stacker='median',
                                                         vstart=vstart, 
                                                         vend=vend)
                    
                    # save vQ: PENDING
                    
                    # HER
                    
                    HERprof = MOT_FFaux.extract_overscan_profiles(kccdobj, 
                                                    [1.E3, 4.E4], 
                                                    direction='serial')
                    
                    # save HERprof: PENDING
                    # save HER value: PENDING
                    
                    # RAMP; bit-histograms extraction
                    
                    bitsmean = bits.get_histo_bits(kccdobj, Q, vstart=vstart, vend=vend)
                    bitsbin = np.arange(0,16)+0.5
                    bitshisto = dict(bins=bitsbin, H=bitsmean)
                    
                    # save bit histograms: PENDING
                    
                    
        
    def meta_analysis(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)
        
        return # TESTS
    