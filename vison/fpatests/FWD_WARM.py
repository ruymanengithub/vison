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
from vison.fpatests import fpatask
from vison.datamodel import inputs

# END IMPORT

#def _basic_onLE1(mode,kwargs_retrieve,kwargs_apply):
#    print(mode)
#    if mode == 'retrieve':
#        
#        queue = kwargs_retrieve['queue']
#        CCDID = kwargs_retrieve['CCDID']
#        LE1 = kwargs_retrieve['LE1']
#        print('processing %s' % CCDID)
#        kccdobj = LE1.get_ccdobj(CCDID)
#        
#        reply = dict(CCDID=CCDID,
#                     Quads=kccdobj.Quads)
#        for Q in kccdobj.Quads:
#            reply[Q] = kccdobj.get_stats(Q,sector='img',statkeys=['mean'])
#        
#        queue.put(reply)
#
#    elif mode == 'apply':
#        
#        Quads = kwargs_apply['Quads']
#        CCDID = kwargs_apply['CCDID']
#        holder = kwargs_apply['holder']
#        print('applying %s' % CCDID)
#        for Q in Quads:
#            holder[CCDID][Q] = reply[Q]

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
    
    def _basic_onLE1(self,**kwargs):
        
        
        CCDID = kwargs['CCDID']
        LE1 = kwargs['LE1']
        vstart = kwargs['vstart']
        vend = kwargs['vend']
        kccdobj = LE1.get_ccdobj(CCDID)
                
        for Q in LE1.Quads:
            
            # vertical profile: RAMP
            
            vQ = kccdobj.get_1Dprofile(Q=Q, orient='ver', area='img', stacker='median',
                                     vstart=vstart, vend=vend)
            
            
            # save vQ: PENDING
            
            # extract and save slope of profile for comparison: PENDING
            
            # HER
            
            HERprofQ = MOT_FFaux.extract_overscan_profiles(kccdobj, 
                                            [1.E3, 4.E4], 
                                            direction='serial')
            
            # save HERprof: PENDING
            # save HER value: PENDING
            
            # RAMP; bit-histograms extraction
            
            bitsmean = bits.get_histo_bits(kccdobj, Q, vstart=vstart, vend=vend)
            bitsbin = np.arange(0,16)+0.5
            bitshistoQ = dict(bins=bitsbin, H=bitsmean)
            
            # save bit histograms: PENDING
            
            stop()
    
    
    def basic_analysis(self):
        """        
        To-Do:
            Extract average profiles along columns (image area)
            Extract HER profiles and metrics
            Extract bits histograms
        """
        
        if self.report is not None:
            self.report.add_Section(keyword='extract', Title='Extraction', level=0)
                
        iObs = 0
        vstart = 0
        vend = 2086
        
        LE1fits = self.dd.mx['File_name'][iObs]

        fullLE1fits = os.path.join(self.inputs['datapath'], LE1fits)
        
        LE1 = self.load_LE1(fullLE1fits)
        
        
        kwargs = dict(vstart=vstart,
                      vend=vend)
        
        
        self.iterate_over_CCDs(LE1, FWD_WARM._basic_onLE1, **kwargs)
        
        
                    
        
    def meta_analysis(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)
        
        return # TESTS
    
