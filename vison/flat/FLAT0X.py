#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FLAT0X

Flat-fields acquisition / analysis script

Created on Tue Aug 29 17:32:52 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
from copy import deepcopy
from collections import OrderedDict

from vison.support import context
from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
#from vison.pipe import FlatFielding as FFing
#from vison.support.report import Report
#from vison.support import files
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
#from vison.pipe.task import Task
from vison.flat.FlatTask import FlatTask
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']

FLAT0X_commvalues = dict(program='CALCAMP',
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,vstart=0,vend=2086,
  exptime=0.,shuttr=1,
  siflush=0,
  wave=4,
  comments='')

FLU_lims = dict(CCD1= dict(
                    col1=0.25 * 2**16 * (1.+np.array([-0.10,0.10])),
                    col2=0.50 * 2**16 * (1.+np.array([-0.10,0.10])),
                    col3=0.75 * 2**16 * (1.+np.array([-0.10,0.10]))))
for i in [2,3]: FLU_lims['CCD%i' % i] = deepcopy(FLU_lims['CCD1'])

class FLATS0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
            ('exptimes',([list],'Exposure times for each fluence.')),
            ('frames',([list],'Number of Frames for each fluence.')),
            ('wavelength',([int],'Wavelength')),
            ])))


class FLAT0X(FlatTask):
    """ """
    
    inputsclass = FLATS0X_inputs
    

    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(FLAT0X,self).__init__(inputs,log,drill,debug)
        self.name = 'FLAT0X'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('indivflats',self.do_indiv_flats),
                    ('masterflat',self.do_master_flat),
                    ('prmask',self.do_prdef_mask)]
        self.HKKeys = HKKeys
        self.figdict = dict() 
        self.inputs['subpaths'] = dict(figs='figs',pickles='ccdpickles')
        
   
    def set_inpdefaults(self,**kwargs):
        """ """
        
        try: wavelength = kwargs['wavelength']
        except KeyError: wavelength = 800
        try: test = kwargs['test']
        except KeyError: test = 'FLAT0X'
        
        t_dummy_F0X = np.array([25.,50.,75])/100.
        exptimesF0X = (ogse.tFWC_flat['nm%i' % wavelength] * t_dummy_F0X).tolist() # s
        framesF0X = [80,60,30]
        
        self.inpdefaults = dict(exptimes=exptimesF0X,
                       frames=framesF0X,
                       wavelength=wavelength,
                       test=test)
        
        
    def set_perfdefaults(self,**kwargs):
        #wavelength = self.inputs['wavelength']        
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)        
        self.perfdefaults['FLU_lims'] = FLU_lims
        
    

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """Builds FLAT0X script structure dictionary.
        
        :param diffvalues: dict, opt, differential values.
        
        """
        
        exptimes = self.inputs['exptimes']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        test = self.inputs['test']
        
        FW_ID = ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])
        
        FLAT0X_commvalues['wave'] = FW_IDX
        FLAT0X_commvalues['test'] = test
                         
        assert len(exptimes) == len(frames)
        #assert len(exptimes) == len(flags)
        
        FLAT0X_sdict = dict()
        for i,exptime in enumerate(exptimes):
            FLAT0X_sdict['col%i' % (i+1,)] = dict(frames=frames[i],exptime=exptimes[i],
                         comments='EXP%.1e' % exptime) #,comments=flags[i])
    
        Ncols = len(FLAT0X_sdict.keys())    
        FLAT0X_sdict['Ncols'] = Ncols
                    
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(FLAT0X_commvalues)
        
        FLAT0X_sdict = sc.update_structdict(FLAT0X_sdict,commvalues,diffvalues)
        
        return FLAT0X_sdict

    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """ """
        wavedkeys = []
        return super(FLAT0X,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys)
    
    
    
    def do_indiv_flats(self):
        """
        
        **METACODE**
        
        ::
        
            Preparation of data for further analysis and 
            produce flat-field for each OBSID.
    
            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        subtract offset
                        opt: [sub bias frame]
                        model 2D fluence distro in image area
                        produce average profile along rows
                        produce average profile along cols
                        
                    save 2D model and profiles in a pick file for each OBSID-CCD
                    divide by 2D model to produce indiv-flat
                    save indiv-Flat to FITS, update add filename
    
            plot average profiles f. each CCD and Q (color coded by time)
        
        """
        
        raise NotImplementedError
        
        
    def do_master_flat(self):
        """ 
        
        **METACODE**
        
        ::
        
            Produces Master Flat-Field
    
            f.e.CCD:
                f.e.Q:
                    stack individual flat-fields by chosen estimator
            save Master FF to FITS
            measure PRNU and 
            report PRNU figures
        
        """
        
        raise NotImplementedError
    
    
    def do_prdef_mask(self):
        """
        **METACODE**
        
        ::
        
            Produces mask of defects in Photo-Response
        
            f.e.CCD:
                f.e.Q:
                    produce mask of PR defects
                    save mask of PR defects
                    count dead pixels / columns 
        
            report PR-defects stats
        
        """    
        
        raise NotImplementedError
    
    
