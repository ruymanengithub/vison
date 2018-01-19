# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: NL01

End-To-End Non-Linearity Curve


Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract offset level.
    - Divide by Flat-field.
    - Synoptic analysis:
        fluence ratios vs. extime ratios >> non-linearity curve
    - extract: Non-Linearity curve for each CCD and quadrant
    - produce synoptic figures
    - Save results.


Created on Mon Apr  3 17:38:00 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from copy import deepcopy
from collections import OrderedDict

from vison.support import context
from vison.pipe import lib as pilib
from vison.ogse import ogse
from vison.point import lib as polib
from vison.datamodel import  scriptic as sc
import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
#from vison.pipe.task import Task
from FlatTask import FlatTask
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T','CCD2_OD_T','CCD3_OD_T','COMM_RD_T',
'CCD2_IG1_T','CCD3_IG1_T','CCD1_TEMP_T','CCD2_TEMP_T','CCD3_TEMP_T',
'CCD1_IG1_T','COMM_IG2_T','FPGA_PCB_TEMP_T','CCD1_OD_B',
'CCD2_OD_B','CCD3_OD_B','COMM_RD_B','CCD2_IG1_B','CCD3_IG1_B','CCD1_TEMP_B',
'CCD2_TEMP_B','CCD3_TEMP_B','CCD1_IG1_B','COMM_IG2_B']


NL01_commvalues = dict(program='CALCAMP',
  test='NL01',
  IPHI1=1,IPHI2=1,IPHI3=1,IPHI4=0,
  rdmode='fwd_bas',
  flushes=7,exptime=0.,shuttr=1,
  siflsh=1,siflsh_p=500,
  wave=6,
  source='flat',
  comments='')

class NL01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs
    manifesto.update(OrderedDict(sorted([
            ('exptimes',([dict,list],'Exposure times for each fluence.')),
            ('exptinter',([float],'Exposure time for interleaved fluence stability control frames.')),
            ('frames',([list],'Number of Frames for each fluence.')),
            ('wavelength',([int],'Wavelength'))
            ])))



class NL01(FlatTask):
    """ """
    
    inputsclass = NL01_inputs
    
    def __init__(self,inputs,log=None,drill=False,debug=False):
        """ """
        super(NL01,self).__init__(inputs,log,drill,debug)
        self.name = 'NL01'
        self.type = 'Simple'
        self.subtasks = [('check',self.check_data),('prep',self.prep_data),
                    ('extract',self.extract_stats),
                    ('NL',self.produce_NLCs),
                    ('satCTE',self.do_satCTE)]
        self.HKKeys = HKKeys
        self.figdict = dict() # B01aux.B01figs
        self.inputs['subpaths'] = dict() # dict(figs='figs',pickles='ccdpickles')
        
        
    def set_inpdefaults(self,**kwargs):

        expts = (np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100. * ogse.tFWC_flat['nm0']).tolist() # ms
        self.inpdefaults = dict(exptimes=expts,
                       exptinter=0.5 * ogse.tFWC_flat['nm0'],
                       frames=(np.ones(11,dtype='int32')*5).tolist(),           
                       wavelength=0,
                       )
        
    def set_perfdefaults(self,**kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)

    def build_scriptdict(self,diffvalues=dict(),elvis=context.elvis):
        """Builds NL01 script structure dictionary.
        
        #:param expts: list of ints [ms], exposure times.
        #:param exptinter: int, ms, exposure time of interleaved source-stability exposures.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 0 (Neutral Density Filter)
        :param diffvalues: dict, opt, differential values.
        """
        
        expts = self.inputs['exptimes']
        exptinter = self.inputs['exptinter']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        
        assert  len(expts) == len(frames)    
        
        FW_ID = ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])
        
        NL01_commvalues['wave'] = FW_IDX
        
        NL01_sdict = dict()
        
        NL01_sdict['col1'] = dict(frames=4,exptime=0,comment='BGD')
        NL01_sdict['col2'] = dict(frames=1,exptime=exptinter,comment='STAB')
        
        for ix,ifra in enumerate(frames):
    
            iexp = expts[ix]
            
            colkeyFlu = 'col%i' % (ix*2+2,)
        
            NL01_sdict[colkeyFlu] = dict(frames=ifra,exptime=iexp,comment='Fluence%i' % (ix+1,))
            
            colkeySta = 'col%i' % (ix*2+2+1,)
        
            NL01_sdict[colkeySta] = dict(frames=1,exptime=exptinter,comment='STA')
        
    
        Ncols = len(NL01_sdict.keys())    
        NL01_sdict['Ncols'] = Ncols
                  
        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(NL01_commvalues)
                
        NL01_sdict = sc.update_structdict(NL01_sdict,commvalues,diffvalues)
        
        return NL01_sdict
    
    
    def filterexposures(self,structure,explogf,datapath,OBSID_lims):
        """Loads a list of Exposure Logs and selects exposures from test PSF0X.
        
        The filtering takes into account an expected structure for the 
        acquisition script.
    
        The datapath becomes another column in DataDict. This helps dealing
        with tests that run overnight and for which the input data is in several
        date-folders.
    
        """
        wavedkeys = []
        return super(NL01,self).filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=True,
                              wavedkeys=wavedkeys)
        
    
    
    def check_data(self):
        """
        
        NL01: Checks that data quality is good enough.
        
        **METACODE**
        
        ::
            
            Check common HK values are within safe / nominal margins
            Check voltages in HK match commanded voltages, within margins
        
            f.e.ObsID:
                f.e.CCD:
                    f.e.Q.:
                        measure offsets/means in pre-, img-, over-
                        measure std in pre-, img-, over-
            assess std in pre- is within allocated margins
            (assess offsets in pre- and over- are equal, within allocated  margins)
            assess image-fluences are within allocated margins for each exposure time
            
            plot fluences vs. exposure time
            plot std-pre vs. time
        
            issue any warnings to log
            issue update to report
        
        """
        
        raise NotImplementedError
    
    def prep_data(self):
        """
        
        Takes Raw Data and prepares it for further analysis. 
        
        **METACODE**
        
        ::
            
            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        subtract offset
                        opt: [sub bias frame] 
        
        """
        raise NotImplementedError
        
    def extract_stats(self):
        """
        
        Performs basic analysis: extracts statistics from 
        image regions to later build NLC.
        
        **METACODE**
        
        ::
            
            create segmentation map given grid parameters    
        
            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        f.e. "img-segment": (done elsewhere)
                            measure central value
                            measure variance
    
        """
        
        raise NotImplementedError
    
    
    def produce_NLCs(self):
        """ 
        
        **METACODE**
        
        ::
        
            Obtains Best-Fit Non-Linearity Curve
        
            f.e. CCD:
                f.e. Q:
            
                    [opt] apply correction for source variability (interspersed exposure 
                      with constant exptime)
                    Build NL Curve (NLC) - use stats and exptimes
                    fit poly. shape to NL curve
        
            plot NL curves for each CCD, Q
            report max. values of NL (table)
        
        """
        
        raise NotImplementedError
        
    def do_satCTE(self):
        """
        
        **METACODE**
        
        ::
        
            select ObsIDs with fluence(exptime) >~ 0.5 FWC
            
            f.e. ObsID: 
                CCD: 
                    Q:
                        measure CTE from amount of charge in over-scan relative to fluence
        
            f.e. CCD: 
                Q:
                    get curve of CTE vs. fluence
                    measure FWC from curve in ADU
            
            report FWCs in electrons [via gain in inputs] f.e. CCD, Q (table)
        
        """
        raise NotImplementedError
        
        

    
