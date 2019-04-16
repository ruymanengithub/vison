# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: MOT_FF

Brighter-Fatter Analysis
   Using data from test PTC01 (via BF01)
Hard Edge Response in serial / parallel
Bit Correlations (ADC health)

Created on Tue Jul 31 18:04:00 2018

:author: raf

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import warnings
import copy
import string as st
from collections import OrderedDict

from vison.pipe.task import HKKeys
from vison.pipe.task import Task
from vison.flat.BF01 import BF01
from vison.other import MOT_FFaux
from vison.datamodel import inputs
from vison.support.files import cPickleRead, cPickleDumpDictionary
from vison.support import utils
from vison.datamodel import cdp
# END IMPORT

isthere = os.path.exists

#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


class MOT_FF_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
        ('Npix', ([int], 'Number of Pixels (linear) to consider for Covariance Matrix')),
        ('surrogate', ([str], 'Test to use as surrogate'))
    ])))


class MOT_FF(BF01):
    """ """

    inputsclass = MOT_FF_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        super(MOT_FF, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug, cleanafter=cleanafter)
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract_COV', self.extract_COV),
                         ('extract_BF', self.extract_BF)]
        self.subtasks += [('extract_ADC', self.extract_ADC),
                         ('extract_HER', self.extract_HER)]
        
        
        self.name = 'MOT_FF'
        #self.type = 'Simple'
        
        #self.HKKeys = HKKeys
        self.figdict = MOT_FFaux.gt_MOT_FF_figs(self.inputs['surrogate'])
        self.inputs['subpaths'] = dict(figs='figs', 
                   ccdpickles='ccdpickles',
                   covariance='covariance',
                   products='products')

    def check_data(self):
        
        kwargs = dict(figkeys=['MOT_FFchecks_offsets', 'MOT_FFchecks_stds',
                                   'MOT_FFchecks_flu', 'MOT_FFchecks_imgstd'])
        
        Task.check_data(self, **kwargs)
    
    def extract_ADC(self):
        
        return
        
    
    def extract_HER(self):
        """ """
        
        # Loop over selected images
        #    extract s-overscan normalized profile (HER)
        #    extract p-overscan normalized profile (HER)
        # average profiles
        # display profiles
        
        if self.report is not None:
            self.report.add_Section(
                keyword='extract_HER', Title='Hard Edge Response (HER) Profiles Extraction', level=0)
        
        #ObsIDs = self.dd.mx['ObsID'][:].copy()
        labels = self.dd.mx['label'][:, 0].copy()
        flu_med_img = self.dd.mx['flu_med_img'][:].copy()
        
        indices = copy.deepcopy(self.dd.indices)
        nObs, nCCD, nQuad = indices.shape
        CCDs = indices.get_vals('CCD')
        Quads = indices.get_vals('Quad')
        
        dpath = self.inputs['subpaths']['ccdpickles']
        prodspath = self.inputs['subpaths']['products']
        
        # The "Hard" work
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        
        selbool = (flu_med_img > 1.E4) & (flu_med_img < 4.E4)
        selix = np.where(np.sum(selbool,axis=(1,2)) == nCCD*nQuad)
        
        
        profiles_ver = cdp.CDP()
        profiles_ver.header = CDP_header.copy()
        profiles_ver.path = prodspath
        profiles_ver.data = OrderedDict()
        
        profiles_ser = cdp.CDP()
        profiles_ser.header = CDP_header.copy()
        profiles_ser.path = prodspath
        profiles_ser.data = OrderedDict()
        
        
        for CCD in CCDs:
            profiles_ver.data[CCD] = OrderedDict()
            profiles_ser.data[CCD] = OrderedDict()
            for Q in Quads:
                profiles_ver.data[CCD][Q] = OrderedDict(
                        x=OrderedDict(),
                        y=OrderedDict())
                profiles_ser.data[CCD][Q] = OrderedDict(
                        x=OrderedDict(),
                        y=OrderedDict())
        
        for iObs in selix[0]:
            
            #ObsID = ObsIDs[iObs]
            label = labels[iObs]
            
            for jCCD, CCDk in enumerate(CCDs):
                
                ccdobj_f = os.path.join(dpath,'%s.pick' % 
                                      self.dd.mx['ccdobj_name'][iObs, jCCD])
                
                ccdobj = copy.deepcopy(cPickleRead(ccdobj_f))
                
                thresholds = [1.e4,5e4]
                jprofiles_ser = MOT_FFaux.extract_overscan_profiles(ccdobj, thresholds, direction='serial')
                jprofiles_ver = MOT_FFaux.extract_overscan_profiles(ccdobj, thresholds, direction='parallel')
                                
                for kQ, Q in enumerate(Quads):
                    
                    profiles_ser.data[CCDk][Q]['x'][label] = jprofiles_ser[Q]['x'].copy()
                    profiles_ser.data[CCDk][Q]['y'][label] = jprofiles_ser[Q]['y'].copy()
                    
                    profiles_ver.data[CCDk][Q]['x'][label] = jprofiles_ver[Q]['x'].copy()
                    profiles_ver.data[CCDk][Q]['y'][label] = jprofiles_ver[Q]['y'].copy()
                    
        # Figures
        
        fdict_S = self.figdict['MOTFF_HER_ser'][1]
        fdict_S['data'] = profiles_ser.data.copy()
        if self.report is not None:
            self.addFigures_ST(figkeys=['MOTFF_HER_ser'], 
                               dobuilddata=False)
        
        fdict_V = self.figdict['MOTFF_HER_ver'][1]
        fdict_V['data'] = profiles_ver.data.copy()
        if self.report is not None:
            self.addFigures_ST(figkeys=['MOTFF_HER_ver'], 
                               dobuilddata=False)
            
        # Saving profiles
        
        profiles_ver.rootname = 'profs_HER_PARAL_MOT_FF'
        profiles_ver.savetopickle()
        
        profiles_ser.rootname = 'profs_HER_SERIAL_MOT_FF'
        profiles_ser.savetopickle()
        
        self.dd.products['profiles_ver_name'] = profiles_ver.rootname
        self.dd.products['profiles_ser_name'] = profiles_ser.rootname

        
