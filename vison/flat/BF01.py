# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BF01

Brighter-Fatter Analysis
   Using data from test PTC01 or PTC02

Created on Wed Mar 7 10:57:00 2018

:author: raf

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
import string as st
from collections import OrderedDict
import pandas as pd

from vison.support import context
#from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
from vison.datamodel import HKtools
from vison.datamodel import core, ccd, cdp
from vison.image import calibration
import ptc as ptclib
from vison.image import performance
#from FlatTask import FlatTask
from PTC0X import PTC0X
from vison.pipe.task import Task
from vison.datamodel import inputs
import BF01aux
from vison.analysis import Guyonnet15 as G15
from vison.image import covariance as covlib

from vison.support import utils
from vison.support.files import cPickleRead, cPickleDumpDictionary
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


class BF01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([dict, list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
        ('Npix', ([int], 'Number of Pixels (linear) to consider for Covariance Matrix')),
        ('surrogate', ([str], 'Test to use as surrogate'))
    ])))


class BF01(PTC0X):
    """ """

    inputsclass = BF01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(BF01, self).__init__(inputs, log, drill, debug)
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract_COV', self.extract_COV),
                         ('extract_BF', self.extract_BF),
                         ('meta', self.meta_analysis)]        
        self.name = 'BF01'
        #self.type = 'Simple'
        
        #self.HKKeys = HKKeys
        self.CDP_lib = BF01aux.CDP_lib.copy()
        self.figdict = BF01aux.gt_BF01figs(self.inputs['surrogate'])
        self.inputs['subpaths'] = dict(figs='figs', 
                   ccdpickles='ccdpickles',
                   covariance='covariance',
                   products='products')

    def set_inpdefaults(self, **kwargs):
        """ """

        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_inpdefaults(**_kwargs)
        self.inpdefaults['test'] = kwargs['test']
        self.inpdefaults['surrogate'] = kwargs['surrogate']
        self.inpdefaults['Npix'] = 5

    def set_perfdefaults(self, **kwargs):
        # maskerading as PTC0X here...
        _kwargs = copy.deepcopy(kwargs)
        _kwargs['test'] = kwargs['surrogate']
        super(BF01, self).set_perfdefaults(**_kwargs)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds PTC0X script structure dictionary.

        #:param exptimes: list of ints [ms], exposure times.
        #:param frames: list of ints, number of frames for each exposure time.
        #:param wavelength: int, wavelength. Default: 800 nm.
        :param diffvalues: dict, opt, differential values.   

        """
        BF01_sdict = super(BF01, self).build_scriptdict(diffvalues=diffvalues, elvis=elvis)
        Ncols = BF01_sdict['Ncols']
        for i in range(1, Ncols+1):
            BF01_sdict['col%i' % i]['test'] = self.inputs['surrogate']
        return BF01_sdict
        #raise NotImplementedError(
        #    "%s: This Task does not build a script, it uses data from another test" % self.name)

    def filterexposures(self, structure, explog, OBSID_lims):
        """

        """
        return Task.filterexposures(self, structure, explog, OBSID_lims,
                                    wavedkeys=['motr_siz'],colorblind=False,
                                    surrogate=self.inputs['surrogate'])
    
    def check_data(self):
        
        kwargs = dict(figkeys=['BF01checks_offsets', 'BF01checks_stds',
                                   'BF01checks_flu', 'BF01checks_imgstd'])
        
        Task.check_data(self, **kwargs)
    
    
    def prepare_images(self):
        Task.prepare_images(self, doExtract=True, doMask=True,
                                         doOffset=True, doBias=False, doFF=False)
        #super(BF01, self).prepare_images(doExtract=True, doMask=True,
        #                                 doOffset=True, doBias=False, doFF=False)

    def extract_COV(self):
        """

        Performs basic analysis of images:
            - extracts COVARIANCE matrix for each fluence

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extractCOV', Title='Covariance-Matrices Extraction', level=0)
        
        Npix = self.inputs['Npix']

        # labels should be the same accross CCDs. PATCH.
        labels = self.dd.mx['label'][:, 0].copy()
        ulabels = np.unique(labels)
        nL = len(ulabels)

        indices = copy.deepcopy(self.dd.indices)
        nObs, nC, nQ = indices.shape
        CCDs = indices.get_vals('CCD')
        Quads = indices.get_vals('Quad')
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        
        dpath = self.inputs['subpaths']['ccdpickles']
        covpath = self.inputs['subpaths']['covariance']
        #prodspath = self.inputs['subpaths']['products']
        
        profscov_1D = cdp.CDP()
        profscov_1D.header = CDP_header.copy()
        profscov_1D.path = covpath
        profscov_1D.data = OrderedDict()
        
        profscov_1D.data['hor'] = OrderedDict()
        profscov_1D.data['ver'] = OrderedDict()
        
        for CCDk in CCDs:
            for tag in ['hor','ver']:
                profscov_1D.data[tag][CCDk] = OrderedDict()            
            for Q in Quads:
                for tag in ['hor','ver']:
                    profscov_1D.data[tag][CCDk][Q] = OrderedDict()
                    profscov_1D.data[tag][CCDk][Q]['x'] = OrderedDict()
                    profscov_1D.data[tag][CCDk][Q]['y'] = OrderedDict()
        
        NP = nC * nQ * nL
        
        COV_dd = OrderedDict()
        COV_dd['CCD'] = np.zeros(NP,dtype='int32')
        COV_dd['Q'] = np.zeros(NP,dtype='int32')
        COV_dd['col'] = np.zeros(NP,dtype='int32')
        COV_dd['av_mu'] = np.zeros(NP,dtype='float32')
        COV_dd['av_var'] = np.zeros(NP,dtype='float32')
        COV_dd['COV_00'] = np.zeros(NP,dtype='float32')
        COV_dd['COV_01'] = np.zeros(NP,dtype='float32')
        COV_dd['COV_10'] = np.zeros(NP,dtype='float32')
        COV_dd['COV_11'] = np.zeros(NP,dtype='float32')
        
        self.dd.products['COV'] = OrderedDict()
        for CCDk in CCDs:
            self.dd.products['COV'][CCDk] = OrderedDict()
        

        if not self.drill:
            
            vfullinpath_adder = utils.get_path_decorator(dpath)
 
            for jCCD, CCDk in enumerate(CCDs):

                for ku, ulabel in enumerate(ulabels):

                    six = np.where(labels == ulabel)
                    
                    ccdobjNamesList = vfullinpath_adder(self.dd.mx['ccdobj_name'][six[0],jCCD],'pick')
                    
                    ccdobjList = [cPickleRead(item) for item in ccdobjNamesList]
                    
                    icovdict = covlib.get_cov_maps(
                        ccdobjList, Npix=Npix, doTest=False)                    
                    
                    self.dd.products['COV'][CCDk][ulabel] = copy.deepcopy(
                            icovdict)
                    
                    for lQ, Q in enumerate(Quads):
                        jj = jCCD * (nQ*nL) + ku * nQ + lQ
                        
                        COV_dd['CCD'][jj] = jCCD+1
                        COV_dd['Q'][jj] = lQ+1
                        COV_dd['col'][jj] = int(st.replace(ulabel,'col',''))
                        COV_dd['av_mu'][jj] = icovdict['av_mu'][Q]
                        COV_dd['av_var'][jj] = icovdict['av_var'][Q]
                        COV_dd['COV_00'][jj] = icovdict['av_covmap'][Q][0,0]
                        COV_dd['COV_01'][jj] = icovdict['av_covmap'][Q][0,1]
                        COV_dd['COV_10'][jj] = icovdict['av_covmap'][Q][1,0]
                        COV_dd['COV_11'][jj] = icovdict['av_covmap'][Q][1,1]
                        
                        profscov_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                                   np.arange(Npix-1)
                        profscov_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                                   icovdict['av_covmap'][Q][1:,0].copy()
                        
                        profscov_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                                   np.arange(Npix-1)
                        profscov_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                                   icovdict['av_covmap'][Q][1:,0].copy()
                        
        else:
            
            for jCCD, CCDk in enumerate(CCDs):

                for ku, ulabel in enumerate(ulabels):
                    
                    for lQ, Q in enumerate(Quads):
            
                        profscov_1D.data['hor'][CCDk][Q]['x'][ulabel] = \
                                               np.arange(Npix-1)
                        profscov_1D.data['hor'][CCDk][Q]['y'][ulabel] = \
                                               np.arange(Npix-1)
                                    
                        profscov_1D.data['ver'][CCDk][Q]['x'][ulabel] = \
                                               np.arange(Npix-1)
                        profscov_1D.data['ver'][CCDk][Q]['y'][ulabel] = \
                                               np.arange(Npix-1)
            

        # PLOTS
        
        for tag in ['hor','ver']:    
            profscov_1D.data[tag]['labelkeys'] = \
                            profscov_1D.data[tag][CCDs[0]][Quads[0]]['x'].keys()
        
        for tag in ['ver','hor']:
        
            fdict_C = self.figdict['BF01_COV_%s' % tag][1]
            fdict_C['data'] = profscov_1D.data[tag].copy()
            if self.report is not None:
                self.addFigures_ST(figkeys=['BF01_COV_%s' % tag], 
                                   dobuilddata=False)
        
        # Saving Profiles
        
        profscov_1D.rootname = 'profs_COV1D_BF01'
        profscov_1D.savetopickle()
        self.dd.products['profscov_1D_name'] = profscov_1D.rootname
        
        # Table of Results
        
        COV_dddf = OrderedDict(COV = pd.DataFrame.from_dict(COV_dd))
        
        covtable_cdp = self.CDP_lib['COVTABLE']
        covtable_cdp.path = self.inputs['subpaths']['products']
        covtable_cdp.ingest_inputs(
                data = COV_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )

        covtable_cdp.init_wb_and_fillAll(header_title='BF01: COVTABLE')
        self.save_CDP(covtable_cdp)
        self.pack_CDP_to_dd(covtable_cdp, 'COVTABLE_CDP')

        if self.report is not None:
            
            fccd = lambda x: CCDs[x-1]
            fq = lambda x: Quads[x-1]
            fcol = lambda x: 'col%i' % x
            fE = lambda x: '%.2E' % x
            
            cov_formatters=[fccd,fq,fcol]+[fE]*6
            
            COVtex = covtable_cdp.get_textable(sheet='COV', caption='BF01: COV',
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=cov_formatters)
            
            
            self.report.add_Text(COVtex)        
        
        

        return None

    def extract_BF(self):
        """ 

        Performs basic analysis of images:
            - extracts BF matrix for each COV matrix

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extractBF', Title='BF-Matrices Extraction', level=0)

        # labels should be the same accross CCDs. PATCH.
        label = self.dd.mx['label'][:, 0].copy()

        indices = copy.deepcopy(self.dd.indices)
        #nObs, nCCD, nQuad = indices.shape
        CCDs = np.array(indices.get_vals('CCD'))
        Quads = np.array(indices.get_vals('Quad'))

        ulabels = np.unique(label)

        # INITIALISATIONS

        kernel_FWHMx = np.zeros(
            (len(ulabels), len(CCDs), len(Quads)), dtype='float32')
        kernel_FWHMy = np.zeros(
            (len(ulabels), len(CCDs), len(Quads)), dtype='float32')
        kernel_e = np.zeros(
            (len(ulabels), len(CCDs), len(Quads)), dtype='float32')

        fluence = np.zeros(len(ulabels), dtype='float32') + np.nan

        self.dd.products['BF'] = OrderedDict()
        self.dd.products['BF']['ulabels'] = ulabels.copy()
        self.dd.products['BF']['CCDs'] = CCDs.copy()
        self.dd.products['BF']['Quads'] = Quads.copy()

        for jCCD, CCDk in enumerate(CCDs):

            self.dd.products['BF'][CCDk] = OrderedDict()

            for kQ, Q in enumerate(Quads):

                self.dd.products['BF'][CCDk][Q] = OrderedDict()

            self.dd.products['BF']['kernel_FWHMx'] = kernel_FWHMx.copy()
            self.dd.products['BF']['kernel_FWHMy'] = kernel_FWHMy.copy()
            self.dd.products['BF']['kernel_e'] = kernel_e.copy()
            self.dd.products['BF']['fluence'] = fluence.copy()

        if not self.drill:
            
            singlepixmap = np.zeros((101, 101), dtype='float32') + 0.01
            singlepixmap[50, 50] = 1.

            for jCCD, CCDk in enumerate(CCDs):
                
                
                for ix, ulabel in enumerate(ulabels):
                    
                    COV_dict = self.dd.products['COV'][CCDk][ulabel].copy()

                    for kQ, Q in enumerate(Quads):


                        COV_mx = COV_dict['av_covmap'][Q].copy()

                        fluence[ix] = COV_dict['av_mu'][Q].copy()
                        
                        try:
                            Asol_Q, psmooth_Q = G15.solve_for_A_linalg(
                                COV_mx, var=1., mu=1., returnAll=True, doplot=False,
                                verbose=False)
                        
                            kernel_Q = G15.degrade_estatic(singlepixmap, Asol_Q)

                            kerQshape = G15.get_cross_shape_rough(
                                    kernel_Q, pitch=12.)
                        
                            self.dd.products['BF'][CCDk][Q][ulabel] = OrderedDict(Asol=Asol_Q.copy(),
                                                                              psmooth=copy.deepcopy(psmooth_Q),
                                                                              kernel=kernel_Q.copy())
                        

                            kernel_FWHMx[ix, jCCD, kQ] = kerQshape['fwhmx']
                            kernel_FWHMy[ix, jCCD, kQ] = kerQshape['fwhmy']
                            kernel_e[ix, jCCD, kQ] = kerQshape['e']
                        except:
                            self.dd.products['BF'][CCDk][Q][ulabel] = OrderedDict()

        self.dd.products['BF']['kernel_FWHMx'] = kernel_FWHMx.copy()
        self.dd.products['BF']['kernel_FWHMy'] = kernel_FWHMy.copy()
        self.dd.products['BF']['kernel_e'] = kernel_e.copy()
        self.dd.products['BF']['fluence'] = fluence.copy()

        # Plots
        # PENDING

        # REPORTING
        # PENDING

    def meta_analysis(self):
        """

        Analyzes the BF results across fluences.


        """
        
        return # TESTS
        raise NotImplementedError
