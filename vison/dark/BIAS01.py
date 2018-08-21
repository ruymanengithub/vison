#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: BIAS01

Bias-structure/RON analysis script

Created on Tue Aug 29 16:53:40 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
import copy
import string as st
from collections import OrderedDict
import pandas as pd
import unittest

from vison.support import context
from vison.pipe import lib as pilib
from vison.datamodel import scriptic as sc
from vison.datamodel import ccd
from vison.image import calibration
from vison.datamodel import core
import B01aux
#from vison.pipe.task import Task
from DarkTask import DarkTask
from vison.image import performance
from vison.datamodel import inputs, cdp
from vison.support.files import cPickleRead, cPickleDumpDictionary
from vison.support import utils
#from vison.support.report import Report


# END IMPORT

isthere = os.path.exists


HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


BIAS01_commvalues = dict(program='CALCAMP', test='BIAS01',
                         flushes=7, exptime=0., shuttr=0,
                         e_shuttr=0, vstart=0, vend=2086,
                         siflsh=1, siflsh_p=500,
                         chinj=0,
                         s_tpump=0,
                         v_tpump=0,
                         motr_on=0,
                         toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                         wave=4,
                         comments='BIAS')


class BIAS01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('N', ([int], 'Number of Frame Acquisitions.')),
    ])))


class BIAS01(DarkTask):
    """ """

    inputsclass = BIAS01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]

        super(BIAS01, self).__init__(inputs, log, drill, debug)
        self.name = 'BIAS01'
        self.type = 'Simple'
        
        
        self.HKKeys = HKKeys        
        self.figdict = B01aux.B01figs.copy()
        self.CDP_lib = B01aux.CDP_lib.copy()
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       profiles='profiles', products='products')
        

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(N=25)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds BIAS01 script structure dictionary.

        ###:param N: integer, number of frames to acquire.
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """

        N = self.inputs['N']
        BIAS01_sdict = dict(col1=dict(frames=N, exptime=0))

        Ncols = len(BIAS01_sdict.keys())
        BIAS01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(BIAS01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        BIAS01_sdict = sc.update_structdict(
            BIAS01_sdict, commvalues, diffvalues)

        return BIAS01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(BIAS01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                   wavedkeys=wavedkeys)

    def prep_data(self):
        """

        BIAS01: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            cosmetics masking

        """
        super(BIAS01, self).prepare_images(
            doExtract=True, doMask=True, doOffset=True)

    def basic_analysis(self):
        """ 

        BIAS01: Basic analysis of data.

        **METACODE**

        ::

            f. e. ObsID:
               f.e.CCD:

                   load ccdobj of ObsID, CCD

                   with ccdobj, f.e.Q:
                       produce a 2D poly model of bias, save coefficients
                       produce average profile along rows
                       produce average profile along cols
                       # save 2D model and profiles in a pick file for each OBSID-CCD
                       measure and save RON after subtracting large scale structure

            plot RON vs. time f. each CCD and Q
            plot average profiles f. each CCD and Q (color coded by time)

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='BIAS01 Extraction', level=0)

        Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
        self.dd.initColumn('profiles1D_name', Cindices,
                           dtype='S100', valini='None')
        self.dd.initColumn('mods2D_name', Cindices,
                           dtype='S100', valini='None')

        DDindices = copy.deepcopy(self.dd.indices)

        self.dd.initColumn('RON', DDindices, dtype='float32', valini=0.)

        nObs, nCCD, nQuad = DDindices.shape
        Quads = DDindices[2].vals
        CCDs = DDindices.get_vals('CCD')

        # The "Hard"-work

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        profs1D2plot = dict()
        profs1D2plot['hor'] = OrderedDict()
        profs1D2plot['ver'] = OrderedDict()

        for CCDk in CCDs:
            for tag in ['hor', 'ver']:
                profs1D2plot[tag][CCDk] = OrderedDict()
            for Q in Quads:
                for tag in ['hor', 'ver']:
                    profs1D2plot[tag][CCDk][Q] = OrderedDict()
                    profs1D2plot[tag][CCDk][Q]['x'] = OrderedDict()
                    profs1D2plot[tag][CCDk][Q]['y'] = OrderedDict()
        
        prof_all = []

        if not self.drill:
            
            CDP_header['DATE'] = self.get_time_tag()

            ccdpicklespath = self.inputs['subpaths']['ccdpickles']
            profilespath = self.inputs['subpaths']['profiles']

            for iObs in range(nObs):

                OBSID = self.dd.mx['ObsID'][iObs]

                for jCCD, CCDk in enumerate(CCDs):

                    vstart = self.dd.mx['vstart'][iObs, jCCD]
                    vend = self.dd.mx['vend'][iObs, jCCD]

                    ccdobj_name = self.dd.mx['ccdobj_name'][iObs, jCCD]

                    if ccdobj_name == 'None':
                        continue  # TESTS

                    fullccdobj_name = os.path.join(
                        ccdpicklespath, '%s.pick' % ccdobj_name)

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    iprofiles1D = cdp.CDP()
                    iprofiles1D.header = CDP_header.copy()
                    iprofiles1D.path = profilespath
                    iprofiles1D.data = OrderedDict()

                    iprofiles2D = cdp.CDP()
                    iprofiles2D.header = CDP_header.copy()
                    iprofiles2D.path = profilespath
                    iprofiles2D.data = OrderedDict()

                    for kQ, Q in enumerate(Quads):

                        iprofiles1D.data[Q] = OrderedDict()

                        # produce a 2D poly model of bias, [save coefficients?]
                        # measure and save RON after subtracting large scale structure

                        mod2D = ccdobj.get_region2Dmodel(Q=Q, area='all', kind='poly2D',
                                        pdegree=5, doFilter=False, doBin=True, binsize=300,
                                        vstart=vstart, vend=vend, canonical=False, extension=-1)
                        
                        onlyRONimg = ccdobj.get_quad(
                            Q, canonical=False) - mod2D.imgmodel
                        _RON = onlyRONimg.std()

                        self.dd.mx['RON'][iObs, jCCD, kQ] = _RON

                        iprofiles2D.data[Q] = OrderedDict()
                        iprofiles2D.data[Q]['polycoeffs'] = copy.deepcopy(
                            mod2D.polycoeffs)
                        iprofiles2D.data[Q]['polyfunc'] = copy.deepcopy(
                            mod2D.polyfunc)

                        # produce average profile along rows

                        hor1Dprof = ccdobj.get_1Dprofile(Q=Q, orient='hor', area='img', stacker='mean',
                                                         vstart=vstart, vend=vend)

                        # produce average profile along cols

                        ver1Dprof = ccdobj.get_1Dprofile(Q=Q, orient='ver', area='img', stacker='mean',
                                                         vstart=vstart, vend=vend)

                        # save profiles in locally for plotting

                        iprofiles1D.data[Q]['hor'] = copy.deepcopy(hor1Dprof)
                        iprofiles1D.data[Q]['ver'] = copy.deepcopy(ver1Dprof)
                        
                        _profs = dict(hor=hor1Dprof,ver=ver1Dprof)

                        for oritag in ['hor', 'ver']:
                            _profdata = _profs[oritag].data.copy()
                            xorder = np.argsort(_profdata['x'])
                            _y = _profdata['y'][xorder]
                            _x = np.arange(len(_y))
                            profs1D2plot[oritag][CCDk][Q]['y']['OBS%i' %
                                                                 OBSID] = _y.copy()
                            profs1D2plot[oritag][CCDk][Q]['x']['OBS%i' %
                                                                 OBSID] = _x.copy()
                            
                            prof_all += _y.tolist()
                            
                            

                    # save (2D model and) profiles in a pick file for each OBSID-CCD

                    iprofiles1D.rootname = 'profs1D_%i_%s_BIAS01' % (
                        OBSID, CCDk)
                    iprofiles2D.rootname = 'mod2D_%i_%s_BIAS01' % (
                        OBSID, CCDk)

                    #fprofilespickf = os.path.join(profilespath,profilespickf)

                    iprofiles1D.savetopickle()
                    iprofiles2D.savetopickle()

                    # cPickleDumpDictionary(profiles,fprofilespickf)

                    self.dd.mx['profiles1D_name'][iObs,
                                                  jCCD] = iprofiles1D.rootname
                    self.dd.mx['mods2D_name'][iObs,
                                              jCCD] = iprofiles2D.rootname

        # PLOTS
        # profiles 1D (Hor & Ver) x (CCD&Q)
        # histograms of RON per CCD&Q
        
        if len(prof_all) != 0:        
             ylim_1D = [np.percentile(prof_all,5)-20.,
                       np.percentile(prof_all,95)+20.]
        else:
            ylim_1D= [0, 2.**16]

        figkeys1 = ['B01basic_prof1D_hor', 'B01basic_prof1D_ver',
                    'B01basic_histosRON']
        self.figdict['B01basic_prof1D_hor'][1]['data'] = profs1D2plot['hor']
        self.figdict['B01basic_prof1D_hor'][1]['meta']['ylim'] = ylim_1D
        self.figdict['B01basic_prof1D_ver'][1]['data'] = profs1D2plot['ver']
        self.figdict['B01basic_prof1D_ver'][1]['meta']['ylim'] = ylim_1D
        
        def _load_histo(histos,data,label):
           
           hist, bin_edges = np.histogram(data,bins=10,normed=False)
           histos['x'][label] = bin_edges.copy()
           histos['y'][label] = hist.copy()
           
           return histos
        
        RON_histos = OrderedDict()
        for jCCD, CCDk in enumerate(CCDs):
            RON_histos[CCDk] = OrderedDict()
            for kQ, Q in enumerate(Quads):
                RON_histos[CCDk][Q] = OrderedDict()
                Rh_cq = RON_histos[CCDk][Q]
                Rh_cq['x'] = OrderedDict()
                Rh_cq['y'] = OrderedDict()
                Rh_cq = _load_histo(Rh_cq,self.dd.mx['RON'][:,jCCD,kQ],'ALL')
                Rh_cq = _load_histo(Rh_cq,self.dd.mx['std_pre'][:,jCCD,kQ],'PRE')
                Rh_cq = _load_histo(Rh_cq,self.dd.mx['std_ove'][:,jCCD,kQ],'OVE')
        
        
        
        self.figdict['B01basic_histosRON'][1]['data'] = RON_histos.copy()
        
        
        if self.report is not None:
            self.addFigures_ST(figkeys=figkeys1, dobuilddata=False)
        #figkeys2 = ['B01basic_histosRON']
        # self.figdict['B01basic_histosRON'][1]['data'] = None # PENDING
        # self.addFigures_ST(figkeys=figkeys2,dobuilddata=True) # PENDING

        # REPORTS

        # Table with median (Best Estimate) values of RON per CCD&Q
        # NEEDS REFACTORING?! (move somewhere else, abstractize it and recycle)
        
        RONdct = OrderedDict()
        RONdct['RON_ALL'] = self.dd.mx['RON'][:].copy()
        RONdct['RON_PRE'] = self.dd.mx['std_pre'][:].copy()
        RONdct['RON_OVE'] = self.dd.mx['std_ove'][:].copy()
        
        ron_cdp = self.CDP_lib['RON']
        ron_cdp.path = self.inputs['subpaths']['products']
        ron_cdp.ingest_inputs(mx_dct=RONdct.copy(),
                              CCDs=CCDs,
                              Quads=Quads,
                              meta=dict(),
                              header=CDP_header.copy())
        
        ron_cdp.init_wb_and_fillAll(header_title='BIAS01: RON')

        self.save_CDP(ron_cdp)
        self.pack_CDP_to_dd(ron_cdp, 'RON_CDP')

        if self.report is not None:
            beRONtex = ron_cdp.get_textable(sheet='RON_ALL', caption='BIAS01: RON (quadrant, minus model)',
                                            longtable=False)
            self.report.add_Text(beRONtex)
            
            bePRERONtex = ron_cdp.get_textable(sheet='RON_PRE', caption='BIAS01: RON (pre-scan)',
                                               longtable=False)
            self.report.add_Text(bePRERONtex)
            
            beOVERONtex = ron_cdp.get_textable(sheet='RON_OVE', caption='BIAS01: RON (over-scan)',
                                               longtable=False)
            self.report.add_Text(beOVERONtex)


    def meta_analysis(self):
        """

        **METACODE**

        ::

            f. each CCD:
               f. e. Q:
                   stack all ObsIDs to produce Master Bias
                   measure average profile along rows
                   measure average profile along cols
            plot average profiles of Master Bias(s) f. each CCD,Q
            produce table(s) with summary of results, include in report
            show Master Bias(s) (3 images), include in report
            save name of MasterBias(s) to DataDict, report

        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='BIAS01 Meta-Analysis', level=0)
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuad = DDindices.shape
        Quads = DDindices[2].vals
        CCDs = DDindices.get_vals('CCD')
        
        # The "Hard"-work

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        if not self.drill:
            
            CDP_header['DATE'] = self.get_time_tag()
            ccdpicklespath = self.inputs['subpaths']['ccdpickles']
            profilespath = self.inputs['subpaths']['profiles']
            
            for jCCD, CCDk in enumerate(CCDs):
                
                
                ccdobj_names = self.dd.mx['ccdobj_name'][:, jCCD]
                fullccdobj_names = map(lambda x: os.path.join('%s.pick' % x), 
                                       ccdobj_names)
                
                ccdobjs_list = [copy.deepcopy(cPickleRead(item)) \
                                for item in fullccdobj_names]
                
                stop()


class Test(unittest.TestCase):
    """
    Unit tests for the BIAS01 class.
    """

    def setUp(self):

        inputs = dict()
        self.b01 = BIAS01(inputs, log=None, drill=True, debug=False)

    def test_check_data(self):
        """

        :return: None
        """
        self.b01.check_data()


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(Test)
    unittest.TextTestRunner(verbosity=3).run(suite)
