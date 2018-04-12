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
                         siflush=0,  # sinvflushp=500,
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
        super(BIAS01, self).__init__(inputs, log, drill, debug)
        self.name = 'BIAS01'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        self.HKKeys = HKKeys
        self.figdict = B01aux.B01figs.copy()
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       profiles='profiles',products='products')

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(N=25)

    def set_perfdefaults(self, **kwargs):
        self.perfdefaults = dict()
        self.perfdefaults.update(performance.perf_rdout)

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

    def filterexposures(self, structure, explogf, datapath, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(BIAS01, self).filterexposures(structure, explogf, datapath, OBSID_lims, colorblind=True,
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
        CCDs = DDindices[DDindices.names.index('CCD')].vals

        # The "Hard"-work

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        profs1D2plot = dict()
        profs1D2plot['hor'] = OrderedDict()
        profs1D2plot['ver'] = OrderedDict()

        for CCD in CCDs:
            profs1D2plot['hor']['CCD%i' % CCD] = OrderedDict()
            profs1D2plot['ver']['CCD%i' % CCD] = OrderedDict()
            for Q in Quads:
                profs1D2plot['hor']['CCD%i' % CCD][Q] = []
                profs1D2plot['ver']['CCD%i' % CCD][Q] = []

        if not self.drill:

            ccdpicklespath = self.inputs['subpaths']['ccdpickles']
            profilespath = self.inputs['subpaths']['profiles']

            for iObs in range(nObs):

                OBSID = self.dd.mx['ObsID'][iObs]

                for jCCD, CCD in enumerate(CCDs):

                    vstart = self.dd.mx['vstart'][iObs, jCCD]
                    vend = self.dd.mx['vend'][iObs, jCCD]

                    ccdobj_name = self.dd.mx['ccdobj_name'][iObs, jCCD]

                    if ccdobj_name == 'None':
                        continue  # TESTS

                    fullccdobj_name = os.path.join(
                        ccdpicklespath, '%s.pick' % ccdobj_name)

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    CDP_header['DATE'] = self.get_time_tag()

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
                        RON = onlyRONimg.std()

                        self.dd.mx['RON'][iObs, jCCD, kQ] = RON

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

                        profs1D2plot['hor']['CCD%i' %
                                            CCD][Q].append(hor1Dprof.data.copy())
                        profs1D2plot['hor']['CCD%i' %
                                            CCD][Q][-1].update(dict(label='OBS%i' % OBSID))
                        profs1D2plot['ver']['CCD%i' %
                                            CCD][Q].append(ver1Dprof.data.copy())
                        profs1D2plot['ver']['CCD%i' %
                                            CCD][Q][-1].update(dict(label='OBS%i' % OBSID))

                    # save (2D model and) profiles in a pick file for each OBSID-CCD

                    iprofiles1D.rootname = 'profs1D_%i_CCD%i_BIAS01' % (
                        OBSID, CCD)
                    iprofiles2D.rootname = 'mod2D_%i_CCD%i_BIAS01' % (
                        OBSID, CCD)

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
        
        
        

        # REPORTS
        # Table with median (Best Estimate) values of RON per CCD&Q
        # NEEDS REFACTORING!! (move somewhere else, abstractize it and recycle)

        RON = self.dd.mx['RON'][:].copy()
        RONmsk = np.zeros_like(RON, dtype='int32')
        RONmsk[np.where((np.isclose(RON, 0.)) | (np.isnan(RON)))] = 1
        beRON = np.ma.median(np.ma.masked_array(
            RON, mask=RONmsk), axis=1).data.copy()

        beRONdict = OrderedDict()
        for jCCD, CCD in enumerate(CCDs):
            beRONdict['CCD%i' % CCD] = OrderedDict()
            for kQ, Q in enumerate(Quads):
                beRONdict['CCD%i' % CCD][Q] = beRON[jCCD, kQ]
        beRONdf = pd.DataFrame.from_dict(beRONdict)  # PENDING

        RON_CDP = cdp.Tables_CDP()
        RON_CDP.rootname = 'RON_BIAS01'
        RON_CDP.path = self.inputs['subpaths']['products']
        RON_CDP.ingest_inputs(data=dict(RON=beRONdf),
                              meta=dict(),
                              header=CDP_header.copy())
        RON_CDP.init_workbook()
        RON_CDP.fill_Header(title='BIAS01: RON')
        RON_CDP.fill_Meta()
        RON_CDP.fill_allDataSheets()
        RON_CDP.savehardcopy()
        RON_CDP.savetopickle()

        self.dd.products['RON_CDP'] = os.path.join(
            RON_CDP.path, '%s.pick' % RON_CDP.rootname)

        if self.report is not None:
            beRONtex = RON_CDP.get_textable(sheet='RON', caption='BIAS01: RON')
            self.report.add_Text(beRONtex)

    def meta_analysis(self):
        """

        **METACODE**

        ::

            f. each CCD:
               f. e. Q:
                   stack all ObsIDs to produce Master Bias
                   measure average profile along rows
                   measure average profile along cols
            plot average profiles of Master Bias f. each Q
            produce table(s) with summary of results, include in report
            show Master Bias (image), include in report
            save name of MasterBias to DataDict, report

        """

        raise NotImplementedError
