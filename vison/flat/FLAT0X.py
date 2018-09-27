#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FLAT0X

Flat-fields acquisition / analysis script

Created on Tue Aug 29 17:32:52 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import datetime
import copy
from collections import OrderedDict

from vison.support import context
from vison.pipe import lib as pilib
from vison.point import lib as polib
from vison.datamodel import scriptic as sc
import FlatFielding as FFing
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
import FLAT0Xaux as FL0Xaux
from vison.support import utils
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

FLAT0X_commvalues = dict(program='CALCAMP',
                         flushes=7, siflsh=1, siflsh_p=500,
                         inisweep=1,
                         vstart=0, vend=2086,
                         exptime=0., shuttr=1,e_shuttr=0,
                         motr_on=0,
                         wave=4,
                         source='flat',
                         comments='')

FLU_lims = dict(CCD1=dict(
    col1=0.25 * 2**16 * (1.+np.array([-0.10, 0.10])),
    col2=0.50 * 2**16 * (1.+np.array([-0.10, 0.10])),
    col3=0.75 * 2**16 * (1.+np.array([-0.10, 0.10]))))
for i in [2, 3]:
    FLU_lims['CCD%i' % i] = copy.deepcopy(FLU_lims['CCD1'])


class FLATS0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimes', ([list], 'Exposure times for each fluence.')),
        ('frames', ([list], 'Number of Frames for each fluence.')),
        ('wavelength', ([int], 'Wavelength')),
    ])))


class FLAT0X(FlatTask):
    """ """

    inputsclass = FLATS0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('indivflats', self.do_indiv_flats),
                         ('masterflat', self.do_master_flat),
                         ('prmask', self.do_prdef_mask)]
        super(FLAT0X, self).__init__(inputs, log, drill, debug)
        self.name = 'FLAT0X'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = FL0Xaux.gt_FL0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       ccdflats='ccdflats',
                                       cdps='cdps')

    def set_inpdefaults(self, **kwargs):
        """ """

        try:
            wavelength = kwargs['wavelength']
        except KeyError:
            wavelength = 800
        try:
            test = kwargs['test']
        except KeyError:
            test = 'FLAT0X'

        t_dummy_F0X = np.array([25., 50., 75])/100.
        tFWCw = self.ogse.profile['tFWC_flat']['nm%i' % wavelength]
        exptimesF0X = (tFWCw * t_dummy_F0X).tolist()  # s
        framesF0X = [80, 60, 30]

        self.inpdefaults = dict(exptimes=exptimesF0X,
                                frames=framesF0X,
                                wavelength=wavelength,
                                test=test)

    def set_perfdefaults(self, **kwargs):
        #wavelength = self.inputs['wavelength']
        super(FLAT0X, self).set_perfdefaults(**kwargs)
        self.perfdefaults['FLU_lims'] = FLU_lims.copy()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds FLAT0X script structure dictionary.

        :param diffvalues: dict, opt, differential values.

        """

        exptimes = self.inputs['exptimes']
        frames = self.inputs['frames']
        wavelength = self.inputs['wavelength']
        test = self.inputs['test']

        FW_ID = self.ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])

        FLAT0X_commvalues['wave'] = FW_IDX
        FLAT0X_commvalues['test'] = test

        assert len(exptimes) == len(frames)
        #assert len(exptimes) == len(flags)

        FLAT0X_sdict = dict()
        for i, exptime in enumerate(exptimes):
            FLAT0X_sdict['col%i' % (i+1,)] = dict(frames=frames[i], exptime=exptimes[i],
                                                  comments='EXP%.1e' % exptime)  # ,comments=flags[i])

        Ncols = len(FLAT0X_sdict.keys())
        FLAT0X_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(FLAT0X_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        FLAT0X_sdict = sc.update_structdict(
            FLAT0X_sdict, commvalues, diffvalues)

        return FLAT0X_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(FLAT0X, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                   wavedkeys=wavedkeys)

    def prepare_images(self):
        """

        FLAT0X: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            [bias structure subtraction, if available]
            cosmetics masking

        """
        super(FLAT0X, self).prepare_images(doExtract=True,
                                           doMask=True, doOffset=True, doBias=True)

    def do_indiv_flats(self):
        """

        **METACODE**

        ::

            Preparation of data for further analysis and 
            produce flat-field for each OBSID.

            f.e. ObsID:
                f.e.CCD:

                    load ccdobj

                    f.e.Q:

                        model 2D fluence distro in image area
                        produce average profile along rows
                        produce average profile along cols

                    save 2D model and profiles in a pick file for each OBSID-CCD
                    divide by 2D model to produce indiv-flat
                    save indiv-Flat to FITS(?), update add filename

            plot average profiles f. each CCD and Q (color coded by time)

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='indivFF', Title='Individual Flat-Fields', level=0)

        # INITIALISATION

        indices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuad = indices.shape

        CCDs = indices.get_vals('CCD')

        ccdindices = copy.deepcopy(self.dd.mx['CCD'].indices)

        self.dd.initColumn('indiv_flats', ccdindices,
                           dtype='S100', valini='None')

        # The heavy-lifting

        if not self.drill:

            dpath = self.inputs['subpaths']['ccdpickles']
            fpath = self.inputs['subpaths']['ccdflats']
                
            vfullinpath_adder = utils.get_path_decorator(dpath)
            
            fccdobj_names = vfullinpath_adder(self.dd.mx['ccdobj_name'][:],'pick')
            
            wavelength = self.inputs['wavelength']

            for iObs in range(nObs):
                
                ObsID = self.dd.mx['ObsID'][iObs]

                for jCCD, CCDk in enumerate(CCDs):
                    ilabel = self.dd.mx['label'][iObs, jCCD]

                    self.dd.mx['indiv_flats'][iObs, jCCD] = 'EUC_FF_%inm_%s_ROE1_%i_%s.fits' %\
                        (wavelength, ilabel, ObsID, CCDk)

            vfulloutpath_adder = utils.get_path_decorator(fpath)
            
            findiv_flats = vfulloutpath_adder(self.dd.mx['indiv_flats'][:])

            ffsettings = dict()
            
            FFing.produce_IndivFlats(fccdobj_names.flatten(), findiv_flats.flatten(),
                                     settings=ffsettings,
                                     runonTests=False,
                                     processes=self.processes)

            # MISSING: 1D profiles. Do as part of produce_IndivFlats, or separate?

            # Show 1D Profiles for each fluence across CCDs:
            #   2 profiles x N-fluences = 2N plots (each with 4Qx3CCD subplots)

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

        if self.report is not None:
            self.report.add_Section(
                keyword='MasterFF', Title='Master Flat-Fields', level=0)

        wavelength = self.inputs['wavelength']
        settings = dict()

        indices = copy.deepcopy(self.dd.indices)

        CCDs = indices.get_vals('CCD')
        
        dpath = self.inputs['subpaths']['ccdflats']
        cdppath = self.inputs['subpaths']['cdps']
        
        vCCDs = self.dd.mx['CCD'][:].copy()

        vlabels = self.dd.mx['label'][:].copy()
        ulabels = np.unique(vlabels)

        if not self.drill:

            PRNU = OrderedDict()

            self.dd.products['MasterFFs'] = OrderedDict()

            for ulabel in ulabels:

                PRNU[ulabel] = OrderedDict()

                self.dd.products['MasterFFs'][ulabel] = OrderedDict()

                for jCCD, CCDk in enumerate(CCDs):

                    PRNU[ulabel][CCDk] = OrderedDict()

                    FFname = 'EUC_FF_%inm_%s_ROE1_%s.fits' % \
                        (wavelength, ulabel, CCDk)

                    FFpath = os.path.join(cdppath, FFname)

                    selix = np.where((vlabels == ulabel) & (vCCDs == CCDk))

                    FFlist = self.dd.mx['indiv_flats'][selix].flatten().copy()
                    
                    vfullinpath_adder = utils.get_path_decorator(dpath)

                    FFlist = vfullinpath_adder(FFlist)

                    # MISSING: proper defects and useful area masking
                    #   (mask-out pre/over scans)
                    

                    FFing.produce_MasterFlat(
                        FFlist, FFpath, mask=None, settings=settings)

                    self.dd.products['MasterFFs'][ulabel][CCDk] = FFpath

                    # Measure Flat-Field (PRNU): per-Q

                    FF = FFing.FlatField(FFpath)

                    iFextension = FF.extnames.index('FLAT')

                    Quads = FF.Quads

                    for Q in Quads:

                        iQ_PRNU = FF.get_stats(Q, sector='img', statkeys=['std'],
                                               ignore_pover=True,
                                               extension=iFextension)[0]

                        PRNU[ulabel][CCDk][Q] = iQ_PRNU

        # REPORT PRNU results

        # SHOW FFs

    def do_prdef_mask(self):
        """
        **METACODE**

        ::

            Produces mask of defects in Photo-Response
            Could use master FF, or a stack of a subset of images (in order
            to produce mask, needed by other tasks, quicker).

            f.e.CCD:
                f.e.Q:
                    produce mask of PR defects
                    save mask of PR defects
                    count dead pixels / columns 

            report PR-defects stats

        """

        raise NotImplementedError
