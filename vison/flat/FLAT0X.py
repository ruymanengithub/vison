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
import string as st
import pandas as pd
from matplotlib.colors import Normalize
from scipy import ndimage

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
        self.CDP_lib = FL0Xaux.CDP_lib.copy()
        self.HKKeys = HKKeys
        self.figdict = FL0Xaux.gt_FL0Xfigs(self.inputs['test'])
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       ccdflats='ccdflats',
                                       products='products')

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
        Quads = indices.get_vals('Quad')

        ccdindices = copy.deepcopy(self.dd.mx['CCD'].indices)

        self.dd.initColumn('indiv_flats', ccdindices,
                           dtype='S100', valini='None')
        
        vCCDs = self.dd.mx['CCD'][:].copy()
        vlabels = self.dd.mx['label'][:].copy()
        ulabels = np.unique(vlabels)
        
        # 1D Profiles
        
        _foox = np.linspace(0,2000,100)
        _fooy = np.ones(10)*1.E4
        
        profs_1D = OrderedDict()
        for dire in ['hor','ver']:
            profs_1D[dire] = OrderedDict()
            for ulabel in ulabels:
                profs_1D[dire][ulabel] = OrderedDict()                
                for CCD in CCDs:
                    profs_1D[dire][ulabel][CCD] = OrderedDict()
                    for Q in Quads:
                        profs_1D[dire][ulabel][CCD][Q] = OrderedDict()
                        profs_1D[dire][ulabel][CCD][Q]['x'] = OrderedDict()
                        profs_1D[dire][ulabel][CCD][Q]['y'] = OrderedDict()
                        
                        for i in range(3):
                            profs_1D[dire][ulabel][CCD][Q]['x']['OBS%i' % i] = _foox.copy()
                            profs_1D[dire][ulabel][CCD][Q]['y']['OBS%i' % i] = _fooy.copy()
        

        # The heavy-lifting
        
        onTests = False # True on TESTS

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
            
            if not onTests:
                FFing.produce_IndivFlats(fccdobj_names.flatten(), findiv_flats.flatten(),
                                         settings=ffsettings,
                                         runonTests=False,
                                         processes=self.processes)

            # 1D profiles.
            
            
            for ulabel in ulabels:
                
                vstart = self.inputs['structure'][ulabel]['vstart']
                truevend = self.inputs['structure'][ulabel]['vend']
                
                iswithpover = truevend==(ccd.NrowsCCD+ccd.voverscan)
                
                vend1D = ccd.NrowsCCD
                
                for jCCD, CCDk in enumerate(CCDs):
                    
                    selix = np.where((vlabels == ulabel) & (vCCDs == CCDk))
                    
                    _indiv_flats = self.dd.mx['indiv_flats'][selix].copy()
                    f_indiv_flats = vfulloutpath_adder(_indiv_flats)
                    
                    _OBSIDS = self.dd.mx['ObsID'][selix[0]].copy()
                    
                    for ix,f_indiv_flat in enumerate(f_indiv_flats):
                        
                        OBSID = _OBSIDS[ix]
                        
                        ccdobj = ccd.CCD(infits=str(f_indiv_flat),getallextensions=True,
                                         withpover=iswithpover)
                        
                        for kQ, Q in enumerate(Quads):
                            
                            _pro1Dhor = ccdobj.get_1Dprofile(Q=Q,orient='hor',area='img',
                                    stacker='mean', vstart=vstart, vend=vend1D,
                                    extension=-1)
                            _pro1Dver = ccdobj.get_1Dprofile(Q=Q,orient='ver',area='img',
                                    stacker='mean',vstart=vstart, vend=vend1D,
                                    extension=-1)
                            
                            profs_1D['hor'][ulabel][CCD][Q]['x']['OBS%i' % OBSID] = \
                                    _pro1Dhor.data['x'].copy()
                            profs_1D['hor'][ulabel][CCD][Q]['y']['OBS%i' % OBSID] = \
                                    _pro1Dhor.data['y'].copy()
                                    
                            profs_1D['ver'][ulabel][CCD][Q]['x']['OBS%i' % OBSID] = \
                                    _pro1Dver.data['x'].copy()
                            profs_1D['ver'][ulabel][CCD][Q]['y']['OBS%i' % OBSID] = \
                                    _pro1Dver.data['y'].copy()           
                

        # Show 1D Profiles (hor/ver) for each OBSID/fluence across CCDs/Quads:
        #   2 profiles x N-fluences = 2N plots (each with 4Qx3CCD subplots)
        
        for ulabel in ulabels:
            for dire in ['hor','ver']:
                
                figkey = 'FL0Xindiv_prof1D_%s_%s' % (dire,ulabel)
                
                self.figdict[figkey]=copy.deepcopy(self.figdict['FL0Xindiv_prof1D_%s_generic' % dire])
                self.figdict[figkey][1]['figname'] = '%s_profs1D_%s_%s.png' % \
                            (self.inputs['test'],dire,ulabel)
                self.figdict[figkey][1]['caption'] = st.replace(self.figdict[figkey][1]['caption'],
                            'PLACEHOLDER',ulabel)
                self.figdict[figkey][1]['meta']['suptitle'] = st.replace(self.figdict[figkey][1]['meta']['suptitle'],
                            'PLACEHOLDER',ulabel)
                
                self.figdict[figkey][1]['data'] = profs_1D[dire][ulabel].copy()
                labelkeys = profs_1D[dire][ulabel]['CCD1'][Q]['x'].keys()
                self.figdict[figkey][1]['data']['labelkeys'] = labelkeys
        
        
        figkeys = []
        for ulabel in ulabels:
            figkeys += [
                    'FL0Xindiv_prof1D_hor_%s' % ulabel,
                    'FL0Xindiv_prof1D_ver_%s' % ulabel,
                    ]

        if self.report is not None:
            self.addFigures_ST(figkeys=figkeys, dobuilddata=False)
            

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
            
        OnTests = True # True on TESTS

        wavelength = self.inputs['wavelength']
        settings = dict(
                WAVEL=wavelength,
                ID=self.ID,
                BLOCKID=self.BLOCKID,
                CHAMBER=self.CHAMBER
                )

        indices = copy.deepcopy(self.dd.indices)

        CCDs = indices.get_vals('CCD')
        nC = len(CCDs)
        Quads = indices.get_vals('Quad')
        nQ = len(Quads)
        
        dpath = self.inputs['subpaths']['ccdflats']
        productspath = self.inputs['subpaths']['products']
        
        vCCDs = self.dd.mx['CCD'][:].copy()

        vlabels = self.dd.mx['label'][:].copy()
        ulabels = np.unique(vlabels)
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        NP = nC * nQ
        
        PRNU_TBs = OrderedDict()
        
        for ulabel in ulabels:
            PRNU_TBs[ulabel] = OrderedDict()
            PRNU_TBs[ulabel]['CCD'] = np.zeros(NP,dtype='int32')
            PRNU_TBs[ulabel]['Q'] = np.zeros(NP,dtype='int32')
            PRNU_TBs[ulabel]['AVFLUENCE'] = np.zeros(NP,dtype='float32')
            PRNU_TBs[ulabel]['PRNU_PC'] = np.zeros(NP,dtype='float32')
            
            
        MFF_2PLOT = OrderedDict()
        for ulabel in ulabels:
            MFF_2PLOT[ulabel] = OrderedDict()
            for CCDk in CCDs:
                MFF_2PLOT[ulabel][CCDk] = OrderedDict()
                for Q in Quads:
                    MFF_2PLOT[ulabel][CCDk][Q] = OrderedDict()
        MFF_p5s = []
        MFF_p95s = []

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

                    FFpath = os.path.join(productspath, FFname)

                    selix = np.where((vlabels == ulabel) & (vCCDs == CCDk))

                    FFlist = self.dd.mx['indiv_flats'][selix].flatten().copy()
                    
                    vfullinpath_adder = utils.get_path_decorator(dpath)

                    FFlist = vfullinpath_adder(FFlist)

                    # MISSING: proper defects and useful area masking
                    #   (mask-out pre/over scans)
                    
                    if OnTests:
                        FFlist = FFlist[0:2]
                    
                    jsettings = copy.deepcopy(settings)
                    
                    for kQ,Q in enumerate(Quads):
                        jsettings['AVFLU_%s' % Q] = \
                                 np.nanmedian(self.dd.mx['flu_med_img'][:,:,kQ][selix])

                    FFing.produce_MasterFlat(
                        FFlist, FFpath, mask=None, settings=jsettings)

                    self.dd.products['MasterFFs'][ulabel][CCDk] = FFpath

                    # Measure Flat-Field (PRNU): per-Q

                    FF = FFing.FlatField(FFpath)

                    iFextension = FF.extnames.index('FLAT')

                    Quads = FF.Quads

                    for jQ, Q in enumerate(Quads):
                        
                        kk = jCCD * nQ + jQ
                        
                        kQ_PRNU = FF.get_stats(Q, sector='img', statkeys=['std'],
                                               ignore_pover=True,
                                               extension=iFextension)[0]
                        
                        PRNU_TBs[ulabel]['CCD'][kk] = jCCD+1
                        PRNU_TBs[ulabel]['Q'][kk] = jQ+1
                        PRNU_TBs[ulabel]['PRNU_PC'][kk] = kQ_PRNU * 100.
                        PRNU_TBs[ulabel]['AVFLUENCE'][kk] = jsettings['AVFLU_%s' % Q]
                        
                        qdata = FF.get_quad(Q,canonical=False,extension=1).copy()
                        sqdata = ndimage.filters.gaussian_filter(qdata,sigma=5.,
                                                                 mode='constant',
                                                                 cval=1.)
                        #MFF_p5s.append(np.percentile(sqdata,5))
                        #MFF_p95s.append(np.percentile(sqdata,95))
                        MFF_2PLOT[ulabel][CCDk][Q]['img'] = sqdata.transpose()

        # REPORT PRNU results
        
        
        PRNU_TB_dddf = OrderedDict()
        for ulabel in ulabels:
            PRNU_TB_dddf['PRNU_%s' % ulabel] = pd.DataFrame.from_dict(PRNU_TBs[ulabel])
        
        prnu_tb_cdp = self.CDP_lib['PRNU_TB']
        prnu_tb_cdp.rootname = prnu_tb_cdp.rootname % \
          (self.inputs['test'],self.inputs['wavelength'])
        prnu_tb_cdp.path = self.inputs['subpaths']['products']
        prnu_tb_cdp.ingest_inputs(
                data = PRNU_TB_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )

        prnu_tb_cdp.init_wb_and_fillAll(header_title='%s (%inm): PRNU TABLE' % \
                (self.inputs['test'],self.inputs['wavelength']))
        self.save_CDP(prnu_tb_cdp)
        self.pack_CDP_to_dd(prnu_tb_cdp, 'PRNU_TB_CDP')
    
        for ulabel in ulabels:
            
            if self.report is not None:
                
                fccd = lambda x: CCDs[x-1]
                fq = lambda x: Quads[x-1]
                ff = lambda x: '%.2f' % x
                fE = lambda x: '%.2E' % x
                
                cov_formatters=[fccd,fq,fE,ff]
                
                caption = '%s (%inm): PRNU TABLE [%s]' % \
                    (self.inputs['test'],self.inputs['wavelength'],ulabel)
                nicecaption = st.replace(caption,'_','\_')
                Ptex = prnu_tb_cdp.get_textable(sheet='PRNU_%s' % ulabel, 
                                                   caption=nicecaption,
                                                   fitwidth=True,
                                                   tiny=True,
                                                   formatters=cov_formatters)
                self.report.add_Text(Ptex)

        # SHOW FFs
        
        for ulabel in ulabels:
            
            self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel] =  copy.deepcopy(self.figdict['FL0Xmeta_MFF_2D'])
            self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['figname'] =\
                        st.replace(self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['figname'],
                                   'PLACEHOLDER',
                                   ulabel)
            self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['caption'] =\
                        st.replace(self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['caption'],
                                   'PLACEHOLDER',
                                   ulabel)
            self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['meta']['suptitle'] =\
                st.replace(self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['meta']['suptitle'],
                                   'PLACEHOLDER',
                                   ulabel)
            
            self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['data'] = MFF_2PLOT[ulabel].copy()
        
            # UPDATING scaling based on data
            
            
            if len(MFF_p5s)>0:
                normfunction = Normalize(vmin=np.nanmin(MFF_p5s),vmax=np.nanmax(MFF_p95s),clip=True)
            else:
                normfunction = Normalize(vmin=0.99,vmax=1.01,clip=False)
            
            self.figdict['FL0Xmeta_MFF_2D_%s' % ulabel][1]['meta']['corekwargs']['norm'] = normfunction
            
        if self.report is not None:
            MFFfigkeys = ['FL0Xmeta_MFF_2D_%s' % ulabel for ulabel in ulabels]
            self.addFigures_ST(figkeys=MFFfigkeys,
                           dobuilddata=False)

            

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
