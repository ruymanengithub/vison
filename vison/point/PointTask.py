#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:55:00 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import copy
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import pandas as pd
import string as st

from vison.image import sextractor as sex
from vison.datamodel import compliance as complimod
from vison.pipe.task import Task
from vison.point import lib as polib
from vison.point import startracker as strackermod
from vison.datamodel import core, ccd, inputs
#from vison.pipe import lib as pilib
from vison.support import context
from vison.support import utils
from vison.xtalk import xtalk, opt_xtalk
from vison.support.files import cPickleRead
# END IMPORT

BGD_lims = OrderedDict(CCD1=OrderedDict(E=[-5., 10.]))
for Q in ['F', 'G', 'H']:
    BGD_lims['CCD1'][Q] = copy.deepcopy(BGD_lims['CCD1']['E'])
for iCCD in [2, 3]:
    BGD_lims['CCD%i' % iCCD] = copy.deepcopy(BGD_lims['CCD1'])


class Point_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('offsetxy', ([list, tuple], 'Point Sources Offset wrt nominal positions (x,y)/[x,y]'))
    ])))


class PointTask(Task):

    stampw = polib.stampw

    def __init__(self, *args, **kwargs):
        super(PointTask, self).__init__(*args, **kwargs)
        
        CCDs = self.ogse.startrackers.keys()
        offsetxy = np.array(self.inputs['offsetxy'])
        if not np.isclose(np.abs(np.sqrt(np.dot(offsetxy,offsetxy))),0.):
            for jCCD, CCDk in enumerate(CCDs):
                simmx = self.ogse.startrackers[CCDk].get_similaritymx(1.0, 
                                     0., offsetxy)
                self.ogse.startrackers[CCDk].apply_patt_transform(simmx)

    def check_data(self, **kwargs):
        """ """
        
        
        
        test = self.inputs['test']
        if 'PSF01' in test:
            _kwargs = dict(figkeys=['PSF0Xchecks_offsets', 'PSF0Xchecks_stds',
                                    'PSF0Xchecks_bgd', 'PSF0Xchecks_fluence',
                                    'PSF0Xchecks_fwhmx', 'PSF0Xchecks_fwhmy'])
        elif 'PSF02' in test:
            _kwargs = dict(figkeys=['PSF0Xchecks_offsets', 'PSF0Xchecks_stds',
                                    'PSF0Xchecks_bgd', 'PSF0Xchecks_fluence',
                                    'PSF0Xchecks_fwhmx', 'PSF0Xchecks_fwhmy'])
        elif 'FOCUS00' in test:
            _kwargs = dict(figkeys=['F00checks_offsets', 'F00checks_stds',
                                    'F00checks_bgd', 'F00checks_fluence',
                                    'F00checks_fwhmx', 'F00checks_fwhmy'])
        elif 'PSFLUX00' in test:
            _kwargs = dict(figkeys=['PSF0Xchecks_offsets', 'PSF0Xchecks_stds',
                                    'PSF0Xchecks_bgd', 'PSF0Xchecks_fluence',
                                    'PSF0Xchecks_fwhmx', 'PSF0Xchecks_fwhmy'])
        kwargs.update(_kwargs)

        Task.check_data(self, **kwargs)

    def get_checkstats_ST(self, **kwargs):
        """ """

        # Initialize new columns

        Qindices = copy.deepcopy(self.dd.indices)
    
        if 'Quad' not in Qindices.names:
            Qindices.append(core.vIndex('Quad', vals=context.Quads))
            
        if 'Spot' in Qindices.names:
            Qindices.pop(Qindices.find('Spot'))

        valini = 0.

        newcolnames_off = ['offset_pre', 'offset_ove']
        
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off, Qindices,
                               dtype='float32', valini=valini)

        self.dd.initColumn('bgd_img', Qindices, dtype='float32', valini=valini)

        newcolnames_std = ['std_pre', 'std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std, Qindices,
                               dtype='float32', valini=valini)

        Sindices = copy.deepcopy(Qindices)
        Sindices.append(core.vIndex(
            'Spot', vals=strackermod.starnames))

        self.dd.initColumn('chk_x', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('chk_x_ccd', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('chk_y', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('chk_y_ccd', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('chk_peak', Sindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('chk_fluence', Sindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('chk_fwhmx', Sindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('chk_fwhmy', Sindices,
                           dtype='float32', valini=valini)

        chkkeycorr = dict(chk_x='x', chk_y='y', 
                          chk_x_ccd='x_ccd', chk_y_ccd = 'y_ccd',
                          chk_peak='peak', chk_fluence='fluence',
                          chk_fwhmx='fwhmx', chk_fwhmy='fwhmy')
        
        nObs = Qindices[0].len
        CCDs = Qindices.get_vals('CCD')
        Quads = Qindices.get_vals('Quad')
        Spots = Sindices.get_vals('Spot')

        # Get statistics in different regions

        if not self.drill:

            strackers = self.ogse.startrackers

            psCCDcoodicts = OrderedDict(names=strackers['CCD1'].starnames)

            for jCCD, CCDk in enumerate(CCDs):
                psCCDcoodicts[CCDk] = strackers[CCDk].get_allCCDcoos(
                    nested=True)

            for iObs in range(nObs):
                
                for jCCD, CCDk in enumerate(CCDs):

                    dpath = self.dd.mx['datapath'][iObs, jCCD]
                    ffits = os.path.join(dpath, '%s.fits' %
                                         self.dd.mx['File_name'][iObs, jCCD])
                    vstart = self.dd.mx['vstart'][iObs][jCCD]
                    vend = self.dd.mx['vend'][iObs][jCCD]

                    ccdobj = ccd.CCD(ffits)

                    for kQ, Quad in enumerate(Quads):

                        for reg in ['pre', 'ove']:
                            stats_bias = ccdobj.get_stats(Quad, sector=reg, statkeys=['median', 'std'], trimscan=[5, 5],
                                                          ignore_pover=True, extension=-1, VSTART=vstart, VEND=vend)
                            self.dd.mx['offset_%s' %
                                       reg][iObs, jCCD, kQ] = stats_bias[0]
                            self.dd.mx['std_%s' %
                                       reg][iObs, jCCD, kQ] = stats_bias[1]

                        # To measure the background we mask out the sources

                        alt_ccdobj = copy.deepcopy(ccdobj)

                        mask_sources = polib.gen_point_mask(
                            CCDk, Quad, width=self.stampw, sources='all',
                            coodict=psCCDcoodicts.copy())

                        alt_ccdobj.get_mask(mask_sources)

                        alt_ccdobj.sub_offset(Quad, method='row', scan='pre', trimscan=[5, 5],
                                              ignore_pover=True, extension=-1)

                        imgstats = alt_ccdobj.get_stats(Quad, sector='img', statkeys=['median'], trimscan=[5, 5],
                                                        ignore_pover=True, extension=-1)

                        self.dd.mx['bgd_img'][iObs, jCCD, kQ] = imgstats[0]

                        alt_ccdobj = None
                        
                        
                        for xSpot, SpotName in enumerate(Spots):

                            #coo = polib.Point_CooNom[CCDk][Quad][SpotName]
                            coo = psCCDcoodicts[CCDk][Quad][SpotName]

                            spot = polib.extract_spot(ccdobj, coo, Quad, log=self.log,
                                                      stampw=self.stampw)

                            try:
                                res_bas = spot.measure_basic(
                                    rap=10, rin=15, rout=-1)
                            except:
                                res_bas = dict(zip(chkkeycorr.values(), np.zeros(
                                    len(chkkeycorr), dtype='float32')))

                            for chkkey in chkkeycorr:
                                self.dd.mx[chkkey][iObs, jCCD, kQ,
                                                   xSpot] = res_bas[chkkeycorr[chkkey]]
                            
                            #if (CCDk == 'CCD2') and (Quad=='H') and (SpotName=='ALPHA'):
                            #    stop()
                            

    def check_stat_perCCDQSpot(self, arr, lims, CCDs=['CCD1', 'CCD2', 'CCD3']):
        """ """
        spotnames = strackermod.starnames
        Qs = ccd.Quads

        if isinstance(lims[CCDs[0]][Qs[0]][spotnames[0]], (dict, OrderedDict)):
            colnames = lims[CCDs[0]][Qs[0]][spotnames[0]].keys()
            indexer = self.dd.mx['label'][:].copy()
        elif isinstance(lims[CCDs[0]][Qs[0]][spotnames[0]], (list, tuple)):
            colnames = None
            indexer = None

        compliance = complimod.ComplianceMX_CCDQColSpot(spotnames,
                                                        colnames=colnames,
                                                        indexer=indexer,
                                                        CCDs=CCDs,
                                                        Qs=Qs,
                                                        lims=lims.copy())

        compliance.check_stat(arr)
        return compliance

    def check_metrics_ST(self, **kwargs):
        """         
        TO-CHECK:
            - offset levels (pre and over-scan), abs. and relative
            - RON in pre and overscan
            - background level in image area
            - spot fluences
            - spot sizes

        """

        #test = self.inputs['test']

        Xindices = self.dd.indices
        CCDs = Xindices.get_vals('CCD')

        if self.report is not None:
            self.report.add_Section(
                keyword='check_metrics', Title='Offsets, RON, Spots metrics', level=1)

        # absolute value of offsets

        offsets_lims = self.perflimits['offsets_lims']
        for reg in ['pre', 'ove']:
            arr = self.dd.mx['offset_%s' % reg]
            _compliance_offsets = self.check_stat_perCCDandQ(
                arr, offsets_lims, CCDs)
            
            self.addComplianceMatrix2Self(_compliance_offsets,'offsets_lims_%s' % reg)

            if not self.IsComplianceMatrixOK(_compliance_offsets):
                self.dd.flags.add('POORQUALDATA')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)

        # cross-check of offsets: referred to pre-scan

        offsets_gradients = self.perflimits['offsets_gradients']
        for ireg, reg in enumerate(['ove']):
            _lims = dict()
            for CCDk in CCDs:
                _lims[CCDk] = offsets_gradients[CCDk][reg]
            arr = self.dd.mx['offset_%s' % reg][:]-self.dd.mx['offset_pre'][:]
            _xcheck_offsets = self.check_stat_perCCDandQ(arr, _lims, CCDs)
            
            self.addComplianceMatrix2Self(_xcheck_offsets,'offsets_gradients_%s' % reg)

            if not self.IsComplianceMatrixOK(_xcheck_offsets):
                self.dd.flags.add('POORQUALDATA')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _xcheck_offsets, label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _xcheck_offsets, label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)

        # absolute value of std

        regs_std = ['pre', 'ove']
        RONs_lims = self.perflimits['RONs_lims']

        for reg in regs_std:
            _compliance_std = self.check_stat_perCCDandQ(
                self.dd.mx['std_%s' % reg], RONs_lims, CCDs)
            
            self.addComplianceMatrix2Self(_compliance_std,'std_%s' % reg)

            if not self.IsComplianceMatrixOK(_compliance_std):
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('RON_OOL')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)

        # Background Level

        BGD_lims = self.perflimits['BGD_lims']  # dict
        _compliance_bgd = self.check_stat_perCCDandQ(
            self.dd.mx['bgd_img'], BGD_lims, CCDs)
        
        self.addComplianceMatrix2Self(_compliance_bgd,'bgd_img')

        if not self.IsComplianceMatrixOK(_compliance_bgd):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('BGD_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_bgd, label='COMPLIANCE BGD:')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_bgd, label='COMPLIANCE BGD:')

        # Spot FWHM-(x**2+y**2_)**0.5

        FWHM_lims = self.perflimits['FWHM_lims']  # dict
        chk_fwhm = (self.dd.mx['chk_fwhmx'][:]**2. +
                    self.dd.mx['chk_fwhmy'][:]**2.)**0.5
        _compliance_fwhm = self.check_stat_perCCDQSpot(
            chk_fwhm, FWHM_lims, CCDs)
        
        self.addComplianceMatrix2Self(_compliance_fwhm,'fwhm')

        if not self.IsComplianceMatrixOK(_compliance_fwhm):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FOCUS_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_fwhm, label='COMPLIANCE FWHM(x2+y2)**0.5:')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_fwhm, label='COMPLIANCE FWHM(x2+y2)**0.5:')

        # Spot Fluence

        Flu_lims = self.perflimits['Flu_lims']  # dict
        _compliance_flu = self.check_stat_perCCDQSpot(
            self.dd.mx['chk_fluence'], Flu_lims, CCDs)
        
        self.addComplianceMatrix2Self(_compliance_flu,'fluence')

        if not self.IsComplianceMatrixOK(_compliance_flu):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCE_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_flu, label='COMPLIANCE FLUENCE:')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_flu, label='COMPLIANCE FLUENCE:')
            
            
        # Peak Saturation Times
        
        fluences = np.mean(self.dd.mx['chk_fluence'][:].copy(),axis=-1)
        exptime = self.dd.mx['exptime'][:].copy()
            
        ixnozero = np.where(exptime[:,0]>0)
            
        _f = np.squeeze(fluences[ixnozero,...])
        _e = np.expand_dims(np.squeeze(exptime[ixnozero,...]),axis=-1)
            
        fluxes = np.nanmean(_f/_e,axis=0) # CRUDE!
        
        pkfluxes = fluxes / (2.*np.pi*(2./2.355)**2.) # assuming a fwhm of 2 pixels at best focus
        pksat_times = 2.**16/pkfluxes
        pksat_times = np.expand_dims(pksat_times,axis=0)      
        
        
        exp_sat_time = self.ogse.profile['tFWC_point']['nm%i' % self.inputs['wavelength']]
        sat_time_lims = dict()
        for CCD in CCDs:
            sat_time_lims[CCD] = (exp_sat_time * np.array([0.9,1.1])).tolist()
            
        _compliance_flux = self.check_stat_perCCDandQ(
                pksat_times, sat_time_lims, CCDs)
        
        self.addComplianceMatrix2Self(_compliance_flux,'peakflux')        

        if not self.IsComplianceMatrixOK(_compliance_flux):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUX_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_flux, label='COMPLIANCE SATURATION TIME (PEAK FLUX)')
        if self.report is not None:
            self.addComplianceMatrix2Report(
            _compliance_flux, label='COMPLIANCE SATURATION TIME (PEAK FLUX)',
            caption='Saturation times in seconds. Assuming a FWHM of 2 pixels (best focus).')
        
    def relock(self,):
        """ """
        
        if 'LOCK_TB_CDP' not in self.dd.products:
            return
        
        ddindices = copy.deepcopy(self.dd.indices)
        CCDs = ddindices.get_vals('CCD')
        
        lock_tb_cdp = cPickleRead(self.dd.products['LOCK_TB_CDP'])
        lock_tb = lock_tb_cdp['data']['LOCK_TB']
        
        meta = lock_tb_cdp['meta']
        scale_lims = meta['scale_lims']
        rot_lims  = meta['rot_lims']
        trans_lims = meta['trans_lims']
        maxSTDpix = meta['maxSTDpix']
        
        for jCCD,CCDk in enumerate(CCDs):
            
            ixselCCD = np.where(lock_tb['CCD'][:]== jCCD+1)[0][0]
            
            ROTATION = lock_tb['ROTATION'][ixselCCD]
            SCALE = lock_tb['SCALE'][ixselCCD]
            TRANS_X = lock_tb['TRANS_X'][ixselCCD]
            TRANS_Y = lock_tb['TRANS_Y'][ixselCCD]
            
            # CHECK COMPLIANCE of transformations against allowed changes in
            #       scale, rotation and translation
            
            rotation_inlims = np.all((rot_lims[0] <= ROTATION) & \
                                 (ROTATION<= rot_lims[1]))
            
            rotation_ok = rotation_inlims and \
                np.std(ROTATION * 4096.) < maxSTDpix
            
            scale_inlims = np.all((scale_lims[0] <= SCALE) & \
                                 (SCALE<= scale_lims[1]))
            
            scale_ok = scale_inlims and \
                np.std(SCALE * 4096.) < maxSTDpix
            
            trans_mod = np.sqrt(TRANS_X**2.+TRANS_Y**2.)
            
            trans_inlims = np.all((trans_mod>=trans_lims[0]) & \
                               (trans_mod<=trans_lims[1]))
            
            trans_ok = trans_inlims and \
                np.std(trans_mod) < maxSTDpix
            
            all_ok = rotation_ok and scale_ok and trans_ok
            
            if all_ok:
                
                if self.log is not None:
                    self.log.info('%s: Relocating Point Sources!' % CCDk)
                
                self.ogse.load_startrackers(withpover=True)
                
                TRANSLATION = (TRANS_X,TRANS_Y)
                
                simmx = self.ogse.startrackers[CCDk].get_similaritymx(SCALE, 
                                     ROTATION, TRANSLATION)
                    
                self.ogse.startrackers[CCDk].apply_patt_transform(simmx)
    
    def lock_on_stars(self, iObs=0, sexconfig=None):
        """ """
                
        devel = False # TESTS
        
        if self.report is not None:
            self.report.add_Section(
                keyword='lock', Title='Stars Locking', level=0)
        
        _sexconfig = dict(MINAREA=2.,
             DET_THRESH=14.,
             MAG_ZERPOINT=20.,
             SATUR_LEVEL=65535.,
             SEEING_FWHM=1.,
             PIXEL_SCALE=1.,
             GAIN=1.
             )
        
        if sexconfig is not None:
            _sexconfig.update(sexconfig)
        
        class _Transf(object):
            rotation = 1.4*np.pi
            scale = 2.0
            translation = [2000.0,2000.0]
        
        Qindices = copy.deepcopy(self.dd.indices)
        CCDs = Qindices.get_vals('CCD')
        
        scale_lims = [0.95,1.05] # adim.
        rot_lims = np.array([-1.,1.])*1./180.*np.pi # radians
        trans_lims = [-1.,1000.] # pixels
        maxSTDpix = 50. # pixel-rms variations in rot/scale/trans across CCDs
        
        # INITIALISATIONS
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        LOCK_TB = OrderedDict()
        NCCDs = len(CCDs)
        LOCK_TB['CCD'] = np.zeros(NCCDs,dtype='int32')
        LOCK_TB['NMATCH'] = np.zeros(NCCDs,dtype='int32')
        LOCK_TB['SCALE'] = np.zeros(NCCDs,dtype='float32')
        LOCK_TB['ROTATION'] = np.zeros(NCCDs,dtype='float32')
        LOCK_TB['TRANS_X'] = np.zeros(NCCDs,dtype='float32')
        LOCK_TB['TRANS_Y'] = np.zeros(NCCDs,dtype='float32')
        
        
        strackers = self.ogse.startrackers
        
        transfs_dict = OrderedDict()
        
        for jCCD, CCDk in enumerate(CCDs):
        
            dpath = self.dd.mx['datapath'][iObs, jCCD]
            ffits = os.path.join(dpath, '%s.fits' %
                                         self.dd.mx['File_name'][iObs, jCCD])
            #vstart = self.dd.mx['vstart'][iObs][jCCD]
            #vend = self.dd.mx['vend'][iObs][jCCD]
            
            # DEFAULTS
            transf = _Transf()
            s_list = np.arange(5).tolist()
            t_list = copy.deepcopy(s_list)
            
            if not devel:
                
                ccdobj = ccd.CCD(ffits)
                
                img = ccdobj.extensions[-1].data.T.copy()
                SExCatroot = 'StarFinder_%s' % CCDk
                SExCat = sex.easy_run_SEx(img, SExCatroot, 
                            sexconfig=_sexconfig,                                      
                            cleanafter=True)
                
                sel = np.where((SExCat['ELONGATION']<2.) &
                               (SExCat['A_IMAGE']<5.) &
                               (SExCat['A_IMAGE']>0.2) &
                               (SExCat['B_IMAGE']<5.) &
                               (SExCat['B_IMAGE']>0.2) &
                               (SExCat['FLUX_AUTO']>5000.) &
                               (SExCat['ISOAREA_IMAGE']>3.) &
                               (SExCat['ISOAREA_IMAGE']<100.)
                               )
                
                if len(sel[0])> 0:
                    
                    X_IMAGE = SExCat['X_IMAGE'][sel].copy() - 1.
                    Y_IMAGE = SExCat['Y_IMAGE'][sel].copy() - 1.
                                    
                    
                    X_PHYS, Y_PHYS = self.ccdcalc.cooconv_CCD_2_Phys(X_IMAGE,Y_IMAGE)                    
                    whatQ = self.ccdcalc.get_Q_of_CCDcoo(X_IMAGE,Y_IMAGE)[0]
                    
                    ixnonan = np.where(~(np.isnan(X_PHYS) | np.isnan(Y_PHYS)) &
                            (whatQ != 'G'))
                    X_PHYS = X_PHYS[ixnonan]
                    Y_PHYS = Y_PHYS[ixnonan]
                    
                    
                    try:
                        transf, (s_list, t_list) = strackers[CCDk].find_patt_transform(X_PHYS,Y_PHYS,
                                  Full=True, discardQ=['G'],debug=False)
                    except:
                    
                        if self.log is not None:
                            self.log.info('Failed Locking Stars on %s' % CCDk)
                else:
                    if self.log is not None:
                            self.log.info('Did not Find enough Stars to Lock-on on %s' % CCDk)
            
            transfs_dict[CCDk] = copy.deepcopy(transf)
            
            LOCK_TB['CCD'][jCCD] = jCCD + 1
            LOCK_TB['NMATCH'][jCCD] = len(s_list)
            LOCK_TB['SCALE'][jCCD] = transf.scale
            LOCK_TB['ROTATION'][jCCD] = transf.rotation
            LOCK_TB['TRANS_X'][jCCD] = transf.translation[0]
            LOCK_TB['TRANS_Y'][jCCD] = transf.translation[1]
            
        
        # REPORT LOCK_TB
        
        LOCK_TB_dddf =OrderedDict(LOCK_TB = pd.DataFrame.from_dict(LOCK_TB))
        
        meta_cdp = OrderedDict()
        meta_cdp['scale_lims'] = scale_lims
        meta_cdp['rot_lims'] = rot_lims 
        meta_cdp['trans_lims'] = trans_lims
        meta_cdp['maxSTDpix'] = maxSTDpix
        meta_cdp.update(_sexconfig)
        
        lock_tb_cdp = self.CDP_lib['LOCK_TB']
        lock_tb_cdp.rootname = lock_tb_cdp.rootname % self.inputs['test']
        lock_tb_cdp.path = self.inputs['subpaths']['products']
        lock_tb_cdp.ingest_inputs(
                data = LOCK_TB_dddf.copy(),
                meta=meta_cdp.copy(),
                header=CDP_header.copy()
                )
        
        
        lock_tb_cdp.init_wb_and_fillAll(header_title='%s: STARS LOCK TABLE' % \
                self.inputs['test'])
        self.save_CDP(lock_tb_cdp)
        self.pack_CDP_to_dd(lock_tb_cdp, 'LOCK_TB_CDP')
        
        if self.report is not None:
            
            fccd = lambda x: CCDs[x-1]
            fi = lambda x: '%i' % x
            ff = lambda x: '%.2f' % x
            fE = lambda x: '%.2E' % x
            
            cov_formatters=[fccd,fi,ff,fE,ff,ff]
            
            caption = '%s: Star Lock Table' % self.inputs['test']
            nicecaption = st.replace(caption,'_','\_')
            Ltex = lock_tb_cdp.get_textable(sheet='LOCK_TB', caption=nicecaption,
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=cov_formatters)
            
            
            self.report.add_Text(Ltex)        
    
        
        
        # CHECK COMPLIANCE of transformations against allowed changes in
        #       scale, rotation and translation
        
        rotation_inlims = np.all((rot_lims[0] <= LOCK_TB['ROTATION']) & \
                             (LOCK_TB['ROTATION']<= rot_lims[1]))
        
        rotation_ok = rotation_inlims and \
            np.std(LOCK_TB['ROTATION'] * 4096.) < maxSTDpix
        
        scale_inlims = np.all((scale_lims[0] <= LOCK_TB['SCALE']) & \
                             (LOCK_TB['SCALE']<= scale_lims[1]))
        
        scale_ok = scale_inlims and \
            np.std(LOCK_TB['SCALE'] * 4096.) < maxSTDpix
        
        trans_mod = np.sqrt(LOCK_TB['TRANS_X']**2.+LOCK_TB['TRANS_Y']**2.)
        
        trans_inlims = np.all((trans_mod>=trans_lims[0]) & \
                           (trans_mod<=trans_lims[1]))
        
        trans_ok = trans_inlims and \
            np.std(trans_mod) < maxSTDpix
        
        all_ok = rotation_ok and scale_ok and trans_ok
        
        if all_ok:
            
            if self.report is not None:
                self.report.add_Text('Updating Point Source Locations!')
            if self.log is not None:
                self.log.info('Updating Point Source Locations!')
            
                            
            for jCCD, CCDk in enumerate(CCDs):
            
                tr = transfs_dict[CCDk]
                
                simmx = strackers[CCDk].get_similaritymx(tr.scale, 
                                 tr.rotation, tr.translation)
                
                strackers[CCDk].apply_patt_transform(simmx)
                
                
        else:
            # raise FLAG
            self.dd.flags.add('STARS_MISSING')
        
        
        
        
        
