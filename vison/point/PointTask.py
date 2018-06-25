#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:55:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import copy
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.datamodel import compliance as complimod
from vison.pipe.task import Task
from vison.point import lib as polib
from vison.point import startracker as strackermod
from vison.datamodel import core, ccd
#from vison.pipe import lib as pilib
from vison.support import context
# END IMPORT

BGD_lims = OrderedDict(CCD1=OrderedDict(E=[-5., 10.]))
for Q in ['F', 'G', 'H']:
    BGD_lims['CCD1'][Q] = copy.deepcopy(BGD_lims['CCD1']['E'])
for iCCD in [2, 3]:
    BGD_lims['CCD%i' % iCCD] = copy.deepcopy(BGD_lims['CCD1'])


class PointTask(Task):

    stampw = polib.stampw

    def __init__(self, *args, **kwargs):
        super(PointTask, self).__init__(*args, **kwargs)

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
        if 'Spot' not in Sindices.names:
            Sindices.append(core.vIndex(
                'Spot', vals=strackermod.starnames))

        self.dd.initColumn('chk_x', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('chk_y', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('chk_peak', Sindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('chk_fluence', Sindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('chk_fwhmx', Sindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('chk_fwhmy', Sindices,
                           dtype='float32', valini=valini)

        chkkeycorr = dict(chk_x='x', chk_y='y', chk_peak='peak', chk_fluence='fluence',
                          chk_fwhmx='fwhmx', chk_fwhmy='fwhmy')

        nObs, _, _ = Qindices.shape
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

        if not self.IsComplianceMatrixOK(_compliance_flu):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCE_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_flu, label='COMPLIANCE FLUENCE:')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_flu, label='COMPLIANCE FLUENCE:')
