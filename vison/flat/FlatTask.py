#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 16:00:10 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
import os

from vison.pipe.task import Task
from vison.datamodel import core
from vison.datamodel import ccd
#from vison.pipe import lib as pilib
from vison.support import context
# END IMPORT


class FlatTask(Task):

    def __init__(self, *args, **kwargs):
        super(FlatTask, self).__init__(*args, **kwargs)

    def check_data(self):
        """ """
        test = self.inputs['test']
        if 'FLAT01'in test:  # AD-HOC modification of test label
            kwargs = dict(figkeys=['FL0Xchecks_offsets', 'FL0Xchecks_deltaoff',
                                   'FL0Xchecks_stds',
                                   'FL0Xchecks_flu',
                                   'FL0Xchecks_imgstd'])
        elif 'FLAT_STB' in test:
            kwargs = dict(figkeys=['FL0Xchecks_offsets', 'FL0Xchecks_deltaoff',
                                   'FL0Xchecks_stds',
                                   'FL0Xchecks_flu',
                                   'FL0Xchecks_imgstd'])
        elif 'FLAT02' in test:
            kwargs = dict(figkeys=['FL0Xchecks_offsets', 'FL0Xchecks_deltaoff',
                                   'FL0Xchecks_stds',
                                   'FL0Xchecks_flu',
                                   'FL0Xchecks_imgstd'])
        elif test == 'PTC01':
            kwargs = dict(figkeys=['PTC0Xchecks_offsets', 'PTC0Xchecks_deltaoff',
                                   'PTC0Xchecks_stds',
                                   'PTC0Xchecks_flu', 'PTC0Xchecks_imgstd'])
        elif 'FLATFLUX00' in test:
            kwargs = dict(figkeys=['PTC0Xchecks_offsets', 'PTC0Xchecks_deltaoff',
                                   'PTC0Xchecks_stds',
                                   'PTC0Xchecks_flu', 'PTC0Xchecks_imgstd'])
        elif 'PTC02' in test:
            kwargs = dict(figkeys=['PTC0Xchecks_offsets', 'PTC0Xchecks_deltaoff',
                                   'PTC0Xchecks_stds',
                                   'PTC0Xchecks_flu', 'PTC0Xchecks_imgstd'])
        elif test in ['NL01', 'NL02'] or 'NL02' in test:
            kwargs = dict(figkeys=['NL01checks_offsets', 'NL01checks_deltaoff',
                                   'NL01checks_stds',
                                   'NL01checks_flu',
                                   'NL01checks_imgstd'])

        Task.check_data(self, **kwargs)

    def get_checkstats_ST(self, **kwargs):
        """ """

        # Initialize new columns

        Xindices = copy.deepcopy(self.dd.indices)

        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad', vals=context.Quads))

        valini = 0.

        newcolnames_off = ['offset_pre', 'offset_ove']
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off, Xindices,
                               dtype='float32', valini=valini)

        newcolnames_deltaoff = ['deltaoff_pre', 'deltaoff_ove']
        for newcolname_deltaoff in newcolnames_deltaoff:
            self.dd.initColumn(newcolname_deltaoff, Xindices,
                               dtype='float32', valini=valini)

        self.dd.initColumn('flu_med_img', Xindices,
                           dtype='float32', valini=valini)
        self.dd.initColumn('flu_std_img', Xindices,
                           dtype='float32', valini=valini)

        newcolnames_std = ['std_pre', 'std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std, Xindices,
                               dtype='float32', valini=valini)

        nObs = Xindices.shape[0]
        CCDs = Xindices.get_vals('CCD')
        Quads = Xindices.get_vals('Quad')

        # Get statistics in different regions

        trimscans = dict(pre=[25, 5],
                         img=[5, 5],
                         ove=[5, 5])

        if not self.drill:
            # if 3==2: # DEBUG!

            for iObs in range(nObs):
                # for iObs in range(3): # TESTS
                #    print 'WARNING: reduced count to just 3 frames for TESTS!'

                print('get_checkstats_ST: Getting metrics from frame %i/%i' % (iObs + 1, nObs,))

                for jCCD, CCDk in enumerate(CCDs):
                    dpath = self.dd.mx['datapath'][iObs, jCCD]
                    ffits = os.path.join(dpath, '%s.fits' %
                                         self.dd.mx['File_name'][iObs, jCCD])
                    ccdobj = ccd.CCD(ffits)
                    vstart = self.dd.mx['vstart'][iObs][jCCD]
                    vend = self.dd.mx['vend'][iObs][jCCD]

                    for kQ, Quad in enumerate(Quads):

                        for reg in ['pre', 'ove']:
                            stats_bias = ccdobj.get_stats(
                                Quad,
                                sector=reg,
                                statkeys=[
                                    'median',
                                    'std'],
                                trimscan=trimscans[reg],
                                ignore_pover=True,
                                extension=-1,
                                VSTART=vstart,
                                VEND=vend)
                            self.dd.mx['offset_%s' %
                                       reg][iObs, jCCD, kQ] = stats_bias[0]
                            self.dd.mx['std_%s' %
                                       reg][iObs, jCCD, kQ] = stats_bias[1]

                        stats_img = ccdobj.get_stats(
                            Quad,
                            sector='img',
                            statkeys=[
                                'median',
                                'std'],
                            trimscan=trimscans['img'],
                            ignore_pover=True,
                            extension=-1)

                        self.dd.mx['flu_med_img'][iObs,
                                                  jCCD, kQ] = stats_img[0] - stats_bias[0]
                        self.dd.mx['flu_std_img'][iObs,
                                                  jCCD, kQ] = stats_img[1]

            for jCCD, CCDk in enumerate(CCDs):

                for kQ, Quad in enumerate(Quads):

                    for reg in ['pre', 'ove']:

                        _meanoff = np.nanmean(self.dd.mx['offset_%s' % reg][:, jCCD, kQ])

                        for iObs in range(nObs):

                            self.dd.mx['deltaoff_%s' % reg][iObs, jCCD, kQ] = \
                                self.dd.mx['offset_%s' % reg][iObs, jCCD, kQ] - _meanoff

    def check_metrics_ST(self, **kwargs):
        """

        TODO:
            - offset levels (pre and over-scan), abs. and relative
            - RON in pre and overscan
            - fluence in image area [script-column-dependent]
            - variance in image area [script-column-dependent]

        """

        #test = self.inputs['test']

        Xindices = self.dd.indices
        CCDs = Xindices.get_vals('CCD')

        if self.report is not None:
            self.report.add_Section(
                keyword='check_metrics', Title='Offsets, RON, Fluences', level=1)

        # absolute value of offsets

        offsets_lims = self.perflimits['offsets_lims']
        for reg in ['pre', 'ove']:
            arr = self.dd.mx['offset_%s' % reg]
            _compliance_offsets = self.check_stat_perCCDandQ(
                arr, offsets_lims, CCDs)

            self.addComplianceMatrix2Self(_compliance_offsets, 'offset_%s' % reg)

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
            arr = self.dd.mx['offset_%s' % reg][:] - self.dd.mx['offset_pre'][:]
            _xcheck_offsets = self.check_stat_perCCDandQ(arr, _lims, CCDs)

            self.addComplianceMatrix2Self(_xcheck_offsets, 'offsets_grad_%s' % reg)

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

            self.addComplianceMatrix2Self(_compliance_std, 'std_%s' % reg)

            if not self.IsComplianceMatrixOK(_compliance_std):
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('RON_OOL')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)

        if not ('NL02' in self.inputs['test']):

            # IMG FLUENCES
            FLU_lims = self.perflimits['FLU_lims']  # dict

            _compliance_flu = self.check_stat_perCCDandCol(
                self.dd.mx['flu_med_img'], FLU_lims, CCDs)

            self.addComplianceMatrix2Self(_compliance_flu, 'fluence')

            # IMG FLUXES

            fluences = self.dd.mx['flu_med_img'][:].copy()
            exptime = self.dd.mx['exptime'][:].copy()

            ixnozero = np.where(exptime[:, 0] > 0)

            _f = np.squeeze(fluences[ixnozero, ...])

            _f[np.where(_f >= 0.8 * 2**16)] = np.nan  # mask-out saturations
            _e = np.expand_dims(np.squeeze(exptime[ixnozero, ...]), axis=-1)

            fluxes = np.nanmean(_f / _e, axis=0)  # CRUDE!
            sat_times = 2.**16 / fluxes
            sat_times = np.expand_dims(sat_times, axis=0)

            exp_sat_time = self.ogse.profile['tFWC_flat']['nm%i' % self.inputs['wavelength']]
            sat_time_lims = dict()
            for CCD in CCDs:
                sat_time_lims[CCD] = (exp_sat_time * np.array([0.9, 1.1])).tolist()

            _compliance_flux = self.check_stat_perCCDandQ(
                sat_times, sat_time_lims, CCDs)

            self.addComplianceMatrix2Self(_compliance_flux, 'flux')

            if not self.IsComplianceMatrixOK(_compliance_flux):
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('FLUX_OOL')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_flux, label='COMPLIANCE SATURATION TIME (FLUX)')
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_flux, label='COMPLIANCE SATURATION TIME (FLUX)',
                    caption='Saturation times in seconds.')

            if not self.IsComplianceMatrixOK(_compliance_flu):
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('FLUENCE_OOL')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_flu, label='COMPLIANCE FLUENCE:')
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_flu, label='COMPLIANCE FLUENCE:')
