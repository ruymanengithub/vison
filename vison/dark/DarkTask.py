#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:54:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
import os

from vison.datamodel import core, ccd
#from vison.pipe import lib as pilib
from vison.pipe.task import Task
from vison.support import context
# END IMPORT


class DarkTask(Task):

    from vison.dark.darkaux import get_DarkDefectsMask_CDP

    def __init__(self, *args, **kwargs):
        super(DarkTask, self).__init__(*args, **kwargs)

    def check_data(self):
        """ """
        test = self.inputs['test']
        if test == 'BIAS01':
            kwargs = dict(figkeys=['B01checks_offsets', 'B01checks_stds'])
        elif test == 'DARK01':
            kwargs = dict(figkeys=['D01checks_offsets', 'D01checks_stds',
                                   'D01checks_flu'])
        Task.check_data(self, **kwargs)

    def get_checkstats_ST(self, **kwargs):
        """ """

        test = self.inputs['test']

        # Initialize new columns

        Xindices = copy.deepcopy(self.dd.indices)

        if 'Quad' not in Xindices.names:
            Xindices.append(core.vIndex('Quad', vals=context.Quads))

        valini = 0.

        newcolnames_off = ['offset_pre', 'offset_img', 'offset_ove']
        for newcolname_off in newcolnames_off:
            self.dd.initColumn(newcolname_off, Xindices,
                               dtype='float32', valini=valini)

        newcolnames_std = ['std_pre', 'std_img', 'std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std, Xindices,
                               dtype='float32', valini=valini)

        if test == 'DARK01':
            self.dd.initColumn('chk_flu_img', Xindices,
                               dtype='float32', valini=valini)

        nObs, _, _ = Xindices.shape
        CCDs = Xindices.get_vals('CCD')
        Quads = Xindices.get_vals('Quad')

        # Get statistics in different regions

        if not self.drill:

            for iObs in range(nObs):
                for jCCD, CCDk in enumerate(CCDs):
                    dpath = self.dd.mx['datapath'][iObs, jCCD]
                    ffits = os.path.join(dpath, '%s.fits' %
                                         self.dd.mx['File_name'][iObs, jCCD])
                    ccdobj = ccd.CCD(ffits)
                    vstart = self.dd.mx['vstart'][iObs][jCCD]
                    vend = self.dd.mx['vend'][iObs][jCCD]

                    for kQ, Quad in enumerate(Quads):

                        for reg in ['pre', 'ove', 'img']:
                            stats = ccdobj.get_stats(Quad, sector=reg, statkeys=['median', 'std'], trimscan=[5, 5],
                                                     ignore_pover=True, extension=-1, VSTART=vstart, VEND=vend)
                            self.dd.mx['offset_%s' %
                                       reg][iObs, jCCD, kQ] = stats[0]
                            self.dd.mx['std_%s' %
                                       reg][iObs, jCCD, kQ] = stats[1]

                            if test == 'DARK01' and reg == 'img':
                                offset_cbe = np.mean([self.dd.mx['offset_pre'][iObs, jCCD, kQ],
                                                      self.dd.mx['offset_ove'][iObs, jCCD, kQ]])
                                self.dd.mx['chk_flu_%s' % reg][iObs,
                                                               jCCD, kQ] = stats[0] - offset_cbe

    def check_metrics_ST(self, **kwargs):
        """ 

        """

        test = self.inputs['test']

        Xindices = self.dd.indices
        CCDs = Xindices.get_vals('CCD')

        if self.report is not None:
            self.report.add_Section(
                keyword='check_ronoffset', Title='Offsets and RON', level=1)

        # absolute value of offsets

        offsets_lims = self.perflimits['offsets_lims']
        
        
        if test == 'BIAS01':
            regs_off = ['pre', 'img', 'ove']
        elif test == 'DARK01':
            regs_off = ['pre', 'ove']

        for reg in regs_off:
            arr = self.dd.mx['offset_%s' % reg]
            _compliance_offsets = self.check_stat_perCCDandQ(
                arr, offsets_lims, CCDs)

            #if not self.IsComplianceMatrixOK(_compliance_offsets):
            if not _compliance_offsets.IsCompliant():
                self.dd.flags.add('POORQUALDATA')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)

        # cross-check of offsets: referred to pre-scan

        offsets_gradients = self.perflimits['offsets_gradients']

        if test == 'BIAS01':
            regs_grad = ['img', 'ove']
        elif test == 'DARK01':
            regs_grad = ['ove']

        for ireg, reg in enumerate(regs_grad):
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

        if test == 'BIAS01':
            regs_std = ['pre', 'img', 'ove']
        elif test == 'DARK01':
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

        # Fluence (only DARK01)

        if test == 'DARK01':

            Flu_lims = self.perflimits['Flu_lims']

            arr = self.dd.mx['chk_flu_img'][:].copy()
            _compliance_Flu = self.check_stat_perCCDandQ(arr, Flu_lims, CCDs)

            if not self.IsComplianceMatrixOK(_compliance_Flu):
                self.dd.flags.add('POORQUALDATA')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_Flu, label='COMPLIANCE FLUENCE [img]:')
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_Flu, label='COMPLIANCE FLUENCE [img]:')
