#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:56:00 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import copy
import os
from collections import OrderedDict
import pandas as pd
import sys

#from vison.support import context
from vison.inject import lib as ilib
from vison.pipe.task import Task
from vison.datamodel import core, ccd
#from vison.pipe import lib as pilib
from vison.support import context
#from vison.pipe.task import Task
from vison.datamodel import cdp
from vison.support import utils, files
# END IMPORT

lineoffsets = ilib.lineoffsets


def _get_CCDhalf(Q):
    if Q in ['E', 'F']:
        return 'B'
    elif Q in ['G', 'H']:
        return 'T'


class InjTask(Task):

    def __init__(self, *args, **kwargs):
        super(InjTask, self).__init__(*args, **kwargs)

    def check_data(self, **kwargs):
        """ """
        test = self.inputs['test']
        if test == 'CHINJ01':
            _kwargs = dict(figkeys=['CH01checks_offsets', 'CH01checks_deltaoff',
                                    'CH01checks_stds',
                                    'CH01checks_injlevel', 'CH01checks_injstd'])
            kwargs.update(_kwargs)
        elif test == 'CHINJ02':
            _kwargs = dict(figkeys=['CH02checks_offsets', 'CH02checks_deltaoff',
                                    'CH02checks_stds',
                                    'CH02checks_injlevel', 'CH02checks_injstd'])
            kwargs.update(_kwargs)

        Task.check_data(self, **kwargs)

    def prepare_images(self, doExtract=True, doBadPixels=True, doMask=True, doOffset=True,
                       doBias=True, doFF=False):
        """

        InjTask: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            [bias structure subtraction, if available]
            cosmetics masking

        """
        super(InjTask, self).prepare_images(
            doExtract=doExtract, doBadPixels=doBadPixels,
            doMask=doMask, doOffset=doOffset, doBias=doBias,
            doFF=doFF)

    def predict_expected_injlevels(self, teststruct):
        """ """
        CCDs = ['CCD1', 'CCD2', 'CCD3']
        Quads = ['E', 'F', 'G', 'H']

        SecTags = dict(E='B', F='B', G='T', H='T')

        Ncols = teststruct['Ncols']

        expectation = OrderedDict()

        for jCCD, CCDkey in enumerate(CCDs):

            expectation[CCDkey] = OrderedDict()

            for Q in Quads:
                expectation[CCDkey][Q] = OrderedDict()
                sectag = SecTags[Q]

                for icol in range(1, Ncols + 1):
                    coldict = teststruct['col%03i' % icol]
                    chinj = coldict['chinj']
                    IGs = (coldict['IG1_%i_%s' % (jCCD + 1, sectag)],
                           coldict['IG2_%s' % sectag])
                    ID = (coldict['IDL'], coldict['IDH'])
                    id_timing = (coldict['id_wid'], coldict['id_dly'])
                    toi_ch = coldict['toi_ch']

                    if chinj == 0:
                        _inj = 0.
                    else:
                        _inj = ilib.predict_inj_level(ID, IGs, id_timing,
                                                      toi_ch, sectag)

                    expectation[CCDkey][Q]['col%03i' % icol] = _inj

        return expectation

    def get_FluenceAndGradient_limits(self):
        """ """

        tmpstructure = self.build_scriptdict(diffvalues={}, elvis=self.elvis)
        inj_exp = self.predict_expected_injlevels(tmpstructure)

        Flu_lims = OrderedDict()
        FluGrad_lims = OrderedDict()

        for CCDkey in list(inj_exp.keys()):
            Flu_lims[CCDkey] = OrderedDict()
            FluGrad_lims[CCDkey] = OrderedDict()

            for Q in list(inj_exp[CCDkey].keys()):
                Flu_lims[CCDkey][Q] = OrderedDict()
                FluGrad_lims[CCDkey][Q] = OrderedDict()

                for colkey in list(inj_exp[CCDkey][Q].keys()):
                    _inj = inj_exp[CCDkey][Q][colkey]

                    if np.isnan(_inj):
                        Flu_lims[CCDkey][Q][colkey] = [-10., 1.01 * 2.**16]
                        FluGrad_lims[CCDkey][Q][colkey] = [0., 1.E4]
                    else:
                        Flu_lims[CCDkey][Q][colkey] = _inj * \
                            (1. + np.array([-0.5, 0.5]))
                        FluGrad_lims[CCDkey][Q][colkey] = _inj * \
                            0.3 * (1. + np.array([-0.9, 0.5]))

        return Flu_lims, FluGrad_lims

    def get_checkstats_ST(self, **kwargs):
        """ """

        #test = self.inputs['test']

        if 'pattern' in kwargs:
            pattern = kwargs['pattern']

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

        newcolnames_std = ['std_pre', 'std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std, Xindices,
                               dtype='float32', valini=valini)

        for chk_inj_col in ['chk_mea_inject', 'chk_med_inject', 'chk_std_inject']:
            self.dd.initColumn(chk_inj_col, Xindices,
                               dtype='float32', valini=valini)

        nObs, _, _ = Xindices.shape
        CCDs = Xindices.get_vals('CCD')
        Quads = Xindices.get_vals('Quad')

        # Get statistics in different regions

        trimscans = dict(pre=[25, 5],
                         img=[5, 5],
                         ove=[5, 5])

        if not self.drill:

            for iObs in range(nObs):

                if self.debug:
                    print('InjTask.get_checkstats_ST: processing ObsID %i/%i' % (iObs + 1, nObs))

                for jCCD, CCDk in enumerate(CCDs):
                    dpath = self.dd.mx['datapath'][iObs, jCCD]
                    ffits = os.path.join(dpath, '%s.fits' %
                                         self.dd.mx['File_name'][iObs, jCCD])
                    ccdobj = ccd.CCD(ffits)

                    vstart = self.dd.mx['vstart'][iObs][jCCD]
                    vend = self.dd.mx['vend'][iObs][jCCD]
                    dochinj = self.dd.mx['chinj'][iObs][jCCD]
                    if 'pattern' not in locals():
                        non = self.dd.mx['chinj_on'][iObs][jCCD]
                        noff = self.dd.mx['chinj_of'][iObs][jCCD]
                        nrep = (vend - vstart) // (non + noff) + 1
                        pattern = (non, noff, nrep)

                    for kQ, Quad in enumerate(Quads):

                        for reg in ['pre', 'ove']:
                            stats = ccdobj.get_stats(
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
                                       reg][iObs, jCCD, kQ] = stats[0]
                            self.dd.mx['std_%s' %
                                       reg][iObs, jCCD, kQ] = stats[1]

                        if dochinj:

                            ccdobj.sub_offset(
                                Quad,
                                method='row',
                                scan='pre',
                                trimscan=trimscans['pre'],
                                ignore_pover=True,
                                extension=-1)

                            extract_res = ilib.extract_injection_lines(ccdobj,
                                                                       Quad, pattern, VSTART=vstart,
                                                                       VEND=vend, suboffmean=False)

                            istats = extract_res['stats_injection']

                            self.dd.mx['chk_mea_inject'][iObs, jCCD,
                                                         kQ] = istats['mean']
                            self.dd.mx['chk_med_inject'][iObs, jCCD,
                                                         kQ] = istats['p50']
                            self.dd.mx['chk_std_inject'][iObs, jCCD,
                                                         kQ] = istats['std']

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
            - mean fluence/signal in image area [script-column-dependent]
            - med fluence/signal in image area [script-column-dependent]
            - std in image area [script-column-dependent]


        """
        # test = self.inputs['test']

        Xindices = self.dd.indices
        CCDs = Xindices.get_vals('CCD')

        if self.report is not None:
            self.report.add_Section(
                keyword='check_metrics', Title='Offsets, RON, Injection levels', level=1)

        # absolute value of offsets

        offsets_lims = self.perflimits['offsets_lims']
        for reg in ['pre', 'ove']:
            arr = self.dd.mx['offset_%s' % reg]
            _compliance_offsets = self.check_stat_perCCDandQ(
                arr, offsets_lims, CCDs)

            self.addComplianceMatrix2Self(_compliance_offsets, 'offsets_%s' % reg)

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

        # IMG Signal Levels

        Flu_lims = self.perflimits['Flu_lims']  # dict

        _compliance_flu = self.check_stat_perCCDQandCol(
            self.dd.mx['chk_med_inject'], Flu_lims, CCDs)

        self.addComplianceMatrix2Self(_compliance_flu, 'fluence')

        if not self.IsComplianceMatrixOK(_compliance_flu):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCE_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_flu, label='COMPLIANCE FLUENCE:')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_flu, label='COMPLIANCE FLUENCE:',
                caption='Beware the expected fluences are computed rather coarsely.' +
                ' Use with caution, and also check plots of injection vs. time/CCD voltage(s) ' +
                'to judge quality of data.')

        # IMG std (injection-noise) levels

        FluGrad_lims = self.perflimits['FluGrad_lims']  # dict

        _compliance_flugrad = self.check_stat_perCCDQandCol(
            self.dd.mx['chk_std_inject'], FluGrad_lims, CCDs)

        if not self.IsComplianceMatrixOK(_compliance_flugrad):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('FLUENCEGRAD_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_flugrad, label='COMPLIANCE FLUENCE GRADIENT:')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_flugrad, label='COMPLIANCE FLUENCE GRADIENT:',
                caption='Beware the expected fluence gradients are computed rather coarsely.' +
                ' Use with caution, and also check plots of injection STD vs. time/CCD voltage(s) ' +
                'to judge quality of data.')

    def basic_analysis(self):
        """

        Basic analysis of data.

        **METACODE**

        ::

            f. e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        extract average 2D injection pattern (and save)
                        produce average profile along/across lines
                        measure charge-inj. non-uniformity
                        measure charge spillover into non-injection
                        measure stats of injection (mean, med, std, min/max, percentiles)

            plot average inj. profiles along lines f. each CCD, Q and VOLTAGE
                save as a rationalized set of curves
            plot average inj. profiles across lines f. each CCD, Q and VOLTAGE
                save as a rationalized set of  curves

            Report injection stats as a table/tables

        """
        import warnings
        warnings.warn('WARNING! task InjTask.basic_analysis has not been properly tested!')
        #sys.exit('BROKEN! TODO: address polymorphism correctly, for figures, cdps, etc.')

        testname = self.inputs['test']
        testkey = dict(CHINJ01='CH01',
                       CHINJ02='CH02')[testname]

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='%s Extraction' % testname,
                level=0)

        DDindices = copy.deepcopy(self.dd.indices)

        nObs, nCCD, nQuad = DDindices.shape[0:3]
        Quads = DDindices.get_vals('Quad')
        CCDs = DDindices.get_vals('CCD')

        ccdpicklespath = self.inputs['subpaths']['ccdpickles']
        prodspath = self.inputs['subpaths']['products']

        # Initializing new columns

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        valini = 0.

        # measure charge-inj. non-uniformity
        # measure charge spillover into non-injection
        # measure stats of injection (mean, med, std, min/max, percentiles)

        self.dd.initColumn('chinj_nonuni', DDindices, dtype='float32', valini=valini)
        self.dd.initColumn('chinj_spill', DDindices, dtype='float32', valini=valini)
        statkeys = ['mean', 'std', 'min', 'max', 'p5', 'p25', 'p50', 'p75',
                    'p95']
        for statkey in statkeys:
            self.dd.initColumn('chinj_%s' % statkey, DDindices, dtype='float32', valini=valini)

        # EXTRACTION TABLE

        NP = nObs * nCCD * nQuad

        # OBSID CCD Q IG1 id_dly MEAN MEDIAN NONUNI

        CH0X_dd = OrderedDict()
        CH0X_dd['ObsID'] = np.zeros(NP, dtype='int32')
        CH0X_dd['CCD'] = np.zeros(NP, dtype='int32')
        CH0X_dd['Q'] = np.zeros(NP, dtype='int32')
        if testname == 'CHINJ01':
            CH0X_dd['IG1'] = np.zeros(NP, dtype='float32')
        elif testname == 'CHINJ02':
            CH0X_dd['IDL'] = np.zeros(NP, dtype='float32')
        CH0X_dd['ID_DLY'] = np.zeros(NP, dtype='float32')
        CH0X_dd['MEAN_INJ'] = np.zeros(NP, dtype='float32')
        CH0X_dd['MED_INJ'] = np.zeros(NP, dtype='float32')
        CH0X_dd['NU_INJ'] = np.zeros(NP, dtype='float32')

        # Initializing injection profiles

        prof_alrow_cdp = cdp.CDP()
        prof_alrow_cdp.rootname = 'PROFS_ALROW_%s' % testname
        prof_alrow_cdp.header = CDP_header.copy()
        prof_alrow_cdp.path = prodspath
        prof_alrow_cdp.data = OrderedDict()

        prof_alcol_cdp = cdp.CDP()
        prof_alcol_cdp.rootname = 'PROFS_ALCOL_%s' % testname
        prof_alcol_cdp.header = CDP_header.copy()
        prof_alcol_cdp.path = prodspath
        prof_alcol_cdp.data = OrderedDict()

        xdummy = np.arange(10, dtype='float32')
        ydummy = np.zeros(10, dtype='float32')
        for jCCD, CCDk in enumerate(CCDs):
            prof_alrow_cdp.data[CCDk] = OrderedDict()
            prof_alcol_cdp.data[CCDk] = OrderedDict()

            for kQ, Q in enumerate(Quads):
                prof_alrow_cdp.data[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                           y=OrderedDict())
                prof_alcol_cdp.data[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                           y=OrderedDict())

                for iObs in range(nObs):

                    if testname == 'CHINJ01':
                        IG1_key = 'IG1_%i_%s' % (jCCD + 1, _get_CCDhalf(Q))
                        IG1_val = self.dd.mx[IG1_key][iObs, jCCD]
                        sub_tag = 'IG1_%.2fV' % IG1_val
                    elif testname == 'CHINJ02':
                        IDL_val = self.dd.mx['IDL'][iObs, jCCD]
                        sub_tag = 'IDL_%.2fV' % IDL_val

                    prof_alrow_cdp.data[CCDk][Q]['x'][sub_tag] = xdummy.copy()
                    prof_alrow_cdp.data[CCDk][Q]['y'][sub_tag] = ydummy.copy()
                    prof_alcol_cdp.data[CCDk][Q]['x'][sub_tag] = xdummy.copy()
                    prof_alcol_cdp.data[CCDk][Q]['y'][sub_tag] = ydummy.copy()

        # Images with average profiles of injection stacked

        Profile_Images = OrderedDict()

        Qshape = self.ccdcalc.wQ, self.ccdcalc.hQ

        for jCCD, CCDk in enumerate(CCDs):
            Profile_Images[CCDk] = OrderedDict()
            for kQ, Q in enumerate(Quads):
                Profile_Images[CCDk][Q] = np.zeros(Qshape, dtype='float32')

        # The hardwork

        prof_max_vals = []

        if not self.drill:

            _Quads_dict = dict(bottom=['E', 'F'],
                               top=['G', 'H'])

            for iObs in range(nObs):
                #            for iObs in range(5): # TESTS

                ObsID = self.dd.mx['ObsID'][iObs]

                print('Processing Obsid %i/%i' % (iObs + 1, nObs))

                for jCCD, CCDk in enumerate(CCDs):

                    ccdpickf = os.path.join(ccdpicklespath,
                                            '%s.pick' % self.dd.mx['ccdobj_name'][iObs, jCCD])

                    ccdobj = files.cPickleRead(ccdpickf)

                    vstart = self.dd.mx['vstart'][iObs, jCCD]
                    vend = self.dd.mx['vend'][iObs, jCCD]
                    dochinj = self.dd.mx['chinj'][iObs, jCCD]
                    non = self.dd.mx['chinj_on'][iObs, jCCD]
                    noff = self.dd.mx['chinj_of'][iObs, jCCD]
                    nrep = (vend - vstart) // (non + noff) + 1
                    pattern = (non, noff, nrep)

                    for kQ, Q in enumerate(Quads):

                        ix = iObs * nCCD * nQuad + jCCD * nQuad + kQ

                        if dochinj:

                            if testname == 'CHINJ01':
                                IG1_key = 'IG1_%i_%s' % (jCCD + 1, _get_CCDhalf(Q))
                                IG1_val = self.dd.mx[IG1_key][iObs, jCCD]
                                sub_tag = 'IG1_%.2fV' % IG1_val
                            elif testname == 'CHINJ02':
                                IDL_val = self.dd.mx['IDL'][iObs, jCCD]
                                sub_tag = 'IDL_%.2fV' % IDL_val

                            id_dly = self.dd.mx['id_dly'][iObs, jCCD]
                            toi_ch = self.dd.mx['toi_ch'][iObs, jCCD]

                            if np.isclose(id_dly / toi_ch, 1.5):
                                CCDhalf = 'top'
                            elif np.isclose(id_dly / toi_ch, 2.5):
                                CCDhalf = 'bottom'

                            ext_res = ilib.extract_injection_lines(
                                ccdobj, Q, pattern, VSTART=vstart, VEND=vend, suboffmean=False)

                            stats = ext_res['stats_injection']

                            for skey in statkeys:
                                self.dd.mx['chinj_%s' % skey][iObs, jCCD, kQ] = \
                                    stats[skey]

                            nonuni = (stats['p95'] - stats['p5']) / stats['p50']

                            self.dd.mx['chinj_nonuni'][iObs, jCCD, kQ] = nonuni

                            spill = ilib.get_spill(ext_res['avprof_alcol'], pattern)

                            self.dd.mx['chinj_spill'][iObs, jCCD, kQ] = spill

                            yalrows = ext_res['avprof_alrow'].copy()
                            yalcols = ext_res['avprof_alcol'].copy()

                            if Q in _Quads_dict[CCDhalf]:

                                prof_alrow_cdp.data[CCDk][Q]['y'][sub_tag] = yalrows.copy()
                                prof_alrow_cdp.data[CCDk][Q]['x'][sub_tag] = \
                                    np.arange(len(yalrows), dtype='float32')

                                prof_alcol_cdp.data[CCDk][Q]['y'][sub_tag] = yalcols.copy()
                                prof_alcol_cdp.data[CCDk][Q]['x'][sub_tag] = \
                                    np.arange(len(yalcols), dtype='float32')

                                Profile_Images[CCDk][Q][self.ccdcalc.prescan: - \
                                    self.ccdcalc.overscan, iObs] = yalrows.copy()

                                prof_max_vals.append(np.nanmax([yalrows.max(), yalcols.max()]))

                            CH0X_dd['ObsID'][ix] = ObsID
                            CH0X_dd['CCD'][ix] = jCCD
                            CH0X_dd['Q'][ix] = kQ
                            if testname == 'CHINJ01':
                                CH0X_dd['IG1'][ix] = IG1_val
                            elif testname == 'CHINJ02':
                                CH0X_dd['IDL'][ix] = IDL_val
                            CH0X_dd['ID_DLY'][ix] = id_dly
                            CH0X_dd['MEAN_INJ'][ix] = self.dd.mx['chinj_mean'][iObs, jCCD, kQ]
                            CH0X_dd['MED_INJ'][ix] = self.dd.mx['chinj_p50'][iObs, jCCD, kQ]
                            CH0X_dd['NU_INJ'][ix] = self.dd.mx['chinj_nonuni'][iObs, jCCD, kQ]

                            # if CCDk == 'CCD2' and Q == 'E':
                            #    stop()

        for jCCD, CCDk in enumerate(CCDs):

            tmpccdobj = ccd.CCD(withpover=True)
            tmpccdobj.add_extension(data=np.zeros(tmpccdobj.shape, dtype='float32'))

            for Q in Quads:
                tmpccdobj.set_quad(Profile_Images[CCDk][Q].copy(), Q,
                                   canonical=True, extension=-1)

            profimg = tmpccdobj.extensions[-1].data.copy()

            jprofimgcdp = cdp.CCD_CDP(ID=self.ID,
                                      BLOCKID=self.BLOCKID,
                                      CHAMBER=self.CHAMBER)

            profimgdata = dict(PROFILES=profimg.copy(),
                               labels=['PROFILES'])

            jprofimgcdp.ingest_inputs(data=profimgdata,
                                      meta=None,
                                      header=CDP_header)

            jprofimgcdp.path = prodspath
            jprofimgcdp.rootname = 'EUC_INJPROFIMG_%s_%s' % \
                (CCDk, self.inputs['BLOCKID'])
            self.save_CDP(jprofimgcdp)
            self.pack_CDP_to_dd(jprofimgcdp, 'INJPROFIMG_%s' % CCDk)

        # plot average inj. profiles along/across lines
        # save as a rationalized set of curves

        #maxmedinjection = np.nanmax(self.dd.mx['chinj_p50'][:])
        maxmedinjection = np.nanmax(prof_max_vals)

        prof_alrow_cdp.data['labelkeys'] = list(prof_alrow_cdp.data[CCDs[0]][Quads[0]]['x'].keys())

        fdict_alrow = self.figdict['%s_alrow' % testkey][1]
        fdict_alrow['data'] = prof_alrow_cdp.data.copy()
        fdict_alrow['meta']['ylim'] = [0., maxmedinjection * 1.1]

        prof_alcol_cdp.data['labelkeys'] = list(prof_alcol_cdp.data[CCDs[0]][Quads[0]]['x'].keys())

        fdict_alcol = self.figdict['%s_alcol' % testkey][1]
        fdict_alcol['data'] = prof_alcol_cdp.data.copy()
        fdict_alcol['meta']['ylim'] = [0., maxmedinjection * 1.1]

        self.save_CDP(prof_alrow_cdp)
        self.pack_CDP_to_dd(prof_alrow_cdp, 'PROFS_ALROW')
        self.save_CDP(prof_alcol_cdp)
        self.pack_CDP_to_dd(prof_alcol_cdp, 'PROFS_ALCOL')

        if self.report is not None:
            self.addFigures_ST(figkeys=['%s_alrow' % testkey,
                                        '%s_alcol' % testkey],
                               dobuilddata=False)

        # Report injection stats as a table/tables

        # OBSID CCD Q IG1 id_dly MEAN MEDIAN NONUNI

        EXT_dddf = OrderedDict(EXTRACT=pd.DataFrame.from_dict(CH0X_dd))
        EXT_cdp = self.CDP_lib['EXTRACT']
        EXT_cdp.path = prodspath
        EXT_cdp.ingest_inputs(
            data=EXT_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy())

        EXT_cdp.init_wb_and_fillAll(header_title='%s: EXTRACTION' % testname)
        self.save_CDP(EXT_cdp)
        self.pack_CDP_to_dd(EXT_cdp, 'EXTRACT_CDP')

        if self.report is not None:
            def fi(x): return '%i' % x

            def fccd(x): return CCDs[x]

            def fq(x): return Quads[x]

            def ff(x): return '%.2f' % x

            ext_formatters = [fi, fccd, fq, ff, ff, ff, ff, ff]

            caption = '%s: EXTRACTION TABLE' % testname
            Etex = EXT_cdp.get_textable(sheet='EXTRACT', caption=caption,
                                        fitwidth=True,
                                        tiny=True,
                                        formatters=ext_formatters)

            self.report.add_Text(Etex)
