#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: PERSIST01

CCD Persistence test

Created on Tue Aug 29 17:39:00 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import os
import numpy as np
import copy

from vison.datamodel import core
from vison.datamodel import ccd
from vison.pipe.task import HKKeys
from vison.support import context
from vison.ogse import ogse as ogsemod
from vison.datamodel import scriptic as sc
from vison.pipe.task import Task
from vison.datamodel import inputs
from vison.other import PERSIST01aux as P01aux
# END IMPORT

#HKKeys = []

PER01_commvalues = dict(program='CALCAMP', test='PERSIST01',
                        flushes=7, siflsh=1, siflsh_p=500,
                        inisweep=1,
                        vstart=0, vend=2086,
                        toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                        chinj=0,
                        s_tpump=0,
                        v_tpump=0,
                        exptime=0., shuttr=1, e_shuttr=0,
                        wave=2,
                        mirr_on=1,
                        motr_on=0,
                        source='point',
                        comments='')

BGD_lims = OrderedDict(CCD1=OrderedDict())
BGD_lims['CCD1'] = [0.,200.] # ADU
for i in [2, 3]:
    BGD_lims['CCD%i' % i] = copy.deepcopy(BGD_lims['CCD1'])

SATUR_lims = OrderedDict(CCD1=OrderedDict())
SATUR_lims['CCD1']['E'] = [0.002,0.2]
for Q in ['F','G','H']:
    SATUR_lims['CCD1'][Q] = copy.deepcopy(SATUR_lims['CCD1']['E']) 
for i in [2, 3]:
    SATUR_lims['CCD%i' % i] = copy.deepcopy(SATUR_lims['CCD1'])


class PERSIST01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptSATUR', ([float], 'Exposure times to produce latent.')),
        ('exptLATEN', ([float], 'Exposure times to quantify latent.')),
    ])))


class PERSIST01(Task):
    """ """

    inputsclass = PERSIST01_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(PERSIST01, self).__init__(inputs, log, drill, debug)
        self.name = 'PERSIST01'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = P01aux.get_P01figs()
        self.inputs['subpaths'] = dict(figs='figs')
        

    def set_inpdefaults(self, **kwargs):
        """ """
        wavelength = 590.
        self.inpdefaults = dict(
            exptSATUR=100.*self.ogse.profile['tFWC_point']['nm%i' % wavelength],
            exptLATEN=565. # s
        )
        
    def set_perfdefaults(self, **kwargs):
        super(PERSIST01, self).set_perfdefaults(**kwargs)
        
        
        self.perfdefaults['BGD_lims'] = BGD_lims.copy()
        self.perfdefaults['SATUR'] = SATUR_lims.copy()
        

    

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ 
        Builds PERSISTENCE01 script structure dictionary.

        :param exptSATUR: int, saturation exposure time.
        :param exptLATEN: int, latency exposure time.
        :param diffvalues: dict, opt, differential values.


        """
        exptSATUR = self.inputs['exptSATUR']
        exptLATEN = self.inputs['exptLATEN']
        wave = 2 
        mirr_pos = self.ogse.profile['mirror_nom']['F%i' % wave]-10.

        PER01_sdict = dict(col001=dict(frames=5, exptime=0, shuttr=0, 
                                     wave=4,mirr_pos=self.ogse.profile['mirror_nom']['F4'],
                                     source='flat',
                                     comments='RESET'),
                           col002=dict(frames=3, exptime=exptLATEN, shuttr=0, 
                                     wave=4,mirr_pos=self.ogse.profile['mirror_nom']['F4'],
                                     source='flat',
                                     comments='REFER.'),
                           col003=dict(frames=1, exptime=exptSATUR, 
                                     wave=wave,mirr_pos=mirr_pos,
                                     comments='EXPOSE'),
                           col004=dict(frames=3, exptime=exptLATEN, shuttr=0, 
                                     wave=4,mirr_pos=self.ogse.profile['mirror_nom']['F4'],
                                     source='flat',
                                     comments='LATENT'))

        Ncols = len(PER01_sdict.keys())
        PER01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        PER01_commvalues['mirror_pos'] = self.ogse.profile['mirror_nom']['F%i' % wave]
        commvalues.update(PER01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        PER01_sdict = sc.update_structdict(PER01_sdict, commvalues, diffvalues)

        return PER01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(PERSIST01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                      wavedkeys=wavedkeys)

    def check_data(self):
        """ 

        PERSIST01: Checks quality of ingested data.


        **METACODE**

        ::

            check common HK values are within safe / nominal margins
            check voltages in HK match commanded voltages, within margins

            f.e.ObsID:
                f.e.CCD:
                    f.e.Q.:
                        measure offsets in pre-, over-
                        measure std in pre-, over-
                        measure fluence in apertures around Point Sources

            assess std in pre- (~RON) is within allocated margins
            assess offsets in pre-, and over- are equal, within allocated  margins
            assess fluence is ~expected within apertures (PS) for each frame (pre-satur, satur, post-satur)


            plot point source fluence vs. OBSID, all sources
            [plot std vs. time]

            issue any warnings to log
            issue update to report          


        """
        
        kwargs = dict(figkeys=['P01checks_offsets','P01checks_stds',
                               'P01checks_avgflu'])        
        Task.check_data(self,**kwargs)
    
    def prep_data(self):
        """

        PERSIST01: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            cosmetics masking


        """
        super(PERSIST01, self).prepare_images(
            doExtract=True, doBadPixels=False,
            doMask=True, doOffset=True)
        
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

        newcolnames_std = ['std_pre', 'std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std, Xindices,
                               dtype='float32', valini=valini)

        self.dd.initColumn('chk_avgflu_img', Xindices,
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
                            if reg in ['pre','ove']:
                                self.dd.mx['offset_%s' %
                                           reg][iObs, jCCD, kQ] = stats[0]
                                self.dd.mx['std_%s' %
                                           reg][iObs, jCCD, kQ] = stats[1]
                            
                            
                            offset_cbe = np.mean([self.dd.mx['offset_pre'][iObs, jCCD, kQ],
                                                  self.dd.mx['offset_ove'][iObs, jCCD, kQ]])
                            
                            if reg == 'img':
                                self.dd.mx['chk_avgflu_%s' % reg][iObs,
                                          jCCD, kQ] = stats[0] - offset_cbe

    
    def check_metrics_ST(self, **kwargs):
        """ """

        Xindices = self.dd.indices
        CCDs = Xindices.get_vals('CCD')

        if self.report is not None:
            self.report.add_Section(
                keyword='check_ronoffset', Title='Offsets, RON and Background Levels', level=1)
        
        # absolute value of offsets

        offsets_lims = self.perflimits['offsets_lims']
        
        ixreference = np.where(self.dd.mx['label'][:,0] == 'col001')

        regs_off = ['pre', 'ove']

        for reg in regs_off:
            arr = self.dd.mx['offset_%s' % reg][ixreference[0],...].copy()

            _compliance_offsets = self.check_stat_perCCDandQ(
                arr, offsets_lims, CCDs)
            
            self.addComplianceMatrix2Self(_compliance_offsets,'offsets_%s' % reg)

            # if not self.IsComplianceMatrixOK(_compliance_offsets):
            if not _compliance_offsets.IsCompliant():
                self.dd.flags.add('POORQUALDATA')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg,
                    caption='Only considering the bias frames at beginning of '+\
                    'test in this compliance matrix.')

        # cross-check of offsets: referred to pre-scan

        offsets_gradients = self.perflimits['offsets_gradients']

        regs_grad = ['ove']

        for ireg, reg in enumerate(regs_grad):
            _lims = dict()
            for CCDk in CCDs:
                _lims[CCDk] = offsets_gradients[CCDk][reg]
            arr = self.dd.mx['offset_%s' % reg][ixreference[0],...]-self.dd.mx['offset_pre'][ixreference[0],...]
            _xcheck_offsets = self.check_stat_perCCDandQ(arr, _lims, CCDs)
            
            self.addComplianceMatrix2Self(_xcheck_offsets,'offsets_grad_%s' % reg)

            if not self.IsComplianceMatrixOK(_xcheck_offsets):
                self.dd.flags.add('POORQUALDATA')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _xcheck_offsets, label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _xcheck_offsets, label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg,
                    caption='Only considering the bias frames at beginning of '+\
                    'test in this compliance matrix.')

        # absolute value of std

        regs_std = ['pre', 'ove']

        RONs_lims = self.perflimits['RONs_lims']
        for reg in regs_std:
            _compliance_std = self.check_stat_perCCDandQ(
                self.dd.mx['std_%s' % reg][ixreference[0],...], RONs_lims, CCDs)
            
            self.addComplianceMatrix2Self(_compliance_std,'std_%s' % reg)


            if not self.IsComplianceMatrixOK(_compliance_std):
                self.dd.flags.add('POORQUALDATA')
                self.dd.flags.add('RON_OOL')
            if self.log is not None:
                self.addComplianceMatrix2Log(
                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)
            if self.report is not None:
                self.addComplianceMatrix2Report(
                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg,
                    caption='Only considering the bias frames at beginning of '+\
                    'test in this compliance matrix.')
        
        
        # SATURATIONS
        
        ixsat = np.where(self.dd.mx['label'][:,0] == 'col003')
        
        satur_lims = self.perflimits['SATUR']    
        
        arrSAT = self.dd.mx['chk_NPIXSAT'][ixsat[0],...].copy()
        
        NQactive = self.ccdcalc.NAXIS1 * (self.dd.mx['vend'][ixsat[0],...]-self.dd.mx['vstart'][ixsat[0],...])
        NQactive = np.repeat(NQactive[...,np.newaxis],4,axis=-1) # extend to Quads
        
                      
        FracSat = arrSAT / NQactive
        
        _compliance_sat = self.check_stat_perCCDandQ(FracSat, satur_lims, CCDs)
        
        self.addComplianceMatrix2Self(_compliance_sat,'saturation')
        
        
        if not self.IsComplianceMatrixOK(_compliance_sat):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('NONSATURATION')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_sat, label='COMPLIANCE SATURATION FRACTION')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_sat, label='COMPLIANCE SATURATION FRACTION',
                caption='Fraction (over 1) of CCD image saturated. In this test '+\
                    'saturations have to be extended enough for the test to be useful. '+\
                    'Results obtained on the saturation inducing frame only.'
                    )
        
        # Background levels
        
        BGD_lims = self.perflimits['BGD_lims'] # dict
        
        ixbgd = np.where((self.dd.mx['label'][:,0] == 'col002') |
                        (self.dd.mx['label'][:,0] == 'col004'))
        
        arrbgd = self.dd.mx['chk_avgflu_img'][ixbgd[0],...].copy()

        _compliance_bgd = self.check_stat_perCCDandQ(
                arrbgd, BGD_lims, CCDs)
        
        self.addComplianceMatrix2Self(_compliance_bgd,'bgd')
        
        
        if not self.IsComplianceMatrixOK(_compliance_bgd):
            self.dd.flags.add('POORQUALDATA')
            self.dd.flags.add('BGD_OOL')
        if self.log is not None:
            self.addComplianceMatrix2Log(
                _compliance_bgd, label='COMPLIANCE BACKGROUND LEVELS')
        if self.report is not None:
            self.addComplianceMatrix2Report(
                _compliance_bgd, label='COMPLIANCE BACKGROUND LEVELS',
                caption='Background Levels in ADU (average image area fluences across '+\
                        'control frames pre- \& post- saturation).'
                    )
        

    def basic_analysis(self):
        """

        Basic analysis of data.

        **METACODE**

        ::

            f.e.CCD:
                f.e.Q:
                    use SATURATED frame to generate pixel saturation MASK
                    measure stats in pix satur MASK across OBSIDs 
                     (pre-satur, satur, post-satur)


        """

        raise NotImplementedError

    def meta_analysis(self):
        """

        Meta-analysis of data.

        **METACODE**

        ::


            f.e.CCD:
                f.e.Q:
                    estimate delta-charge_0 and decay tau from time-series

            report:  
                persistence level (delta-charge_0) and time constant


        """
        raise NotImplementedError
