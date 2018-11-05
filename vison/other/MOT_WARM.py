#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: MOT_WARM

Readiness verification: Warm Test before Cooling Down. 


Created on Mon Oct 22 17:11:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import unittest
from matplotlib.colors import Normalize


from vison.datamodel import core
from vison.pipe.task import Task
from vison.point.PointTask import PointTask
from vison.support import context
from vison.datamodel import scriptic as sc
from vison.datamodel import ccd
#import B01aux
from vison.dark.DarkTask import DarkTask
from vison.datamodel import inputs, cdp
from vison.support import utils
from vison.support.files import cPickleRead

# END IMPORT

isthere = os.path.exists


HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


MW_commvalues = dict(program='CALCAMP', test='MOT_WARM', 
                         flushes=7, siflsh=1, siflsh_p=500,
                         inisweep=1,
                         vstart=0, vend=2086,
                         toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                         chinj=0,
                         s_tpump=0,
                         v_tpump=0,
                         exptime=0., 
                         shuttr=0, 
                         e_shuttr=0,
                         mirr_on=0,
                         wave=4,
                         motr_on=0,
                         source='flat',
                         comments='BIAS')


class MOT_WARM_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    #manifesto.update(OrderedDict(sorted([
    #    ('N', ([int], 'Number of Frame Acquisitions.')),
    #])))


class MOT_WARM(DarkTask):
    """ """

    inputsclass = MOT_WARM_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [
                         ('lock', self.lock_on_stars),
                         ('check', self.check_data)]
#                         ('basic', self.basic_analysis)]

        super(MOT_WARM, self).__init__(inputs, log, drill, debug)
        self.name = 'MOT_WARM'
        self.type = 'Simple'
        
        
        self.HKKeys = HKKeys        
        #self.figdict = B01aux.B01figs.copy()
        #self.CDP_lib = B01aux.CDP_lib.copy()
        self.inputs['subpaths'] = dict(figs='figs',
                                       profiles='profiles', 
                                       products='products')
    
    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass()


    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds MOT_WARM script structure dictionary.
        
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """
        
        toi_ro = 250
        
        IDLinj=11.
        IDHinj=18.
        IG1inj=3.5
        IG2inj=7.0
        toi_ch=250.
        
        waveflat = 800
        FW_IDflat = self.ogse.get_FW_ID(waveflat)
        FW_IDflatx = int(FW_IDflat[-1])
        exptimeFL = self.ogse.profile['tFWC_flat']['nm%i' % waveflat]/2.
        
        
        
        MW_sdict = dict(col001=dict(frames=1, exptime=0, rdmode='rwd_bas_vs',
                                    swellw=context.sumwell['rwd_bas_vs'][0],
                                    swelldly=context.sumwell['rwd_bas_vs'][1],
                                    vstart=0,vend=2086,toi_ro=toi_ro,
                                    shuttr=0, mirr_on=0,motr_on=0,
                                    source='flat',
                                    comments='BIAS'),
                        col002=dict(frames=1, exptime=0, rdmode='fwd_bas',
                                    swellw=context.sumwell['fwd_bas'][0],
                                    swelldly=context.sumwell['fwd_bas'][1],
                                    vstart=0,vend=2086,toi_ro=toi_ro,
                                    shuttr=0, mirr_on=0,motr_on=0,
                                    source='flat',
                                    comments='RAMP'),
                        col003=dict(frames=1, 
                                    IDL=IDLinj,IDH=IDHinj,
                                    IG1_1_T=IG1inj, IG1_2_T=IG1inj, IG1_3_T=IG1inj,
                                    IG1_1_B=IG1inj, IG1_2_B=IG1inj, IG1_3_B=IG1inj,
                                    IG2_T=IG2inj,IG2_B=IG2inj,
                                    toi_ch=toi_ch,toi_ro=toi_ro,
                                    id_wid=60, id_dly=toi_ch*2.5,
                                    chin_dly=0,chinj=1,chinj_on=30,chinj_of=50,
                                    exptime=0,shuttr=0,vstart=0,vend=2086,
                                    source='flat',
                                    comments='CHINJ'),
                        col004=dict(frames=1, exptime=exptimeFL, 
                                    wave=FW_IDflatx,
                                    vstart=0,vend=2086,toi_ro=toi_ro,
                                    shuttr=1, mirr_on=0,motr_on=0,
                                    source='flat',
                                    comments='FLAT'))
                        
        for i, wavenm in enumerate([590,730,880]):
            colnr = i+5
            FWIDx = int(self.ogse.get_FW_ID(wavenm)[-1])
            iexptimeps = self.ogse.profile['tFWC_point']['nm%i' % wavenm] * 0.5
            imirror = self.ogse.profile['mirror_nom']['F%i' % FWIDx]
            MW_sdict['col%03i' % colnr] = dict(frames=1,
                    wave=FWIDx,exptime=iexptimeps,
                    vstart=0,vend=2086,toi_ro=toi_ro,
                    shuttr=1,mirr_on=1,mirr_pos=imirror,
                    source='point',
                    comments='PNT')
        
        Ncols = len(MW_sdict.keys())
        MW_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(MW_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        MW_sdict = sc.update_structdict(
            MW_sdict, commvalues, diffvalues)

        return MW_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(MOT_WARM, self).filterexposures(structure, explog, OBSID_lims, 
                                    colorblind=False,
                                    wavedkeys=wavedkeys)
    
    def lock_on_stars(self):
        """ """
        print('MOT_WARM.lock_on_stars: using HARDWIRE value of Observation Index. PLEASE REWORK!')
        iObs = 5 # HARDWIRED
        PointTask.lock_on_stars(self,iObs=iObs)
    
    
    def check_data(self):
        """ """
        kwargs = dict(figkeys=[])        
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

        newcolnames_std = ['std_pre', 'std_ove']
        for newcolname_std in newcolnames_std:
            self.dd.initColumn(newcolname_std, Xindices,
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

                        for reg in ['pre', 'ove']:
                            stats = ccdobj.get_stats(Quad, sector=reg, statkeys=['median', 'std'], trimscan=[5, 5],
                                                     ignore_pover=True, extension=-1, VSTART=vstart, VEND=vend)
                            self.dd.mx['offset_%s' %
                                       reg][iObs, jCCD, kQ] = stats[0]
                            self.dd.mx['std_%s' %
                                       reg][iObs, jCCD, kQ] = stats[1]

    def check_metrics_ST(self, **kwargs):
        """ 

        """

        Xindices = self.dd.indices
        CCDs = Xindices.get_vals('CCD')

        if self.report is not None:
            self.report.add_Section(
                keyword='check_ronoffset', Title='Offsets and RON', level=1)

        # absolute value of offsets

        offsets_lims = self.perflimits['offsets_lims']

        regs_off = ['pre', 'ove']

        for reg in regs_off:
            arr = self.dd.mx['offset_%s' % reg]
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
                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)

        # cross-check of offsets: referred to pre-scan

        offsets_gradients = self.perflimits['offsets_gradients']

        regs_grad = ['ove']

        for ireg, reg in enumerate(regs_grad):
            _lims = dict()
            for CCDk in CCDs:
                _lims[CCDk] = offsets_gradients[CCDk][reg]
            arr = self.dd.mx['offset_%s' % reg][:]-self.dd.mx['offset_pre'][:]
            _xcheck_offsets = self.check_stat_perCCDandQ(arr, _lims, CCDs)
            
            self.addComplianceMatrix2Self(_xcheck_offsets,'offsets_grad_%s' % reg)

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



        

    def basic_analysis(self):
        """ 
        EXPOSURES:
            BIAS, RAMP, CHINJ, FLAT, POINT_590, POINT_730, POINT_880
            
        """
        
        

