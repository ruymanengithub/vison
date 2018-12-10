#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: COSMETICS00

Cosmetics masks for the detectors: defects in darkness and PR


Created on Wed Dec 05 17:59:00 2018

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
from scipy import interpolate

from vison.other import COSaux
from vison.pipe.task import HKKeys
from vison.image import bits
from vison.datamodel import core
from vison.pipe.task import Task
from vison.support import context
from vison.datamodel import scriptic as sc
from vison.datamodel import ccd
from vison.dark.DarkTask import DarkTask
from vison.datamodel import inputs, cdp
from vison.support import utils
from vison.scripts import vis_cosmetics_masker as vcm
# END IMPORT

isthere = os.path.exists


#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


COS_commvalues = dict(program='CALCAMP', test='COSMETICS00', 
                         rdmode='fwd_bas',
                         flushes=7, siflsh=1, siflsh_p=500,
                         swellw= context.sumwell['fwd_bas'][0],
                         swelldly= context.sumwell['fwd_bas'][1],
                         inisweep=1,
                         vstart=0, vend=2086,
                         toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                         chinj=0,
                         s_tpump=0,
                         v_tpump=0,
                         e_shuttr=0,
                         mirr_on=0,
                         wave=4,
                         motr_on=0,
                         source='flat')



class COS_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    #manifesto.update(OrderedDict(sorted([
    #    ('N', ([int], 'Number of Frame Acquisitions.')),
    #])))


class COSMETICS00(DarkTask):
    """ """

    inputsclass = COS_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [
                         ('check', self.check_data),
                         ('masks', self.do_masks)]

        super(COSMETICS00, self).__init__(inputs, log, drill, debug)
        self.name = 'COSMETICS00'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys        
        self.figdict = COSaux.get_COS_figs()
        self.CDP_lib = COSaux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')
    
    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass()


    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds COSMETICS00 script structure dictionary.
        
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """
        
        waveflat = 800
        FW_IDflat = self.ogse.get_FW_ID(waveflat)
        FW_IDflatx = int(FW_IDflat[-1])
        exptimeFL = self.ogse.profile['tFWC_flat']['nm%i' % waveflat]*0.5
        
        
        COS_sdict = OrderedDict()
        
        COS_sdict['col001'] = dict(frames=5, exptime=0, 
                                    shuttr=0, mirr_on=0,motr_on=0,
                                    source='flat',
                                    comments='BIAS')

        
        COS_sdict['col002']=dict(frames=5, exptime=exptimeFL, 
                                    wave=FW_IDflatx,
                                    shuttr=1, mirr_on=0,motr_on=0,
                                    source='flat',
                                    comments='FLAT')
                                
        Ncols = len(COS_sdict.keys())
        COS_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(COS_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()
        
        
        COS_sdict = sc.update_structdict(
            COS_sdict, commvalues, diffvalues)
        

        return COS_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = []
        return super(COSMETICS00, self).filterexposures(structure, explog, OBSID_lims, 
                                    colorblind=True,
                                    wavedkeys=wavedkeys)


    def check_data(self):
        pass

#    def check_data(self):
#        """ """
#        kwargs = dict(figkeys=[])        
#        Task.check_data(self, **kwargs)
        

#    def get_checkstats_ST(self, **kwargs):
#        """ """
#        
#        # Initialize new columns
#
#        Xindices = copy.deepcopy(self.dd.indices)
#
#        if 'Quad' not in Xindices.names:
#            Xindices.append(core.vIndex('Quad', vals=context.Quads))
#
#        valini = 0.
#
#        newcolnames_off = ['offset_pre', 'offset_img', 'offset_ove']
#        for newcolname_off in newcolnames_off:
#            self.dd.initColumn(newcolname_off, Xindices,
#                               dtype='float32', valini=valini)
#
#        newcolnames_std = ['std_pre', 'std_img', 'std_ove']
#        for newcolname_std in newcolnames_std:
#            self.dd.initColumn(newcolname_std, Xindices,
#                               dtype='float32', valini=valini)
#
#        nObs, _, _ = Xindices.shape
#        CCDs = Xindices.get_vals('CCD')
#        Quads = Xindices.get_vals('Quad')
#
#        # Get statistics in different regions
#
#        if not self.drill:
#
#            for iObs in range(nObs):
#                for jCCD, CCDk in enumerate(CCDs):
#                    dpath = self.dd.mx['datapath'][iObs, jCCD]
#                    ffits = os.path.join(dpath, '%s.fits' %
#                                         self.dd.mx['File_name'][iObs, jCCD])
#                    ccdobj = ccd.CCD(ffits)
#                    vstart = self.dd.mx['vstart'][iObs][jCCD]
#                    vend = self.dd.mx['vend'][iObs][jCCD]
#
#                    for kQ, Quad in enumerate(Quads):
#
#                        for reg in ['pre', 'img', 'ove']:
#                            stats = ccdobj.get_stats(Quad, sector=reg, statkeys=['median', 'std'], trimscan=[5, 5],
#                                                     ignore_pover=True, extension=-1, VSTART=vstart, VEND=vend)
#                            self.dd.mx['offset_%s' %
#                                       reg][iObs, jCCD, kQ] = stats[0]
#                            self.dd.mx['std_%s' %
#                                       reg][iObs, jCCD, kQ] = stats[1]
#
#    def check_metrics_ST(self, **kwargs):
#        """ 
#
#        """
#
#        Xindices = self.dd.indices
#        CCDs = Xindices.get_vals('CCD')
#
#        if self.report is not None:
#            self.report.add_Section(
#                keyword='check_ronoffset', Title='Offsets and RON', level=1)
#
#        # absolute value of offsets
#
#        offsets_lims = self.perflimits['offsets_lims']
#        
#        BIAS_ix = ObsIDdict['BIAS']
#
#        regs_off = ['pre', 'ove']
#
#        for reg in regs_off:
#            arr = self.dd.mx['offset_%s' % reg][BIAS_ix+1:,...].copy()
#            _compliance_offsets = self.check_stat_perCCDandQ(
#                arr, offsets_lims, CCDs)
#            
#            self.addComplianceMatrix2Self(_compliance_offsets,'offsets_%s' % reg)
#
#            # if not self.IsComplianceMatrixOK(_compliance_offsets):
#            if not _compliance_offsets.IsCompliant():
#                self.dd.flags.add('POORQUALDATA')
#            if self.log is not None:
#                self.addComplianceMatrix2Log(
#                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)
#            if self.report is not None:
#                self.addComplianceMatrix2Report(
#                    _compliance_offsets, label='COMPLIANCE OFFSETS [%s]:' % reg)
#
#        # cross-check of offsets: referred to pre-scan
#
#        offsets_gradients = self.perflimits['offsets_gradients']
#
#        regs_grad = ['ove']
#
#        for ireg, reg in enumerate(regs_grad):
#            _lims = dict()
#            for CCDk in CCDs:
#                _lims[CCDk] = offsets_gradients[CCDk][reg]
#            arr = self.dd.mx['offset_%s' % reg][BIAS_ix+1:,...]-self.dd.mx['offset_pre'][BIAS_ix+1:,...]
#            _xcheck_offsets = self.check_stat_perCCDandQ(arr, _lims, CCDs)
#            
#            self.addComplianceMatrix2Self(_xcheck_offsets,'offsets_grad_%s' % reg)
#
#            if not self.IsComplianceMatrixOK(_xcheck_offsets):
#                self.dd.flags.add('POORQUALDATA')
#            if self.log is not None:
#                self.addComplianceMatrix2Log(
#                    _xcheck_offsets, label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)
#            if self.report is not None:
#                self.addComplianceMatrix2Report(
#                    _xcheck_offsets, label='OFFSET GRAD [%s-PRE] COMPLIANCE:' % reg)
#
#        # absolute value of std
#
#        regs_std = ['pre', 'ove']
#
#        RONs_lims = self.perflimits['RONs_lims']
#        for reg in regs_std:
#            _compliance_std = self.check_stat_perCCDandQ(
#                self.dd.mx['std_%s' % reg][BIAS_ix+1:,...], RONs_lims, CCDs)
#            
#            self.addComplianceMatrix2Self(_compliance_std,'std_%s' % reg)
#
#
#            if not self.IsComplianceMatrixOK(_compliance_std):
#                self.dd.flags.add('POORQUALDATA')
#                self.dd.flags.add('RON_OOL')
#            if self.log is not None:
#                self.addComplianceMatrix2Log(
#                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)
#            if self.report is not None:
#                self.addComplianceMatrix2Report(
#                    _compliance_std, label='COMPLIANCE RON [%s]:' % reg)

    def do_masks(self):
        """
            
        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='masks', Title='COSMETICS00: MASKS', level=0)
        
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD = DDindices.shape
        #Quads = DDindices.get_vals('Quad')
        CCDs = DDindices.get_vals('CCD')
        
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        
        productspath = self.inputs['subpaths']['products']
        
        dpath = self.dd.mx['datapath'][0,0]
        
        for jCCD, CCDk in enumerate(CCDs):
            
            
            seldk = np.where(self.dd.mx['exptime'][:,0] == 0.)
            OBSID_list_dk = self.dd.mx['ObsID'][seldk].tolist()
            
            selfl = np.where(self.dd.mx['exptime'][:,0] > 0.)
            OBSID_list_fl = self.dd.mx['ObsID'][selfl].tolist()
        
            inputs = OrderedDict()
            
            inputs['sn_ccd'] = self.inputs['diffvalues']['sn_%s' % CCDk.lower()] 
            inputs['CCD'] = jCCD+1
            inputs['elvis'] = self.elvis
            inputs['tag'] = self.BLOCKID
            inputs['doDarkMask'] = True
            inputs['doFlatMask'] = True
            inputs['doMerge'] = True
            inputs['outpath'] = productspath
            inputs['inputs_dkmask'] = OrderedDict()
            inputs['inputs_dkmask']['datapath'] = dpath
            inputs['inputs_dkmask']['OBSID_list'] = OBSID_list_dk
            inputs['inputs_dkmask']['thresholds'] = [-1000.,50.0]
            inputs['inputs_flmask'] = OrderedDict()
            inputs['inputs_flmask']['datapath'] = dpath
            inputs['inputs_flmask']['OBSID_list'] = OBSID_list_fl
            inputs['inputs_flmask']['thresholds'] = [0.5,1.5]
            
            vcm.run_maskmaker(inputs)
            
        
        



