#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (vertical)

Created on Tue Aug 29 17:37:00 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import pandas as pd
import string as st
from matplotlib.colors import Normalize

from vison.datamodel import cdp
from vison.support import utils
from vison.pipe.task import HKKeys
from vison.pipe import lib as pilib
from vison.support import context
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
from vison.pipe.task import Task
from PumpTask import PumpTask
from vison.image import performance
from vison.datamodel import inputs
from vison.support.files import cPickleRead
import TP01aux
import tptools
# END IMPORT

isthere = os.path.exists

#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

IDL=11.
IDH=18.
IG1 = 4.
IG2 = 6.0

TP01_commvalues = dict(program='CALCAMP', test='TP01',
                       IDL=IDL, IDH=IDH,
                       IG1_1_T=IG1, IG1_2_T=IG1, IG1_3_T=IG1,
                       IG1_1_B=IG1, IG1_2_B=IG1, IG1_3_B=IG1,
                       IG2_T=IG2, IG2_B=IG2,
                       flushes=7, siflsh=1, siflsh_p=500,
                       inisweep=1,
                       vstart=0, vend=2086,
                       toi_fl=143., toi_ro=1000., toi_chinj=500,
                       chinj=1, chinj_on=2066, chinj_of=0,
                       id_wid=60,
                       v_tpump=1, v_tp_cnt=5000,                       
                       s_tpump=0,
                       exptime=0., shuttr=0,e_shuttr=0, 
                       mirr_on=0,
                       wave=4,
                       motr_on=0,
                       source='flat',
                       comments='')


def _thin_down(cat,N):
    """ """
    ckeys = cat.keys()
    catlen = len(cat[ckeys[0]])
    if N>=catlen:
        return cat.copy()
    else:
        thcat = OrderedDict()
        ixsel = (np.random.choice(np.arange(catlen),N),)
        
        for key in ckeys:
            thcat[key] = cat[key][ixsel].copy()
        return thcat

class TP01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
        ('Nshuffles_V',
         ([int], 'Number of Shuffles, Vertical/Parallel Pumping.')),
        ('id_delays',
         ([list], 'Injection Drain Delays [2, one per CCDs section].')),
        ('toi_tpv', ([list], 'Vector of TOI TP-V values.')),
        ('vpumpmodes',
         ([list], 'Vertical/Parallel Pumping Starting points.'))
    ])))




class TP01(PumpTask):
    """ """

    inputsclass = TP01_inputs
    contrast_threshold = 0.01

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('injection', self.charact_injection),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('debugtask',self.debugtask),
                         ('meta', self.meta_analysis)]
        super(TP01, self).__init__(inputs=inputs, log=log, 
            drill=drill, debug=debug, cleanafter=cleanafter)
        self.name = 'TP01'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = TP01aux.get_TP01figs()
        self.CDP_lib = TP01aux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       products='products')
        

    def set_inpdefaults(self, **kwargs):
        """ """

        toi_chinjTP01 = 250
        self.inpdefaults = dict(toi_chinj=toi_chinjTP01,
                                Nshuffles_V=5000,
                                id_delays=np.array([2.5, 1.5]) * toi_chinjTP01,
                                toi_tpv=[200, 1000, 2000, 4000, 8000],
                                vpumpmodes=[123, 234, 341, 412])

    def set_perfdefaults(self, **kwargs):
        super(TP01, self).set_perfdefaults(**kwargs)

        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()

        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ """

        Nshuffles_V = self.inputs['Nshuffles_V']
        toi_tpv = self.inputs['toi_tpv']
        id_delays = self.inputs['id_delays']
        vpumpmodes = self.inputs['vpumpmodes']
        toi_chinj = self.inputs['toi_chinj']

        assert len(id_delays) == 2

        TP01_sdict = dict()

        TP01_commvalues['v_tp_cnt'] = Nshuffles_V

        # First Injection Drain Delay

        TP01_sdict['col001'] = dict(frames=1, v_tpump=0, comments='BGD',
                                  id_dly=id_delays[0], toi_ch=toi_chinj)

        colcounter = 2
        for i, toi_tp in enumerate(toi_tpv):
            for k, vpumpmode in enumerate(vpumpmodes):
                colkey = 'col%03i' % colcounter
                TP01_sdict[colkey] = dict(frames=1, toi_tp=toi_tp,
                                          id_dly=id_delays[0], v_tpmod=vpumpmode,
                                          toi_ch=toi_chinj)

                colcounter += 1

        # Second Injection Drain Delay

        TP01_sdict['col%03i' % colcounter] = dict(frames=1, v_tpump=0, comments='BGD',
                                                id_dly=id_delays[1], toi_ch=toi_chinj)
        colcounter += 1

        for j, toi_tp in enumerate(toi_tpv):

            for k, vpumpmode in enumerate(vpumpmodes):

                colkey = 'col%03i' % colcounter
                #print colkey
                TP01_sdict[colkey] = dict(frames=1, toi_tp=toi_tp,
                                          id_dly=id_delays[1], v_tpmod=vpumpmode,
                                          toi_ch=toi_chinj)

                colcounter += 1

        Ncols = len(TP01_sdict.keys())
        TP01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        TP01_sdict = sc.update_structdict(TP01_sdict, commvalues, diffvalues)

        return TP01_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(TP01, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                 wavedkeys=wavedkeys)

    
    def prepare_images(self):
        super(TP01, self).prepare_images(doExtract=True, 
             doBadPixels=True,
             doMask=True, # False ON TESTS!
             doOffset=True, 
             doBias=False, 
             doFF=False)
    
    def charact_injection(self):
        """Characterises Charge Injection."""
        
        
        if self.report is not None:
            self.report.add_Section(
                keyword='charact', Title='TP01: Charge Injection Characterisation', level=0)
        
        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        Quads = DDindices.get_vals('Quad')

        ccdpicklespath = self.inputs['subpaths']['ccdpickles']
        productspath = self.inputs['subpaths']['products']
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        id_dlys = np.unique(self.dd.mx['id_dly'][:, 0])
        
        char_res_dict = OrderedDict()
        chinjnoise_dd = OrderedDict()
        
        if not self.drill:
            
            for id_dly in id_dlys:
            
                for jCCD, CCDk in enumerate(CCDs):
                    
                    chinjnoise_dd[CCDk] = OrderedDict()
                    
                    ixsel = np.where((self.dd.mx['id_dly'][:] == id_dly) & 
                        (self.dd.mx['v_tpump'][:] == 0) & 
                        (self.dd.mx['CCD'][:] == CCDk))
                    
                    iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ixsel][0]
                    iccdobj = cPickleRead(os.path.join(ccdpicklespath,iccdobj_f))
                    
                    char_res_dict[CCDk] = tptools.charact_injection(iccdobj)
                    
                    for Q in Quads:                    
                        chinjnoise_dd[CCDk][Q] = np.nanmedian(char_res_dict[CCDk][Q]['injnoise'])                    
        
        # Saving the charge injection characterisation results
        
        chchar_cdp = self.CDP_lib['CHINJCHARACT']
        
        chchar_cdp.header = CDP_header.copy()
        chchar_cdp.meta = dict()
        chchar_cdp.path = productspath
        chchar_cdp.data = char_res_dict.copy()
        
        self.save_CDP(chchar_cdp)
        self.pack_CDP_to_dd(chchar_cdp,'CHINJCHARACT')
        
        # REPORTS
        
        # Table with injection noise per CCD/Quadrant
                
        chinjnoise_cdp = self.CDP_lib['CHINJNOISE']
        chinjnoise_cdp.path = self.inputs['subpaths']['products']
        chinjnoise_dddf = OrderedDict(CHINJNOISE=pd.DataFrame.from_dict(chinjnoise_dd))
        chinjnoise_cdp.ingest_inputs(data=chinjnoise_dddf.copy(),                             
                              meta=dict(),
                              header=CDP_header.copy())
        
        chinjnoise_cdp.init_wb_and_fillAll(header_title='%s: CHARGE INJECTION NOISE' % self.inputs['test'])
        self.save_CDP(chinjnoise_cdp)
        self.pack_CDP_to_dd(chinjnoise_cdp, 'CHINJNOISE_CDP')

        if self.report is not None:
            
            ff = lambda x: '%.2f' % x
            
            formatters = [ff,ff,ff]
            
            CHINJNOISEtex = chinjnoise_cdp.get_textable(sheet='CHINJNOISE', 
                                            caption='%s: Charge Injection Noise (ADU, rms).' % \
                                                             self.inputs['test'],
                                            longtable=False, 
                                            fitwidth=True,
                                            index=True,
                                            formatters=formatters)
            self.report.add_Text(CHINJNOISEtex)
            

        

    def extract(self):
        """ 

        Obtain maps of dipoles.

        **METACODE**

        ::

            f.e. id_delay (there are 2):
                f.e. CCD:
                    f.e. Q:
                        produce reference non-pumped injection map

            f. e. ObsID:
                f.e. CCD:

                    load ccdobj                    
                    f.e.Q.:
                        divide ccdobj.Q by injection map

                    save dipole map and store reference


        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='TP01 Extraction', level=0)
        
        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')

        ccdpicklespath = self.inputs['subpaths']['ccdpickles']
        productspath = self.inputs['subpaths']['products']

        # Initialisations

        self.dd.initColumn('dipoles_raw', self.dd.mx['ccdobj_name'].indices,
                           dtype='S100', valini='None')

        id_dlys = np.unique(self.dd.mx['id_dly'][:, 0])

        if not self.drill:

            # Computing maps of relative amplitude of dipoles

            for id_dly in id_dlys:
                
                for jCCD, CCDk in enumerate(CCDs):

                    ixsel = np.where((self.dd.mx['id_dly'][:] == id_dly) & (
                        self.dd.mx['v_tpump'][:] != 0))

                    for ix in ixsel[0]:
                        ObsID = self.dd.mx['ObsID'][ix]
                        vstart = self.dd.mx['vstart'][ix, jCCD]
                        vend = self.dd.mx['vend'][ix,jCCD]


                        ioutf = 'TP01_rawmap_%i_IDDLY_%i_ROE1_%s' % (
                            ObsID, id_dly, CCDk)

                        iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ix, jCCD]

                        try:
                            iccdobj = cPickleRead(
                                os.path.join(ccdpicklespath, iccdobj_f))
                            irawmap = copy.deepcopy(iccdobj)

                            irawmap = tptools.gen_raw_dpmap_vtpump(
                                irawmap, Navgrows=-1, vstart=vstart, vend=vend)

                            irawmap.writeto(os.path.join(productspath, \
                                                         '%s.fits' % ioutf))

                            self.dd.mx['dipoles_raw'][ix, jCCD] = ioutf

                        except:  # TESTS
                            pass

        if self.report is not None:
            self.report.add_Text('All Done!')

    def basic_analysis(self):
        """

        Basic analysis of data.

        **METACODE**

        ::

            f. e. ObsID [there are different TOI_TP and TP-patterns]:
                f.e.CCD:
                    f.e.Q:
                        load "map of relative pumping"
                        find_dipoles:
                            x, y, rel-amplitude, orientation

            produce & report:  
                map location of dipoles
                PDF of dipole amplitudes (for N and S)
                Counts of dipoles (and N vs. S)

        """
        
        
        if self.report is not None:
            self.report.add_Section(
                keyword='basic', Title='TP01 Basic Analysis', level=0)
        
        threshold = self.contrast_threshold
        #CCDhalves = ['top','bottom']
        _Quads_dict = dict(bottom = ['E','F'],
                           top = ['G','H'])
        
        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        allQuads = DDindices.get_vals('Quad')

        productspath = self.inputs['subpaths']['products']


        id_dlys = np.unique(self.dd.mx['id_dly'][:, 0])
        
        mods = np.unique(self.dd.mx['v_tp_mod'][:,0])
        modkeys = ['m%i' % item for item in mods]
        
        tois = np.unique(self.dd.mx['toi_tp'][:,0])
        toikeys = ['u%04i' % item for item in tois]
        
        # initialisation
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        masterdict = OrderedDict()
        
        for CCDk in CCDs:
            masterdict[CCDk] = OrderedDict()
            for Q in allQuads:
                masterdict[CCDk][Q] = OrderedDict()
                for modkey in modkeys:
                    masterdict[CCDk][Q][modkey] = OrderedDict()
                    for toikey in toikeys:
                        masterdict[CCDk][Q][modkey][toikey] = None
        
        
        if not self.drill:

            # Getting dipole catalogues
            
            ontests = False
            print('WARNING: TP01.basic_analysis incomplete, TESTS')
            
            if not ontests:

                for id_dly in id_dlys:
                    
                    print('idl_dly = %s' % id_dly)
    
                    for jCCD, CCDk in enumerate(CCDs):
                        
                        print('CCD=%s' % CCDk)
    
                        ixsel = np.where((self.dd.mx['id_dly'][:,0] == id_dly) & (
                            self.dd.mx['v_tpump'][:,0] != 0))
                        
                        
                        for ix in ixsel[0]:
                            
                            ObsID = self.dd.mx['ObsID'][ix]                            
                            vstart = self.dd.mx['vstart'][ix, jCCD]
                            vend = self.dd.mx['vend'][ix,jCCD]
                            toi_ch = float(self.dd.mx['toi_ch'][ix,jCCD])
                            v_tp_mod = self.dd.mx['v_tp_mod'][ix,jCCD]
                            toi_tp = self.dd.mx['toi_tp'][ix,jCCD]
                            
                            modkey = 'm%i' % v_tp_mod
                            toikey = 'u%04i' % toi_tp
                            
                            if np.isclose(id_dly/toi_ch,1.5):
                                CCDhalf = 'top'
                            elif np.isclose(id_dly/toi_ch,2.5):
                                CCDhalf = 'bottom'
                                
                            Quads = _Quads_dict[CCDhalf]
    
    
                            imapf = 'TP01_rawmap_%i_IDDLY_%i_ROE1_%s.fits' % (
                                ObsID, id_dly, CCDk)
                            
                            imapccdobj = ccd.CCD(os.path.join(productspath,imapf))
                            
                            print('OBSID=%s, %s' % (ObsID,imapf))
    
                            
                            for iQ, Q in enumerate(Quads):
                               
                                    
                                idd = tptools.find_dipoles_vtpump(imapccdobj,threshold,
                                      Q,vstart=vstart,vend=vend,
                                      extension=-1)
                                 
                                masterdict[CCDk][Q][modkey][toikey] = idd.copy()
                                
                                    
            
            df = tptools._aggregate(masterdict,CCDs,allQuads,modkeys,toikeys,'toi')
            
            
            summaryMean = df.groupby(level=['CCD','Q','mod']).mean()
            summaryTot = df.groupby(level=['CCD','Q','mod']).sum()
        
            summary = summaryTot.copy()
            summary['R'] = summaryMean['R'].copy()
            summary['A'] = summaryMean['A'].copy()
            summary.columns = ['N','<R>','<A>']
        
            summtxtreport = tptools._get_txt(summary)
            
        
            if self.report is not None:
                self.report.add_Text('Aggregated Dipole Statistics: Number, \<Ratio N/S\>, \<Amplitude\>')
                self.report.add_Text(summtxtreport)
            
            
            # Saving the MasterCat
            
            MCmeta = OrderedDict()
            MCmeta['THRESHOLD'] = threshold
            MCmeta['Quads'] = allQuads
            MCmeta['modkeys'] = modkeys
            MCmeta['toikeys'] = toikeys
            
            
            for CCDk in CCDs:
            
                kmastercat = self.CDP_lib['MASTERCAT_%s' % CCDk]
                kmastercat.header = CDP_header.copy()
                kmastercat.meta = MCmeta
                kmastercat.path = productspath
                kmastercat.data = masterdict[CCDk].copy()
                
                self.save_CDP(kmastercat)
                self.pack_CDP_to_dd(kmastercat,'MASTERCAT_%s' % CCDk)
                
                
    
    def debugtask(self):
        
        if self.report is not None:
            self.report.add_Section(
                keyword='debug', Title='TP01 DEBUG', level=0)
        
        CCDs = ['CCD1','CCD2','CCD3']
        
        masterdict = OrderedDict()
        for CCDk in CCDs:
            kmastercatpick = self.dd.products['MASTERCAT_%s' % CCDk]
            masterdict[CCDk] = cPickleRead(kmastercatpick)['data'].copy()
                
        allQuads = ['E','F','G','H']
        modkeys = ['m123', 'm234', 'm341', 'm412']
        toikeys = ['u0200', 'u1000', 'u2000', 'u4000', 'u8000']
              
        df =_aggregate(masterdict,CCDs,allQuads,modkeys,toikeys)
        
        
        
        summaryMean = df.groupby(level=['CCD','Q','mod']).mean()
        summaryTot = df.groupby(level=['CCD','Q','mod']).sum()
        
        summary = summaryTot.copy()
        summary['R'] = summaryMean['R'].copy()
        summary['A'] = summaryMean['A'].copy()
        summary.columns = ['N','<R>','<A>']
        
        txtreport = _get_txt(summary)
            
        
        if self.report is not None:
            self.report.add_Text('Aggregated Dipole Statistics: Number, Ratio N/S, \<Amplitude\>')
            self.report.add_Text(txtreport)
        

    def meta_analysis(self):
        """

        Meta-analysis of data:

            Try to identify tau and pixel-phase location for each trap.
            Need to associate dipoles across TOI_TPs and TP-patterns


        **METACODE**

        ::

            across TOI_TP, patterns:
                
                build catalog of traps: x,y, tp-mode, tau, Pc
                tau, Pc = f({A,TOI})

            Report on :
                Histogram of Taus
                Histogram of Pc (capture probability)
                Histogram of I-phases (larger phases should have more traps, 
                                  statistically) -> check

                Total Count of Traps
                
                

        """
        
        
        
        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='TP01 Meta Analysis', level=0)
        
        
        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        allQuads = DDindices.get_vals('Quad')

        productspath = self.inputs['subpaths']['products']
        Nshuffles = self.inputs['Nshuffles_V']
        
        
        # initialisation
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        mods = np.unique(self.dd.mx['v_tp_mod'][:,0])
        modkeys = ['m%i' % item for item in mods]
        
        tois = np.unique(self.dd.mx['toi_tp'][:,0])
        toikeys = ['u%04i' % item for item in tois]
        Ampcols = ['A_%s' % toikey for toikey in toikeys]
        
        
        onTests = False
        
        if not onTests:
        
            if not self.drill:
                
       
                mergecat = OrderedDict()
                
                
                print('\nMerging Dipole Catalogs...\n')
                
                for jCCD, CCDk in enumerate(CCDs):
                    
                    
                    mastercatpick = self.dd.products['MASTERCAT_%s' % CCDk]
        
                    masterdata = cPickleRead(mastercatpick)['data'].copy()
                    
                    mergecat[CCDk] = OrderedDict()
                    
                    for iQ, Q in enumerate(allQuads):
                        
                        mergecat[CCDk][Q] = OrderedDict()
                        
                        rawcatCQ = masterdata[Q].copy()
                        
                        # Pandizing the toi catalogs
                        
                        for modkey in modkeys:
                            
                            for toikey in toikeys:
                                
                                if onTests:
                                    rawcatCQ[modkey][toikey] = _thin_down(rawcatCQ[modkey][toikey],1000)
                                    
                                rawcatCQ[modkey][toikey] =\
                                        pd.DataFrame.from_dict(rawcatCQ[modkey][toikey])
                        
                        # Merging the toi catalogs
                        
                        for modkey in modkeys:
                            
                            print('%s%s, %s...' % (CCDk,Q,modkey))
                            
#                            if onTests:
                                
#                                mergecat[CCDk][Q][modkey] = cPickleRead('vtpcat_CCD1_E_m123.pick')
#                                N = len(mergecat[CCDk][Q][modkey])
#                                ixsel = np.random.choice(np.arange(N),100)
#                                
#                                
#                                mergecat[CCDk][Q][modkey] = mergecat[CCDk][Q][modkey].iloc[ixsel]
#                                amplitudes = mergecat[CCDk][Q][modkey][Ampcols]
#                                
#                                Pc, tau = tptools.batch_fit_PcTau_vtp(amplitudes,tois,Nshuffles)
#                                
#                                mergecat[CCDk][Q][modkey]['Pc'] = pd.Series(Pc,
#                                        index=mergecat[CCDk][Q][modkey].index)
#                                
#                                mergecat[CCDk][Q][modkey]['tau'] = pd.Series(tau,
#                                        index=mergecat[CCDk][Q][modkey].index)
                                
                                
                            if not onTests:
                                
                                kqkmerged = tptools.merge_vtp_dipole_cats_bypos(
                                        rawcatCQ[modkey].copy(),
                                        toikeys[1:],toikeys[0])
                                
                                cols2drop = []
                                for toikey in toikeys:
                                    cols2drop += ['X_%s' % toikey,'Y_%s' % toikey,'S_%s' % toikey]
                                kqkmerged.drop(cols2drop,axis=1)
                                
    
                                amplitudes = kqkmerged[Ampcols]
                                
    
                                Pc, tau = tptools.batch_fit_PcTau_vtp(amplitudes,tois,Nshuffles)
                                
                                kqkmerged['Pc'] = pd.Series(Pc,
                                        index=kqkmerged.index)
                                
                                kqkmerged['tau'] = pd.Series(tau,
                                        index=kqkmerged.index)
                                
                                
                                mergecat[CCDk][Q][modkey] = kqkmerged
                                        
                                
            
        # Store output catalog(s) as a CDP
        
        if not onTests:
          
            for CCDk in CCDs:
                
                kmergedata = mergecat[CCDk].copy()
                colnames = kmergedata[allQuads[0]][modkeys[0]].keys().tolist()
                
                for Q in allQuads:
                    for modkey in modkeys:
                        kmergedata[Q][modkey] = kmergedata[Q][modkey].as_matrix()
                
                kmergecat = self.CDP_lib['MERGEDCAT_%s' % CCDk]
                kmergecat.header = CDP_header.copy()
                kmergecat.meta = OrderedDict(
                        Quadrants=allQuads,
                        modes=modkeys,
                        colnames=colnames)
                kmergecat.path = productspath
                kmergecat.data = kmergedata.copy()
                
                self.save_CDP(kmergecat)
                self.pack_CDP_to_dd(kmergecat,'MERGEDCAT_%s' % CCDk)
            
#        else:
#            
#            mergecat = OrderedDict()
#            
#            for CCDk in CCDs:
#                
#                kpick = os.path.join(productspath,'TP01_MergedCat_%s.pick' % CCDk)
#                kmergedcat = cPickleRead(kpick)
#                colnames = kmergedcat['meta']['colnames']
#                
#                mergecat[CCDk] = OrderedDict()
#                
#                for Q in allQuads:
#                    mergecat[CCDk][Q] = OrderedDict()
#                    
#                    for modkey in modkeys:
#                        mergecat[CCDk][Q][modkey] = pd.DataFrame(data=kmergedcat['data'][Q][modkey],index=np.arange(100),columns=colnames)
        
        # Produce Pc, tau heatmaps for each tp mode across CCD beam
            
        for modkey in modkeys:
            
            pltfig = self.figdict['TP01meta_%s' % modkey]
            
            pldata = OrderedDict()
            
            HeatmapPeaks = []
            
            ixtau = colnames.index('tau')
            ixPc = colnames.index('Pc')
            
            for CCDk in CCDs:
                pldata[CCDk] = OrderedDict()
                for Q in allQuads:
                    pldata[CCDk][Q] = OrderedDict()
                    
                    logtau = np.log10(mergecat[CCDk][Q][modkey][:,ixtau].copy())                    
                    logPc = np.log10(mergecat[CCDk][Q][modkey][:,ixPc].copy())
                    
                    Heatmap, xedges, yedges = np.histogram2d(logtau, logPc, bins=(100,50), 
                range=[[1.5,5.5],[-6,-3]]) 
                    
                    if CCDk == CCDs[0] and Q == allQuads[0]:
                        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    
                    HeatmapPeaks.append(np.nanmax(Heatmap))
                    
                    
                    pldata[CCDk][Q]['img'] = Heatmap.copy()
            
            if onTests:
                for CCDk in ['CCD2','CCD3']:
                    pldata[CCDk] = pldata['CCD1'].copy()
            
            pltfig[1]['data'] = pldata.copy()
            
            normfunction = Normalize(vmin=1,vmax=max(HeatmapPeaks)/2.)
                    
            pltfig[1]['meta']['corekwargs']['norm'] = normfunction
            pltfig[1]['meta']['corekwargs']['extent'] = extent
                  
        if self.report is not None:
            Mfigkeys = ['TP01meta_%s' % mkey for mkey in modkeys]
            self.addFigures_ST(figkeys=Mfigkeys,
                           dobuilddata=False)
        
        self.canbecleaned = True