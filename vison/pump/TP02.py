#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP02

Trap-Pumping calibration (serial)

Created on Tue Aug 29 17:38:00 2017

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

from vison.pipe.task import HKKeys
from vison.support import utils
from vison.pipe import lib as pilib
from vison.support import context
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
#from vison.pipe.task import Task
from PumpTask import PumpTask
from vison.image import performance
from vison.datamodel import inputs
from vison.support.files import cPickleRead
import TP02aux
import tptools
# END IMPORT

isthere = os.path.exists

#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

IDL = 11.
IDH = 18.
IG1 = 4.
IG2 = 6.0

TP02_commvalues = dict(program='CALCAMP', test='TP02',
                       IDL=IDL, IDH=IDH,
                       IG1_1_T=IG1, IG1_2_T=IG1, IG1_3_T=IG1,
                       IG1_1_B=IG1, IG1_2_B=IG1, IG1_3_B=IG1,
                       IG2_T=IG2, IG2_B=IG2,
                       flushes=7, siflsh=1, siflsh_p=500,
                       inisweep=1,
                       vstart=0, vend=100,
                       toi_fl=143., toi_ro=1000., toi_chinj=500,
                       chinj=1, chinj_on=2066, chinj_of=0,
                       id_wid=60,
                       v_tpump = 0,
                       s_tpump=1, s_tp_cnt=5000,
                       exptime=0., shuttr=0, e_shuttr=0,
                       mirr_on=0,
                       wave=4,
                       motr_on=0,
                       source='flat',
                       comments='')




class TP02_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
        ('Nshuffles_H',
         ([int], 'Number of Shuffles, Horizontal/Serial Pumping.')),
        ('dwell_sv', ([list], 'Dwell Times list [serial].')),
        ('id_delays',
         ([list], 'Injection Drain Delays [2, one per CCDs section].')),
        ('spumpmodes',
         ([list], 'Horizontal/Serial Pumping Starting points.'))
    ])))


class TP02(PumpTask):
    """ """

    inputsclass = TP02_inputs
    contrast_threshold = 0.01

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(TP02, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug, 
                    cleanafter=cleanafter)
        self.name = 'TP02'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = TP02aux.get_TP02figs()
        self.CDP_lib = TP02aux.get_CDP_lib()
        self.inputs['subpaths'] = dict(figs='figs', 
                   ccdpickles='ccdpickles',
                   products='products')

    def set_inpdefaults(self, **kwargs):
        """ """
        toi_chinj = 250
        self.inpdefaults = dict(toi_chinj=toi_chinj,
                                Nshuffles_H=5000,
                                dwell_sv=[0., 4.75, 14.3, 28.6],
                                id_delays=np.array([2.5, 1.5])*toi_chinj,
                                spumpmodes=[23, 31])

    def set_perfdefaults(self, **kwargs):
        super(TP02, self).set_perfdefaults(**kwargs)

        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()

        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ """

        Nshuffles_H = self.inputs['Nshuffles_H']
        dwell_sv = self.inputs['dwell_sv']
        id_delays = self.inputs['id_delays']
        toi_chinj = self.inputs['toi_chinj']
        spumpmodes = self.inputs['spumpmodes']

        assert len(id_delays) == 2

        TP02_sdict = dict()

        TP02_commvalues['s_tp_cnt'] = Nshuffles_H

        # First Injection Drain Delay

        TP02_sdict['col001'] = dict(frames=1, v_tpump=0, s_tpump=0,
                                  comments='BGD', id_dly=id_delays[0], toi_ch=toi_chinj)

        colcounter = 2
        for i, dwell_s in enumerate(dwell_sv):

            for k, sermode in enumerate(spumpmodes):
                colkey = 'col%03i' % colcounter
                TP02_sdict[colkey] = dict(frames=1, dwell_s=dwell_s,
                                          v_tpump=0,s_tpump=1,
                                          id_dly=id_delays[0], s_tpmod=sermode, 
                                          toi_ch=toi_chinj)

                colcounter += 1

        # Second Injection Drain Delay

        TP02_sdict['col%03i' % colcounter] = dict(frames=1, v_tpump=0, s_tpump=0,
                                                comments='BGD', 
                                                id_dly=id_delays[1], toi_ch=toi_chinj)
        colcounter += 1

        for j, dwell_s in enumerate(dwell_sv):

            for k, sermode in enumerate(spumpmodes):

                colkey = 'col%03i' % colcounter
                #print colkey
                TP02_sdict[colkey] = dict(frames=1, dwell_s=dwell_s,
                                          v_tpump=0, s_tpump=1,
                                          id_dly=id_delays[1], s_tpmod=sermode, 
                                          toi_ch=toi_chinj)

                colcounter += 1

        Ncols = len(TP02_sdict.keys())
        TP02_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP02_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        TP02_sdict = sc.update_structdict(TP02_sdict, commvalues, diffvalues)

        return TP02_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(TP02, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                 wavedkeys=wavedkeys)
    
    def prepare_images(self):
        super(TP02, self).prepare_images(doExtract=True, 
             doBadPixels=True,
             doMask=False, # ON TESTS!
             doOffset=True, doBias=False, doFF=False)
    
    def extract(self):
        """ 
        Obtain Maps of Serial Dipoles.
        
        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='TP02 Extraction', level=0)
        
        
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
                    
                    ixsel = np.where((self.dd.mx['id_dly'][:,jCCD] == id_dly) & (
                        self.dd.mx['s_tpump'][:,jCCD] != 0))
                    ixref = np.where((self.dd.mx['id_dly'][:,jCCD] == id_dly) & (
                        self.dd.mx['s_tpump'][:,jCCD] == 0))
                    
                    Rvstart = self.dd.mx['vstart'][ixref[0][0], jCCD]
                    Rvend = self.dd.mx['vend'][ixref[0][0], jCCD]
                    
                    ccdref_f = '%s.pick' % self.dd.mx['ccdobj_name'][ixref[0][0], jCCD]
                    ccdref = cPickleRead( os.path.join(ccdpicklespath, ccdref_f))
                    
                    InjProfiles = OrderedDict()
                    
                    for Q in ccdref.Quads:
                        
                        InjProfiles[Q] = tptools.get_InjProfile(ccdref, 
                                   Q, Navgrows=-1, vstart=Rvstart, vend=Rvend, 
                                   extension=-1)

                    for ix in ixsel[0]:
                        
                        ObsID = self.dd.mx['ObsID'][ix]
                        vstart = self.dd.mx['vstart'][ix, jCCD]
                        vend = self.dd.mx['vend'][ix,jCCD]
                        s_tp_mod = self.dd.mx['s_tp_mod'][ix,jCCD]
                        
                        
                        ioutf = 'TP02_rawmap_%i_IDDLY_%i_ROE1_%s' % (
                            ObsID, id_dly, CCDk)

                        iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ix, jCCD]

                        iccdobj = cPickleRead(
                            os.path.join(ccdpicklespath, iccdobj_f))
                        irawmap = copy.deepcopy(iccdobj)

                        sh_injprofiles = OrderedDict()
                        for Q in ccdref.Quads:
                        
                            #if s_tp_mod == 31:
                            #    sh_injprofiles[Q] = np.roll(InjProfiles[Q],-1,0)
                            #elif s_tp_mod == 23:
                            #    sh_injprofiles[Q] = InjProfiles[Q].copy()
                        
                            sh_injprofiles[Q] = np.roll(InjProfiles[Q],-1,0)
                        
                        irawmap = tptools.gen_raw_dpmap_stpump(
                                irawmap, sh_injprofiles, vstart=vstart, vend=vend)

                        irawmap.writeto(os.path.join(productspath, \
                                '%s.fits' % ioutf),clobber=True)

                        self.dd.mx['dipoles_raw'][ix, jCCD] = ioutf
        
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
                        load raw 1D map of relative pumping (from extract_data)
                        identify dipoles:
                            x, rel-amplitude, orientation (E or W)

            produce & report:  
                map location of dipoles
                PDF of dipole amplitudes (for E and W)
                Counts of dipoles (and E vs. W)

        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='basic', Title='TP02 Basic Analysis', level=0)
        
        threshold = self.contrast_threshold
        #CCDhalves = ['top','bottom']
        _Quads_dict = dict(bottom = ['E','F'],
                           top = ['G','H'])
                
        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        allQuads = DDindices.get_vals('Quad')

        productspath = self.inputs['subpaths']['products']
        
        id_dlys = np.unique(self.dd.mx['id_dly'][:, 0])
        
        mods = np.unique(self.dd.mx['s_tp_mod'][:,0])
        modkeys = ['m%i' % item for item in mods]
        
        
        dwells = np.unique(self.dd.mx['dwell_s'][:,0])
        dwellkeys = ['u%04i' % item for item in dwells]
        
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
                    for dwellkey in dwellkeys:
                        masterdict[CCDk][Q][modkey][dwellkey] = None
        
        
        if not self.drill:

            # Getting dipole catalogues
            
            ontests = False
            print 'WARNING: TP02.basic_analysis incomplete, TESTS'
            
            if not ontests:
                
                for id_dly in id_dlys:
                    
                    print('idl_dly = %s' % id_dly)
    
                    for jCCD, CCDk in enumerate(CCDs):
                        
                        print('CCD=%s' % CCDk)
    
                        ixsel = np.where((self.dd.mx['id_dly'][:,0] == id_dly) & (
                            self.dd.mx['s_tpump'][:,0] != 0))
    
                        for ix in ixsel[0]:
                            ObsID = self.dd.mx['ObsID'][ix]
                            vstart = self.dd.mx['vstart'][ix, jCCD]
                            vend = self.dd.mx['vend'][ix,jCCD]
                            toi_ch = float(self.dd.mx['toi_ch'][ix,jCCD])
                            s_tp_mod = self.dd.mx['s_tp_mod'][ix,jCCD]
                            dwell = self.dd.mx['dwell_s'][ix,jCCD]
                            
                            modkey = 'm%i' % s_tp_mod
                            dwellkey = 'u%04i' % dwell
                            
                            if np.isclose(id_dly/toi_ch,1.5):
                                CCDhalf = 'top'
                            elif np.isclose(id_dly/toi_ch,2.5):
                                CCDhalf = 'bottom'
                                
                            Quads = _Quads_dict[CCDhalf]
                            
                            imapf = 'TP02_rawmap_%i_IDDLY_%i_ROE1_%s.fits' % (
                                ObsID, id_dly, CCDk)
                            
                            imapccdobj = ccd.CCD(os.path.join(productspath,imapf))    
                            
                            print('OBSID=%s, %s' % (ObsID,imapf))
                            
                            for iQ, Q in enumerate(Quads):
                                   
                                idd = tptools.find_dipoles_stpump(imapccdobj,threshold,
                                      Q,vstart=vstart,vend=vend,extension=-1)
                                 
                                masterdict[CCDk][Q][modkey][dwellkey] = idd.copy()
            

            df = tptools._aggregate(masterdict,CCDs,allQuads,modkeys,dwellkeys,'dwell')
            
            
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
            MCmeta['dwellkeys'] = dwellkeys
            
            
            for CCDk in CCDs:
            
                kmastercat = self.CDP_lib['MASTERCAT_%s' % CCDk]
                kmastercat.header = CDP_header.copy()
                kmastercat.meta = MCmeta
                kmastercat.path = productspath
                kmastercat.data = masterdict[CCDk].copy()
                
                self.save_CDP(kmastercat)
                self.pack_CDP_to_dd(kmastercat,'MASTERCAT_%s' % CCDk)
                

    def meta_analysis(self):
        """

        Meta-analysis of data:

            Try to identify tau and pixel-phase location for each trap.
            Need to associate dipoles across TOI_TPs and TP-patterns        


        **METACODE**

        ::

            across TOI_TP, patterns:
                build catalog of traps: x,y,R-phase, amp(dwell)
                from Amp(dwell) -> tau, Pc

            Report on :
               Histogram of Taus
               Histogram of Pc (capture probability)
               Histogram of R-phases

               Total Count of Traps



        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='TP02 Meta Analysis', level=0)
            
        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        allQuads = DDindices.get_vals('Quad')

        productspath = self.inputs['subpaths']['products']
        Nshuffles = self.inputs['Nshuffles_H']
        
        polarities = OrderedDict()
        polarities['North']=1.
        polarities['South']=0.
        
        # initialisation
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        mods = np.unique(self.dd.mx['s_tp_mod'][:,0])
        modkeys = ['m%i' % item for item in mods]
        
        dwells = np.unique(self.dd.mx['dwell_s'][:,0])
        dwellkeys = ['u%04i' % item for item in dwells]
        
        Ampcols = ['A_%s' % dwellkey for dwellkey in dwellkeys]
        
        
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
                        
                        # "Pandizing" the toi catalogs
                        
                        for modkey in modkeys:
                            
                            for dwellkey in dwellkeys:
                                                                    
                                rawcatCQ[modkey][dwellkey] =\
                                        pd.DataFrame.from_dict(rawcatCQ[modkey][dwellkey])
                        
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
                                
                                kqkmerged = tptools.merge_stp_dipole_cats_bypos(
                                        rawcatCQ[modkey].copy(),
                                        dwellkeys[1:],dwellkeys[0])
                                
                                cols2drop = []
                                for dwellkey in dwellkeys:
                                    cols2drop += ['X_%s' % dwellkey,'S_%s' % dwellkey]
                                kqkmerged.drop(cols2drop,axis=1)
                                
    
                                amplitudes = kqkmerged[Ampcols]
                                
    
                                Pc, tau = tptools.batch_fit_PcTau_stp(amplitudes,dwells,Nshuffles)
                                
                                
                                kqkmerged['Pc'] = pd.Series(Pc,index=kqkmerged.index)
                                
                                kqkmerged['tau'] = pd.Series(tau,index=kqkmerged.index)
                                
                                
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
            
        
        # Produce Pc, tau heatmaps for each tp mode across CCD beam
        
        for modkey in modkeys:
            
            pltfig = self.figdict['TP02meta_%s' % modkey]
            
            pldata = OrderedDict(labelkeys=polarities.keys())
            
            ixcoltau = colnames.index('tau')
            ixcolPc = colnames.index('Pc')
            ixcolS = colnames.index('uS')
            
            for CCDk in CCDs:
                pldata[CCDk] = OrderedDict()
                for Q in allQuads:
                    pldata[CCDk][Q] = OrderedDict(x=OrderedDict(),
                                                  y=OrderedDict())
                    
                    logtau = np.log10(mergecat[CCDk][Q][modkey][:,ixcoltau].copy())                    
                    logPc = np.log10(mergecat[CCDk][Q][modkey][:,ixcolPc].copy())
                    
                    S = mergecat[CCDk][Q][modkey][:,ixcolS].copy()
                    
                    for pkey in polarities.keys():
                        ixselS = np.where(S == polarities[pkey])
                        
                        if len(ixselS[0])>0:
                            pldata[CCDk][Q]['x'][pkey] = logtau[ixselS].copy()
                            pldata[CCDk][Q]['y'][pkey] = logPc[ixselS].copy()
                        else:
                            pldata[CCDk][Q]['x'][pkey] = []
                            pldata[CCDk][Q]['y'][pkey] = []
            
            if onTests:
                for CCDk in ['CCD2','CCD3']:
                    pldata[CCDk] = pldata['CCD1'].copy()
            
            pltfig[1]['data'] = pldata.copy()                    

                  
        if self.report is not None:
            Mfigkeys = ['TP02meta_%s' % mkey for mkey in modkeys]
            self.addFigures_ST(figkeys=Mfigkeys,
                           dobuilddata=False)
        
