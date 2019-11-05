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
from skimage import exposure

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
from vison.datamodel import inputs as inputsmod
from vison.support.files import cPickleRead
import TP01aux
from vison.pump import tptools
# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

IDL = 11.
IDH = 18.
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
                       exptime=0., shuttr=0, e_shuttr=0,
                       mirr_on=0,
                       wave=4,
                       motr_on=0,
                       source='flat',
                       comments='')


class TP01_inputs(inputsmod.Inputs):
    manifesto = inputsmod.CommonTaskInputs.copy()
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
    #contrast_threshold = 0.01
    contrast_thresholdfactor = 5.

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('injection', self.charact_injection),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('debugtask', self.debugtask),
                         ('meta', self.meta_analysis)]
        self.commvalues = TP01_commvalues.copy()
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

        self.commvalues['v_tp_cnt'] = Nshuffles_V

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
        commvalues.update(self.commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
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
                                         doMask=True,  # False ON TESTS!
                                         doOffset=True,
                                         doBias=False,
                                         doFF=False)

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

        testname = self.inputs['test']

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='%s Extraction' % testname, level=0)

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
                        vend = self.dd.mx['vend'][ix, jCCD]

                        ioutf = 'TP01_rawmap_%i_IDDLY_%i_ROE1_%s' % (
                            ObsID, id_dly, CCDk)

                        iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ix, jCCD]

                        try:
                            iccdobj = cPickleRead(
                                os.path.join(ccdpicklespath, iccdobj_f))
                            irawmap = copy.deepcopy(iccdobj)

                            irawmap = tptools.gen_raw_dpmap_vtpump(
                                irawmap, Navgrows=-1, vstart=vstart, vend=vend)

                            irawmap.writeto(os.path.join(productspath,
                                                         '%s.fits' % ioutf))

                            self.dd.mx['dipoles_raw'][ix, jCCD] = ioutf

                        except BaseException:  # TESTS
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

        testname = self.inputs['test']

        if self.report is not None:
            self.report.add_Section(
                keyword='basic', Title='%s Basic Analysis' % testname, level=0)

        #threshold = self.contrast_threshold
        threshfactor = self.contrast_thresholdfactor

        #CCDhalves = ['top','bottom']
        _Quads_dict = dict(bottom=['E', 'F'],
                           top=['G', 'H'])

        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        allQuads = DDindices.get_vals('Quad')

        productspath = self.inputs['subpaths']['products']

        id_dlys = np.unique(self.dd.mx['id_dly'][:, 0])

        mods = np.unique(self.dd.mx['v_tp_mod'][:, 0])
        modkeys = ['m%i' % item for item in mods]

        tois = np.unique(self.dd.mx['toi_tp'][:, 0])
        toikeys = ['u%04i' % item for item in tois]

        # initialisation

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        masterdict = OrderedDict()
        threshold_dict = OrderedDict()

        for CCDk in CCDs:
            masterdict[CCDk] = OrderedDict()
            threshold_dict[CCDk] = OrderedDict()
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

                        ixsel = np.where((self.dd.mx['id_dly'][:, 0] == id_dly) & (
                            self.dd.mx['v_tpump'][:, 0] != 0))

                        for ix in ixsel[0]:

                            ObsID = self.dd.mx['ObsID'][ix]
                            vstart = self.dd.mx['vstart'][ix, jCCD]
                            vend = self.dd.mx['vend'][ix, jCCD]
                            toi_ch = float(self.dd.mx['toi_ch'][ix, jCCD])
                            v_tp_mod = self.dd.mx['v_tp_mod'][ix, jCCD]
                            toi_tp = self.dd.mx['toi_tp'][ix, jCCD]

                            modkey = 'm%i' % v_tp_mod
                            toikey = 'u%04i' % toi_tp

                            if np.isclose(id_dly / toi_ch, 1.5):
                                CCDhalf = 'top'
                            elif np.isclose(id_dly / toi_ch, 2.5):
                                CCDhalf = 'bottom'

                            Quads = _Quads_dict[CCDhalf]

                            imapf = 'TP01_rawmap_%i_IDDLY_%i_ROE1_%s.fits' % (
                                ObsID, id_dly, CCDk)

                            imapccdobj = ccd.CCD(os.path.join(productspath, imapf))

                            print('OBSID=%s, %s' % (ObsID, imapf))

                            for iQ, Q in enumerate(Quads):

                                _med = self.dd.products['chinj_mx'][CCDk][Q]
                                _sig = self.dd.products['chinjnoise_mx'][CCDk][Q]

                                threshold = threshfactor * _sig / _med

                                threshold_dict[CCDk][Q] = threshold

                                idd = tptools.find_dipoles_vtpump(imapccdobj, threshold,
                                                                  Q, vstart=vstart, vend=vend,
                                                                  extension=-1)

                                masterdict[CCDk][Q][modkey][toikey] = idd.copy()

            df = tptools._aggregate_CQMT(masterdict, CCDs, allQuads, modkeys, toikeys, 'toi')

            self.dd.products['master_df'] = df.copy()

            #summaryMean = df.groupby(level=['CCD','Q','mod']).mean()
            #summaryTot = df.groupby(level=['CCD','Q','mod']).sum()

            summary = df.groupby(level=['CCD', 'Q', 'mod']).mean().copy()
            #summary['R'] = summaryMean['R'].copy()
            #summary['A'] = summaryMean['A'].copy()
            summary.columns = ['<N>', '<R>', '<A>']

            if self.report is not None:
                summtxtreport = tptools._get_txt(summary)
                self.report.add_Text(["Aggregated Dipole Statistics: $<$Number$>$, $<$Ratio N/S$>$, $<$Amplitude$>$",
                                      "\nAverages accross toi's"])
                self.report.add_Text(summtxtreport)

            # Saving the MasterCat

            MCmeta = OrderedDict()
            MCmeta['THRESHOLDs'] = threshold_dict
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
                self.pack_CDP_to_dd(kmastercat, 'MASTERCAT_%s' % CCDk)

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

        testname = self.inputs['test']

        if self.report is not None:
            self.report.add_Section(
                keyword='meta', Title='%s Meta Analysis' % testname, level=0)

        DDindices = copy.deepcopy(self.dd.indices)
        CCDs = DDindices.get_vals('CCD')
        allQuads = DDindices.get_vals('Quad')

        productspath = self.inputs['subpaths']['products']
        Nshuffles = self.inputs['Nshuffles_V']

        # initialisation

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        mods = np.unique(self.dd.mx['v_tp_mod'][:, 0])
        modkeys = ['m%i' % item for item in mods]

        tois = np.unique(self.dd.mx['toi_tp'][:, 0])
        toikeys = ['u%04i' % item for item in tois]
        Ampcols = ['A_%s' % toikey for toikey in toikeys]

        onTests = False

#        if not onTests:

        if not self.drill:

            mergecat = OrderedDict()

            print('\nMerging Dipole Catalogs...\n')

            for jCCD, CCDk in enumerate(CCDs):
                #            for jCCD, CCDk in enumerate(['CCD1']):

                mastercatpick = self.dd.products['MASTERCAT_%s' % CCDk]

                masterdata = cPickleRead(mastercatpick)['data'].copy()

                mergecat[CCDk] = OrderedDict()

                for iQ, Q in enumerate(allQuads):
                    #                for iQ, Q in enumerate(['E']):

                    mergecat[CCDk][Q] = OrderedDict()

                    rawcatCQ = masterdata[Q].copy()

                    # Pandizing the toi catalogs

                    for modkey in modkeys:

                        for toikey in toikeys:

                            if onTests:
                                rawcatCQ[modkey][toikey] = tptools._thin_down(
                                    rawcatCQ[modkey][toikey], 1000)

                            rawcatCQ[modkey][toikey] =\
                                pd.DataFrame.from_dict(rawcatCQ[modkey][toikey])

                    # Merging the toi catalogs

                    for modkey in modkeys:

                        print('%s%s, %s...' % (CCDk, Q, modkey))

                        kqkmerged = tptools.merge_vtp_dipole_cats_bypos(
                            rawcatCQ[modkey].copy(),
                            toikeys[1:], toikeys[0])

                        cols2drop = []
                        for toikey in toikeys:
                            cols2drop += ['X_%s' % toikey, 'Y_%s' % toikey, 'S_%s' % toikey]
                        kqkmerged.drop(cols2drop, axis=1)

                        amplitudes = kqkmerged[Ampcols]

                        scaled_Pc, tau = tptools.batch_fit_PcTau_vtp(amplitudes, tois, Nshuffles)

                        Pc = scaled_Pc * self.dd.products['chinj_mx'][CCDk][Q]

                        kqkmerged['Pc'] = pd.Series(Pc, index=kqkmerged.index)

                        kqkmerged['tau'] = pd.Series(tau, index=kqkmerged.index)

                        mergecat[CCDk][Q][modkey] = kqkmerged

        # Store output catalog(s) as a CDP

#        if not onTests:

        for CCDk in CCDs:
            #        for CCDk in ['CCD1']:

            kmergedata = mergecat[CCDk].copy()
            colnames = kmergedata[allQuads[0]][modkeys[0]].keys().tolist()

            for Q in allQuads:
                #            for Q in ['E']:
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
            self.pack_CDP_to_dd(kmergecat, 'MERGEDCAT_%s' % CCDk)

        # SUMMARY MX/TABLE

        ixtau = colnames.index('tau')
        ixPc = colnames.index('Pc')

        sumtable = OrderedDict()

        for CCDk in CCDs:
            sumtable[CCDk] = OrderedDict()
            for Q in allQuads:
                sumtable[CCDk][Q] = OrderedDict()
                for modkey in modkeys:

                    tau = mergecat[CCDk][Q][modkey][:, ixtau].copy()
                    Pc = mergecat[CCDk][Q][modkey][:, ixPc].copy()
                    ixnonan = np.where(~np.isnan(tau) & ~np.isnan(Pc))
                    if len(ixnonan[0]) > 0:
                        avtau = np.average(tau[ixnonan], weights=Pc[ixnonan])
                    else:
                        avtau = np.nan

                    sumtable[CCDk][Q][modkey] = [len(Pc[ixnonan]), avtau]

        reform = {(level1_key, level2_key, level3_key): value
                  for level1_key, level2_dict in sumtable.items()
                  for level2_key, level3_dict in level2_dict.items()
                  for level3_key, value in level3_dict.items()}

        mergeL1df = pd.DataFrame(reform).T
        colnames = dict()
        for i, c in enumerate(['N', '<tau>']):
            colnames[i] = c
        mergeL1df.rename(columns=colnames, inplace=True)
        names = ['CCD', 'Q', 'mod']
        mergeL1df.index.set_names(names, inplace=True)

        mergeL2df = mergeL1df.groupby(level=['CCD', 'Q']).sum().copy()
        tausL2df = mergeL1df.groupby(level=['CCD', 'Q']).mean().copy()
        mergeL2df['<tau>'] = tausL2df['<tau>'].copy()

        self.dd.products['mergedL1_df'] = mergeL1df.copy()
        self.dd.products['mergedL2_df'] = mergeL2df.copy()

        if self.report is not None:

            txtreportML1 = tptools._get_txt(mergeL1df)
            self.report.add_Text(["Dipole Tau Statistics: Number, $<$tau [us]$>$",
                                  "\nModel-fit dipoles only."])
            self.report.add_Text(txtreportML1)

            txtreportML2 = tptools._get_txt(mergeL2df)
            self.report.add_Text(["Dipole Tau Statistics (mode aggregated): Number, $<$tau [us]$>$",
                                  "\nModel-fit dipoles only, aggregated across pumping modes."])
            self.report.add_Text(txtreportML2)

        # Produce Pc, tau heatmaps for each tp mode across CCD beam

        for modkey in modkeys:

            pltfig = self.figdict['TP01meta_%s' % modkey]

            pldata = OrderedDict()

            #ixcoltau = colnames.index('tau')
            #ixcolPc = colnames.index('Pc')
            #ixcolS = colnames.index('uS')

            for CCDk in CCDs:
                pldata[CCDk] = OrderedDict()
                for Q in allQuads:

                    pldata[CCDk][Q] = OrderedDict()

                    # if CCDk!= 'CCD1' or Q != 'E':
                    #    pldata[CCDk][Q]['img'] = np.zeros((15,15))
                    #    continue

                    logtau = np.log10(mergecat[CCDk][Q][modkey][:, ixtau].copy())
                    logPc = np.log10(mergecat[CCDk][Q][modkey][:, ixPc].copy())

                    Heatmap, xedges, yedges = np.histogram2d(logtau, logPc, bins=(15, 15),
                                                             range=[[1.5, 5], [-3., 0.5]])

                    Heatmap /= np.nanmax(Heatmap)

                    if CCDk == CCDs[0] and Q == allQuads[0]:
                        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

                    # HeatmapPeaks.append(np.nanmax(Heatmap))

                    heqHeatmap = exposure.equalize_hist(Heatmap, nbins=256)

                    pldata[CCDk][Q]['img'] = heqHeatmap.copy()

            # if onTests:
            #    for CCDk in ['CCD2','CCD3']:
            #        pldata[CCDk] = pldata['CCD1'].copy()

            pltfig[1]['data'] = pldata.copy()

            #normfunction = Normalize(vmin=1,vmax=np.mean(HeatmapPeaks))

            #pltfig[1]['meta']['corekwargs']['norm'] = normfunction
            pltfig[1]['meta']['corekwargs']['extent'] = extent

        if self.report is not None:
            Mfigkeys = ['TP01meta_%s' % mkey for mkey in modkeys]
            self.addFigures_ST(figkeys=Mfigkeys,
                               dobuilddata=False)

        self.canbecleaned = True

    def debugtask(self):

        if self.report is not None:
            self.report.add_Section(
                keyword='debug', Title='TP01 DEBUG', level=0)

        CCDs = ['CCD1', 'CCD2', 'CCD3']

        masterdict = OrderedDict()
        for CCDk in CCDs:
            kmastercatpick = self.dd.products['MASTERCAT_%s' % CCDk]
            masterdict[CCDk] = cPickleRead(kmastercatpick)['data'].copy()

        allQuads = ['E', 'F', 'G', 'H']
        modkeys = ['m123', 'm234', 'm341', 'm412']
        toikeys = ['u0050', 'u0200', 'u1000', 'u2000', 'u4000', 'u8000']

        dfL1 = tptools._aggregate_CQMT(masterdict, CCDs, allQuads, modkeys, toikeys, 'toi')

        summaryCat = dfL1.groupby(level=['CCD', 'Q', 'mod']).mean().copy()
        summaryCat.columns = ['<N>', '<R>', '<A>']

        if self.report is not None:

            txtreport = tptools._get_txt(summaryCat)
            self.report.add_Text(["Aggregated Dipole Statistics: $<$Number$>$, $<$Ratio N/S$>$, $<$Amplitude$>$",
                                  "\nAverages accross toi's"])
            self.report.add_Text(txtreport)

        # MERGED

        rawmergedict = OrderedDict()
        for ik, CCDk in enumerate(CCDs):
            kmergedpick = self.dd.products['MERGEDCAT_%s' % CCDk]
            kmerged = cPickleRead(kmergedpick)
            if ik == 0:
                mcolnames = kmerged['meta']['colnames']
            rawmergedict[CCDk] = kmerged['data'].copy()

        def get_vals(arrayCQM, colnames):
            """ """
            Pc = arrayCQM[:, colnames.index('Pc')]
            tau = arrayCQM[:, colnames.index('tau')]
            ixnonan = np.where(~np.isnan(Pc) & ~np.isnan(tau))
            Pc = Pc[ixnonan]
            tau = tau[ixnonan]
            if len(ixnonan[0]) > 0:
                values = [len(Pc), np.average(tau, weights=Pc)]
            else:
                values[0, np.nan]

            return values

        mergedict = OrderedDict()
        for CCDk in CCDs:
            mergedict[CCDk] = OrderedDict()
            for Q in allQuads:
                mergedict[CCDk][Q] = OrderedDict()
                for mk in modkeys:
                    mergedict[CCDk][Q][mk] = get_vals(rawmergedict[CCDk][Q][mk], mcolnames)

        reform = {(level1_key, level2_key, level3_key): value
                  for level1_key, level2_dict in mergedict.items()
                  for level2_key, level3_dict in level2_dict.items()
                  for level3_key, value in level3_dict.items()}

        mergeL1df = pd.DataFrame(reform).T
        colnames = dict()
        for i, c in enumerate(['N', '<tau>']):
            colnames[i] = c
        mergeL1df.rename(columns=colnames, inplace=True)
        names = ['CCD', 'Q', 'mod']
        mergeL1df.index.set_names(names, inplace=True)

        mergeL2df = mergeL1df.groupby(level=['CCD', 'Q']).sum().copy()
        tausL2df = mergeL1df.groupby(level=['CCD', 'Q']).mean().copy()
        mergeL2df['<tau>'] = tausL2df['<tau>'].copy()

        if self.report is not None:

            txtreportML1 = tptools._get_txt(mergeL1df)
            self.report.add_Text(["Dipole Tau Statistics: Number, $<$tau [us]$>$",
                                  "\nModel-fit dipoles only."])
            self.report.add_Text(txtreportML1)

            txtreportML2 = tptools._get_txt(mergeL2df)
            self.report.add_Text(["Dipole Tau Statistics (mode aggregated): Number, $<$tau [us]$>$",
                                  "\nModel-fit dipoles only, aggregated across pumping modes."])
            self.report.add_Text(txtreportML2)
