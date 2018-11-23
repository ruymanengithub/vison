#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: BIAS0X

Bias-structure/RON analysis script

Created on Tue Aug 29 16:53:40 2017

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


from vison.pipe.task import HKKeys
from vison.support import context
from vison.datamodel import scriptic as sc
from vison.datamodel import ccd
import B0Xaux
from DarkTask import DarkTask
from vison.datamodel import inputs, cdp
from vison.support import utils
from vison.support.files import cPickleRead

# END IMPORT

isthere = os.path.exists


BIAS01_commvalues = dict(program='CALCAMP', test='BIAS01',
                         flushes=7, siflsh=1, siflsh_p=500,
                         inisweep=1,
                         vstart=0, vend=2086,
                         toi_fl=143., toi_tp=1000., toi_ro=1000., toi_ch=1000.,
                         chinj=0,
                         rdmode='fwd_bas',
                         swellw=context.sumwell['fwd_bas'][0],
                         swelldly=context.sumwell['fwd_bas'][1],
                         s_tpump=0,
                         v_tpump=0,
                         exptime=0., shuttr=0, e_shuttr=0,
                         mirr_on=0,
                         wave=4,
                         motr_on=0,
                         source='flat',
                         comments='BIAS')

BIAS02_commvalues = BIAS01_commvalues.copy()
BIAS02_commvalues['test'] = 'BIAS02'
BIAS02_commvalues['rdmode'] = 'rwd_bas_vs'
BIAS02_commvalues['swellw'] = context.sumwell['rwd_bas_vs'][0]
BIAS02_commvalues['swelldly'] = context.sumwell['rwd_bas_vs'][1]


class BIAS0X_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('N', ([int], 'Number of Frame Acquisitions.')),
    ])))


class BIAS0X(DarkTask):
    """ """

    inputsclass = BIAS0X_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data), ('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]

        super(BIAS0X, self).__init__(inputs, log, drill, debug)
        self.name = 'BIAS0X'
        self.type = 'Simple'
                
        self.HKKeys = HKKeys        
        self.figdict = B0Xaux.B0Xfigs.copy()
        self.CDP_lib = B0Xaux.CDP_lib.copy()
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles',
                                       profiles='profiles', products='products')
        

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(N=25)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds BIAS0X script structure dictionary.

        ###:param N: integer, number of frames to acquire.
        :param diffvalues: dict, opt, differential values.
        :param elvis: char, ELVIS version.

        """

        N = self.inputs['N']
        BIAS0X_sdict = dict(col001=dict(frames=N, exptime=0))

        Ncols = len(BIAS0X_sdict.keys())
        BIAS0X_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        if self.inputs['test'] == 'BIAS01':
            commvalues.update(BIAS01_commvalues)
        elif self.inputs['test'] == 'BIAS02':
            commvalues.update(BIAS02_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        BIAS0X_sdict = sc.update_structdict(
            BIAS0X_sdict, commvalues, diffvalues)

        return BIAS0X_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(BIAS0X, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                   wavedkeys=wavedkeys)

    def prep_data(self):
        """

        BIAS0X: Preparation of data for further analysis.
        Calls task.prepare_images().

        Applies:
            offset subtraction
            cosmetics masking

        """
        super(BIAS0X, self).prepare_images(
            doExtract=True, doMask=True, doOffset=True)

    def basic_analysis(self):
        """ 

        BIAS0X: Basic analysis of data.

        **METACODE**

        ::

            f. e. ObsID:
               f.e.CCD:

                   load ccdobj of ObsID, CCD

                   with ccdobj, f.e.Q:
                       produce a 2D poly model of bias, save coefficients
                       produce average profile along rows
                       produce average profile along cols
                       # save 2D model and profiles in a pick file for each OBSID-CCD
                       measure and save RON after subtracting large scale structure

            plot RON vs. time f. each CCD and Q
            plot average profiles f. each CCD and Q (color coded by time)

        """

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', 
                Title='%s Extraction' % self.inputs['test'], 
                level=0)

        Cindices = copy.deepcopy(self.dd.mx['File_name'].indices)
        self.dd.initColumn('profiles1D_name', Cindices,
                           dtype='S100', valini='None')
        self.dd.initColumn('mods2D_name', Cindices,
                           dtype='S100', valini='None')

        DDindices = copy.deepcopy(self.dd.indices)

        self.dd.initColumn('RON', DDindices, dtype='float32', valini=0.)

        nObs, nCCD, nQuad = DDindices.shape
        Quads = DDindices[2].vals
        CCDs = DDindices.get_vals('CCD')

        # The "Hard"-work

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))

        profs1D2plot = dict()
        profs1D2plot['hor'] = OrderedDict()
        profs1D2plot['ver'] = OrderedDict()

        for CCDk in CCDs:
            for tag in ['hor', 'ver']:
                profs1D2plot[tag][CCDk] = OrderedDict()
            for Q in Quads:
                for tag in ['hor', 'ver']:
                    profs1D2plot[tag][CCDk][Q] = OrderedDict()
                    profs1D2plot[tag][CCDk][Q]['x'] = OrderedDict()
                    profs1D2plot[tag][CCDk][Q]['y'] = OrderedDict()
        
        prof_all = []

        if not self.drill:
            
            CDP_header['DATE'] = self.get_time_tag()

            ccdpicklespath = self.inputs['subpaths']['ccdpickles']
            profilespath = self.inputs['subpaths']['profiles']

            for iObs in range(nObs):

                OBSID = self.dd.mx['ObsID'][iObs]

                for jCCD, CCDk in enumerate(CCDs):

                    vstart = self.dd.mx['vstart'][iObs, jCCD]
                    vend = self.dd.mx['vend'][iObs, jCCD]

                    ccdobj_name = self.dd.mx['ccdobj_name'][iObs, jCCD]

                    #if ccdobj_name == 'None':
                    #    continue  # TESTS

                    fullccdobj_name = os.path.join(
                        ccdpicklespath, '%s.pick' % ccdobj_name)

                    ccdobj = copy.deepcopy(cPickleRead(fullccdobj_name))

                    iprofiles1D = cdp.CDP()
                    iprofiles1D.header = CDP_header.copy()
                    iprofiles1D.path = profilespath
                    iprofiles1D.data = OrderedDict()

                    iprofiles2D = cdp.CDP()
                    iprofiles2D.header = CDP_header.copy()
                    iprofiles2D.path = profilespath
                    iprofiles2D.data = OrderedDict()

                    for kQ, Q in enumerate(Quads):

                        iprofiles1D.data[Q] = OrderedDict()

                        # produce a 2D poly model of bias, [save coefficients?]
                        # measure and save RON after subtracting large scale structure

                        imgmod2D = ccdobj.get_region2Dmodel(Q=Q, area='img', kind='poly2D',
                                        pdegree=5, doFilter=False, doBin=True, binsize=300,
                                        vstart=vstart, vend=vend, canonical=True, extension=-1)
                        
                        qimg = ccdobj.get_quad(Q, canonical=True)
                        
                        onlyRONimg = qimg[ccdobj.prescan:-ccdobj.overscan,vstart:vend]-\
                                         imgmod2D.imgmodel
                        
                        
                        #onlyRONimg = ccdobj.get_quad(
                        #    Q, canonical=False) - mod2D.imgmodel
                        _RON = onlyRONimg.std()

                        self.dd.mx['RON'][iObs, jCCD, kQ] = _RON

                        iprofiles2D.data[Q] = OrderedDict()
                        iprofiles2D.data[Q]['polycoeffs'] = copy.deepcopy(
                            imgmod2D.polycoeffs)
                        iprofiles2D.data[Q]['polyfunc'] = copy.deepcopy(
                            imgmod2D.polyfunc)

                        # produce average profile along rows
                        
                        

                        hor1Dprof = ccdobj.get_1Dprofile(Q=Q, orient='hor', area='all', stacker='mean',
                                                         vstart=vstart, vend=vend)

                        # produce average profile along cols

                        ver1Dprof = ccdobj.get_1Dprofile(Q=Q, orient='ver', area='all', stacker='mean',
                                                         vstart=vstart, vend=vend)

                        # save profiles in locally for plotting

                        iprofiles1D.data[Q]['hor'] = copy.deepcopy(hor1Dprof)
                        iprofiles1D.data[Q]['ver'] = copy.deepcopy(ver1Dprof)
                        
                        _profs = dict(hor=hor1Dprof,ver=ver1Dprof)

                        for oritag in ['hor', 'ver']:
                            _profdata = _profs[oritag].data.copy()
                            xorder = np.argsort(_profdata['x'])
                            _y = _profdata['y'][xorder]
                            _x = np.arange(len(_y))
                            profs1D2plot[oritag][CCDk][Q]['y']['OBS%i' %
                                                                 OBSID] = _y.copy()
                            profs1D2plot[oritag][CCDk][Q]['x']['OBS%i' %
                                                                 OBSID] = _x.copy()
                            
                            prof_all += _y.tolist()
                            
                            

                    # save (2D model and) profiles in a pick file for each OBSID-CCD

                    iprofiles1D.rootname = 'profs1D_%i_%s_%s' % (
                        OBSID, CCDk, self.inputs['test'])
                    iprofiles2D.rootname = 'mod2D_%i_%s_%s' % (
                        OBSID, CCDk, self.inputs['test'])

                    #fprofilespickf = os.path.join(profilespath,profilespickf)

                    iprofiles1D.savetopickle()
                    iprofiles2D.savetopickle()

                    # cPickleDumpDictionary(profiles,fprofilespickf)

                    self.dd.mx['profiles1D_name'][iObs,
                                                  jCCD] = iprofiles1D.rootname
                    self.dd.mx['mods2D_name'][iObs,
                                              jCCD] = iprofiles2D.rootname

        # PLOTS
        # profiles 1D (Hor & Ver) x (CCD&Q)
        # histograms of RON per CCD&Q
        
        if len(prof_all) != 0:        
             ylim_1D = [np.percentile(prof_all,5)-20.,
                       np.percentile(prof_all,95)+20.]
        else:
            ylim_1D= [0, 2.**16]

        figkeys1 = ['B0Xbasic_prof1D_hor', 'B0Xbasic_prof1D_ver',
                    'B0Xbasic_histosRON']
        
        for tag in ['hor','ver']:    
            profs1D2plot[tag]['labelkeys'] = \
                        profs1D2plot[tag][CCDs[0]][Quads[0]]['x'].keys()
        
        self.figdict['B0Xbasic_prof1D_hor'][1]['data'] = profs1D2plot['hor']
        self.figdict['B0Xbasic_prof1D_hor'][1]['meta']['ylim'] = ylim_1D
        self.figdict['B0Xbasic_prof1D_ver'][1]['data'] = profs1D2plot['ver']
        self.figdict['B0Xbasic_prof1D_ver'][1]['meta']['ylim'] = ylim_1D
        
        def _load_histo(histos,data,label,ronrange):
           
           hist, bin_edges = np.histogram(data,bins=10,range=ronrange,
                                          normed=False)
           histos['x'][label] = bin_edges.copy()
           histos['y'][label] = hist.copy()
           
           return histos
        
        #medron = np.nanmedian(self.dd.mx['std_pre'][:])
        #minron = np.int((medron-1)/0.5)*0.5
        #maxron = np.int((medron+1)/0.5)*0.5
        #ronrange = [minron,maxron]        
        ronrange = [0.5,2.]
        
        RON_histos = OrderedDict()
        for jCCD, CCDk in enumerate(CCDs):
            RON_histos[CCDk] = OrderedDict()
            for kQ, Q in enumerate(Quads):
                RON_histos[CCDk][Q] = OrderedDict()
                Rh_cq = RON_histos[CCDk][Q]
                Rh_cq['x'] = OrderedDict()
                Rh_cq['y'] = OrderedDict()
                Rh_cq = _load_histo(Rh_cq,self.dd.mx['RON'][:,jCCD,kQ],'IMG', 
                                    ronrange)
                Rh_cq = _load_histo(Rh_cq,self.dd.mx['std_pre'][:,jCCD,kQ],'PRE', 
                                    ronrange)
                Rh_cq = _load_histo(Rh_cq,self.dd.mx['std_ove'][:,jCCD,kQ],'OVE', 
                                    ronrange)
        
        RON_histos['labelkeys'] = ['IMG','PRE','OVE']
        
        
        self.figdict['B0Xbasic_histosRON'][1]['data'] = RON_histos.copy()
        self.figdict['B0Xbasic_histosRON'][1]['meta']['xlim'] = ronrange
        
        if self.report is not None:
            self.addFigures_ST(figkeys=figkeys1, dobuilddata=False)

        # REPORTS

        # Table with median (Best Estimate) values of RON per CCD&Q
        
                            
        RONdct = OrderedDict()
        RONdct['RON_IMG'] = self.dd.mx['RON'][:].copy()
        RONdct['RON_PRE'] = self.dd.mx['std_pre'][:].copy()
        RONdct['RON_OVE'] = self.dd.mx['std_ove'][:].copy()
        
        ron_cdp = self.CDP_lib['RON']
        ron_cdp.path = self.inputs['subpaths']['products']
        ron_cdp.ingest_inputs(mx_dct=RONdct.copy(),
                              CCDs=CCDs,
                              Quads=Quads,
                              meta=dict(),
                              header=CDP_header.copy())
        
        ron_cdp.init_wb_and_fillAll(header_title='%s: RON' % self.inputs['test'])
        self.save_CDP(ron_cdp)
        self.pack_CDP_to_dd(ron_cdp, 'RON_CDP')

        if self.report is not None:
            
            beRONtex = ron_cdp.get_textable(sheet='RON_IMG', 
                                            caption='%s: RON (quadrant-img, minus model)' % \
                                                             self.inputs['test'],
                                            longtable=False, fitwidth=True)
            self.report.add_Text(beRONtex)
            
            bePRERONtex = ron_cdp.get_textable(sheet='RON_PRE', 
                                               caption='%s: RON (pre-scan)' %\
                                                    self.inputs['test'],
                                               longtable=False,fitwidth=True)
            self.report.add_Text(bePRERONtex)
            
            beOVERONtex = ron_cdp.get_textable(sheet='RON_OVE', 
                            caption='%s: RON (over-scan)' % self.inputs['test'],
                            longtable=False, fitwidth=True)
            self.report.add_Text(beOVERONtex)
            
        # OFFSETS

        OFFdct = OrderedDict()
        OFFdct['OFF_PRE'] = self.dd.mx['offset_pre'][:].copy()
        OFFdct['OFF_IMG'] = self.dd.mx['offset_img'][:].copy()
        OFFdct['OFF_OVE'] = self.dd.mx['offset_ove'][:].copy()
        
        off_cdp = self.CDP_lib['OFF']
        off_cdp.path = self.inputs['subpaths']['products']
        off_cdp.ingest_inputs(mx_dct=OFFdct.copy(),
                              CCDs=CCDs,
                              Quads=Quads,
                              meta=dict(),
                              header=CDP_header.copy())
        
        off_cdp.init_wb_and_fillAll(header_title='%s: OFFSETS' % self.inputs['test'])

        self.save_CDP(off_cdp)
        self.pack_CDP_to_dd(off_cdp, 'OFF_CDP')

        if self.report is not None:
            PREOFFtex = off_cdp.get_textable(sheet='OFF_PRE', 
                                caption='%s: Offsets, pre-scan.' % self.inputs['test'],
                                            longtable=False, fitwidth=True)
            self.report.add_Text(PREOFFtex)
            
            IMGOFFtex = off_cdp.get_textable(sheet='OFF_IMG', 
                            caption='%s: Offsets, image area.' % self.inputs['test'],
                                               longtable=False, fitwidth=True)
            self.report.add_Text(IMGOFFtex)
            
            OVEOFFtex = off_cdp.get_textable(sheet='OFF_OVE', 
                            caption='%s: Offsets, over-scan.' % self.inputs['test'],
                                               longtable=False, fitwidth=True)
            self.report.add_Text(OVEOFFtex)



    def meta_analysis(self):
        """

        **METACODE**

        ::

            f. each CCD:
               stack all ObsIDs to produce Master Bias
               f. e. Q:    
                   measure average profile along rows
                   measure average profile along cols
            plot average profiles of Master Bias(s) f. each CCD,Q
            (produce table(s) with summary of results, include in report)
            save Master Bias(s) (3 images) to FITS CDPs
            show Master Bias(s) (3 images) in report
            save name of MasterBias(s) CDPs to DataDict, report
        
        """
        
        if self.report is not None:
            self.report.add_Section(
                keyword='meta', 
                Title='%s Meta-Analysis' % self.inputs['test'], 
                level=0)
        
        DDindices = copy.deepcopy(self.dd.indices)
        
        nObs, nCCD, nQuad = DDindices.shape
        Quads = DDindices[2].vals
        CCDs = DDindices.get_vals('CCD')
        
        # The "Hard"-work

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        # 1D Profiles of Master Bias
        
        profs1D2plot = OrderedDict()
        profs1D2plot['hor'] = OrderedDict()
        profs1D2plot['ver'] = OrderedDict()
        
        for CCDk in CCDs:
            for tag in ['hor', 'ver']:
                profs1D2plot[tag][CCDk] = OrderedDict()
            for Q in Quads:
                for tag in ['hor', 'ver']:
                    profs1D2plot[tag][CCDk][Q] = OrderedDict()
                    profs1D2plot[tag][CCDk][Q]['x'] = np.arange(100,dtype='float32')
                    profs1D2plot[tag][CCDk][Q]['y'] = np.zeros(100,dtype='float32')
        
        
        # Master Bias Data to be plotted
        
        MB_2PLOT = OrderedDict()
        
        for CCDk in CCDs:
            MB_2PLOT[CCDk] = OrderedDict()
            for Q in Quads:
                MB_2PLOT[CCDk][Q] = OrderedDict()
        MB_p5s = []
        MB_p95s = []

        if not self.drill:
            
            CDP_header['DATE'] = self.get_time_tag()
            ccdpicklespath = self.inputs['subpaths']['ccdpickles']
            profilespath = self.inputs['subpaths']['profiles']
            prodspath = self.inputs['subpaths']['products']
            
            
            def _pack_profs(CQdict,prof):
                """ """
                _x = prof.data['x'].copy()
                xorder = np.argsort(prof.data['x'])
                _y = prof.data['y'][xorder].copy()
                _x = np.arange(len(_y))
                
                CQdict['x'] = _x.copy()
                CQdict['y'] = _y.copy()
                
                return CQdict
            
            vfullinpath_adder = utils.get_path_decorator(ccdpicklespath)
            
            for jCCD, CCDk in enumerate(CCDs):
                
                sn_ccd = self.inputs['diffvalues']['sn_%s' % CCDk.lower()]                
                vstart = self.dd.mx['vstart'][0, jCCD]
                vend = self.dd.mx['vend'][0, jCCD]
                
                ccdobj_names = self.dd.mx['ccdobj_name'][:, jCCD]
                fullccdobj_names = vfullinpath_adder(ccdobj_names,extension='pick')
                
                
                devel = False # TESTS
                
                if not devel:
                
                    ccdobjList = [copy.deepcopy(cPickleRead(item)) \
                                    for item in fullccdobj_names]
                    
                    ccdpile = ccd.CCDPile(ccdobjList=ccdobjList,
                                          withpover=ccdobjList[0].withpover)
                    
                    stackimg, stackstd = ccdpile.stack(method='median',dostd=True)
                
                else:
                    
                    stackimg = np.zeros((4238,4172),dtype='float32')
                    stackstd = np.ones((4238,4172),dtype='float32')
                
                if isinstance(stackimg,np.ma.MaskedArray):
                    mask = stackimg.mask.copy()
                    stackimg = stackimg.data.copy()
                    stackstd = stackstd.data.copy()
                else:
                    mask = np.zeros_like(stackimg)
                    
                
                #stackccdobj = ccd.CCD(withpover=ccdpile.withpover)
                stackccdobj = ccd.CCD(withpover=True) # TESTS
                
                stackccdobj.add_extension(stackimg, label='STACK')
                stackccdobj.add_extension(stackstd, label='STD')
                stackccdobj.get_mask(mask)
                
                for kQ, Q in enumerate(Quads):
                    
                    hor1Dprof = stackccdobj.get_1Dprofile(Q=Q, orient='hor',
                                    area='all',stacker='mean',
                                    vstart=vstart,vend=vend,extension=0)
                    
                    profs1D2plot['hor'][CCDk][Q] = _pack_profs(profs1D2plot['hor'][CCDk][Q],hor1Dprof)
                    
                    ver1Dprof = stackccdobj.get_1Dprofile(Q=Q, orient='ver',
                                    area='all',stacker='mean',
                                    vstart=vstart,vend=vend,extension=0)
                    
                    profs1D2plot['ver'][CCDk][Q] = _pack_profs(profs1D2plot['ver'][CCDk][Q],ver1Dprof)
                    
                
                # SAVING Master Bias to a CDP
                
                mb_data = OrderedDict()
                mb_data['STACK'] = stackimg.copy()
                mb_data['STD'] = stackstd.copy()
                mb_data['MASK'] = mask.astype('int32').copy()
                mb_data['labels'] = ['STACK','STD','MASK']
                mb_meta = OrderedDict()
                mb_meta['CCD_SN']= sn_ccd
    
                masterbiascdp = cdp.CCD_CDP(ID=self.ID,
                                      BLOCKID=self.BLOCKID,
                                      CHAMBER=self.CHAMBER)
    
                masterbiascdp.ingest_inputs(data=mb_data, meta=mb_meta, header=CDP_header)
                
                masterbiascdp.path = prodspath
                masterbiascdp.rootname = 'EUC_MASTERBIAS_%s_SN_%s_%s' % \
                                           (CCDk, sn_ccd, self.inputs['BLOCKID'])
                
                
                for Q in Quads:
                    qdata = masterbiascdp.ccdobj.get_quad(Q,canonical=False,extension=1).copy()
                    MB_p5s.append(np.percentile(qdata,5))
                    MB_p95s.append(np.percentile(qdata,95))
                    MB_2PLOT[CCDk][Q]['img'] = qdata.transpose()
                
                self.save_CDP(masterbiascdp)
                self.pack_CDP_to_dd(masterbiascdp,'MASTERBIAS_%s' % CCDk)
                
        
        # PLOTTING 1D PROFILES OF MASTER BIAS
        
        self.figdict['B0Xmeta_prof1D_hor'][1]['data'] = profs1D2plot['hor'].copy()
        self.figdict['B0Xmeta_prof1D_ver'][1]['data'] = profs1D2plot['ver'].copy()
        
        if self.report is not None:
            self.addFigures_ST(figkeys=['B0Xmeta_prof1D_hor',
                                        'B0Xmeta_prof1D_ver'],
                               dobuilddata=False)
        
        # SAVING 1D PROFILES OF MASTER BIAS as a CDP
        
        MB_profiles_cdp = cdp.CDP()
        MB_profiles_cdp.header = CDP_header.copy()
        MB_profiles_cdp.path = profilespath
        MB_profiles_cdp.data = profs1D2plot.copy()
        
        self.save_CDP(MB_profiles_cdp)
        self.pack_CDP_to_dd(MB_profiles_cdp, 'MB_PROFILES')
        
        
        # DISPLAYING THE MASTER BIAS FRAMES
        
        self.figdict['B0Xmeta_MasterBias_2D'][1]['data'] = MB_2PLOT.copy()
        
        # UPDATING scaling based on data
        
        if len(MB_p5s)>0:
            normfunction = Normalize(vmin=np.min(MB_p5s),vmax=np.max(MB_p95s),clip=False)
        else:
            normfunction = False
        
        self.figdict['B0Xmeta_MasterBias_2D'][1]['meta']['corekwargs']['norm'] = normfunction
        
        if self.report is not None:
            self.addFigures_ST(figkeys=['B0Xmeta_MasterBias_2D'],
                               dobuilddata=False)


class Test(unittest.TestCase):
    """
    Unit tests for the BIAS0X class.
    """

    def setUp(self):

        inputs = dict(test='BIAS01')
        self.b01 = BIAS0X(inputs, log=None, drill=True, debug=False)

    def test_check_data(self):
        """

        :return: None
        """
        self.b01.check_data()


if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(Test)
    unittest.TextTestRunner(verbosity=3).run(suite)
