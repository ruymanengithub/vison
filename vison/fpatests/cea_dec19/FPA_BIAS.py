# -*- coding: utf-8 -*-
"""

TEST: BIAS (Multiple Readout Modes)

Created on Wed Nov 13 16:10:10 2019

@author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import sys
from collections import OrderedDict
import pandas as pd
import string as st

from vison.support import utils, files
from vison.image import bits
from vison.other import MOT_FFaux
from vison.fpatests import fpatask
from vison.datamodel import inputs
from vison.datamodel import cdp as cdpmod
from vison.fpatests.cea_dec19 import BIAS_aux
from vison.datamodel import cdp

# END IMPORT

class FPA_BIAS_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('readmode', ([str], 'Readout Mode')),
        ('temperature', ([str], 'Temperature at which the test is performed')),        
    ])))

class FPA_BIAS(fpatask.FpaTask):

    inputsclass = FPA_BIAS_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.subtasks = [('check', self.check_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis),
                         ('appendix', self.appendix)]

        super(FPA_BIAS, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                                       cleanafter=cleanafter)

        self.name = 'FPA_BIAS'
        self.type = 'Simple'
        self.readmode = self.inputs['readmode']
        self.temperature = self.inputs['temperature']

        self.figdict = BIAS_aux.get_Bfigs(self.readmode,self.temperature)
        self.CDP_lib = BIAS_aux.get_CDP_lib(self.readmode,self.temperature)

        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')
        if self.inputs['inCDPs'] is not None:
            if 'gain' in self.inputs['inCDPs']:
                allgains = files.cPickleRead(self.inputs['inCDPs']['gain'])

                self.incdps['GAIN'] = OrderedDict()
                
                for block in self.fpa.flight_blocks:
                    self.incdps['GAIN'][block] = allgains[block]['PTC01'].copy()


    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(
            temperature='cold',
            readmode='fwd',
            preprocessing=dict(
                offsetkwargs=dict(
                    ignore_pover=True,
                    trimscan=[25, 5],
                    method='median',
                    extension=-1,
                    scan='pre'
                    )))

    def load_references(self):
        """ """

        readmode = self.inputs['readmode']
        temperature = self.inputs['temperature']
        
        testkey = '%s_%s' % (readmode.upper(),temperature.upper())

        lookup_OFF_refkeys = dict(
            RWDVS_WARM = ('offsets_rwdvs_warm','ove'),
            RWDVS_COLD = ('offsets_rwdvs_warm','ove'),
            RWDV_WARM = ('offsets_rwdv_warm','ove'),
            RWDV_COLD = ('offsets_fwd_cold','ove'),
            FWD_WARM = ('offsets_fwd_cold','ove'),
            FWD_COLD= ('offsets_fwd_cold','ove'),
            )

        refoffkey, offreg = lookup_OFF_refkeys[testkey]

        refOFF_incdp = cdpmod.Json_CDP()
        refOFF_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refoffkey])
        

        def _get_ref_OFFs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][offreg]
        
        
        RefOFFsMap = self.get_FPAMAP(refOFF_incdp.data.copy(),
                                        extractor=_get_ref_OFFs_MAP)

        self.dd.products['REF_OFFs'] = RefOFFsMap.copy()

        refronkey = 'rons'
        ronreg = 'img'

        refRON_incdp = cdpmod.Json_CDP()
        refRON_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refronkey])
        
        def _get_ref_RONs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][ronreg]
        
        
        RefRONsMap = self.get_FPAMAP(refRON_incdp.data.copy(),
                                        extractor=_get_ref_RONs_MAP)

        self.dd.products['REF_RONs'] = RefRONsMap.copy()



    def _basic_onLE1(self, **kwargs):
    	""" """

    	CCDID = kwargs['CCDID']
        LE1 = kwargs['LE1']
        vstart = kwargs['vstart']
        vend = kwargs['vend']

        kccdobj = LE1.get_ccdobj(CCDID)
        # print(CCDID)
        debug = kwargs['debug']        

        trimscans = dict(pre=[25, 5],
                         img=[5, 5],
                         ove=[5, 5])

        if not debug:

            for Q in self.Quads:

                # produce average profile along rows
                hor1Dprof = kccdobj.get_1Dprofile(
                                Q=Q, orient='hor', area='all', stacker='mean', vstart=vstart, vend=vend)

                hQdata = hor1Dprof.data.copy()
                hx = hQdata['x'].copy()
                hy = hQdata['y'].copy()
                hx -= hx.min()
                if Q in ['F','G']:
                    hx = kccdobj.NAXIS1/2-hx

                ixsort = np.argsort(hx)
                hx = hx[ixsort]
                hy = hy[ixsort]

                hQdata = dict(x=hx.copy(),y=hy.copy())

                self.dd.products['profiles_H_1D'][CCDID][Q] = hQdata.copy()


                for reg in ['pre','img','ove']:

                    # produce average profile along cols
                    ver1Dprof = kccdobj.get_1Dprofile(
                                Q=Q, orient='ver', area=reg, stacker='mean', 
                                vstart=vstart, vend=vend)
                    
                    vQdata = ver1Dprof.data.copy()
                    vx = vQdata['x'].copy()
                    vy = vQdata['y'].copy()
                    vx -= vx.min()
                    if Q in ['E','F']:
                        vx = kccdobj.NAXIS2/2-vx

                    ixsort = np.argsort(vx)
                    vx = vx[ixsort]
                    vy = vy[ixsort]

                    vQdata = dict(x=vx.copy(),y=vy.copy())

                    self.dd.products['profiles_V_1D_%s' % reg][CCDID][Q] = vQdata.copy()


                for reg in ['pre','img','ove']:

                    stats = kccdobj.get_stats(
                                    Q,sector=reg,
                                    statkeys=[
                                        'median','std'],
                                    trimscan=trimscans[reg],
                                    ignore_pover=True,
                                    extension=-1,
                                    VSTART=vstart,
                                    VEND=vend)

                    self.dd.products['MED_%s' % reg.upper()][CCDID][Q] = stats[0]
                    self.dd.products['STD_%s' % reg.upper()][CCDID][Q] = stats[1]


    def basic_analysis(self):
        """
        To-Do:
            Extract average profiles along columns
            Extract average profiles along rows 
            Extract STDs and MEDs in pre/img/ove
            
        """

        if self.report is not None:
            self.report.add_Section(keyword='extract', Title='Extraction', level=0)

        iObs = 0
        vstart = 0
        vend = 2086

        #readmode = self.inputs['readmode']
        #temperature = self.inputs['temperature']

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        LE1fits = self.dd.mx['File_name'][iObs]

        fullLE1fits = os.path.join(self.inputs['datapath'], LE1fits)

        LE1 = self.load_LE1(fullLE1fits)

        # DISPLAY IMAGE

        RawImgDict = self.get_ImgDictfromLE1(LE1, doequalise=True)

        self.figdict['raw_img'][1]['data'] = RawImgDict
        self.figdict['raw_img'][1]['meta']['plotter'] = self.metacal.plot_ImgFPA

        if self.report is not None:
            self.addFigures_ST(figkeys=['raw_img'], dobuilddata=False)

        prodkeys = []
        for reg in ['pre','img','ove']:
            prodkeys.append('profiles_V_1D_%s' % reg)

        prodkeys += ['profiles_H_1D', 
            'STD_PRE', 'STD_IMG','STD_OVE',
            'MED_PRE', 'MED_IMG', 'MED_OVE',
            'DARK']

        for prodkey in prodkeys:

            self.dd.products[prodkey] = dict()

            for jY in range(1, self.NSLICES_FPA + 1):
                for iX in range(1, self.NCOLS_FPA + 1):
                    CCDID = 'C_%i%i' % (jY, iX)
                    self.dd.products[prodkey][CCDID] = dict()

        Bkwargs = dict(vstart=vstart,
                       vend=vend,
                       debug=False)

        self.iterate_over_CCDs_inLE1(LE1, FPA_BIAS._basic_onLE1, **Bkwargs)

        if self.report is not None:
            for prodkey in prodkeys:
                nprodkey = prodkey.replace('_','\_')
                self.report.add_Text('product: %s, all extracted!' % nprodkey)

        # Matrix of offsets (over-scan)

        off_tb_cdp = cdp.Tables_CDP()
        offcdpdict = dict(
            caption='Offsets in the over-scan region.',
            valformat='%.2f')

        def _getOFF_OVEval(self, Ckey, Q):
            return self.dd.products['MED_OVE'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getOFF_OVEval,
                                    cdp=off_tb_cdp,
                                    cdpdict=offcdpdict)


        # Save offsets to a json: GIVES PROBLEMS because of serialisation of numpy.floats!

        #OFF_hdr = OrderedDict()
        #OFF_hdr['title'] = 'OFFSETS'
        #OFF_hdr.update(CDP_header)

        #OFF_cdp = self.CDP_lib['OFF_JSON']
        #OFF_cdp.path=self.inputs['subpaths']['products']

        #OFFdict = OrderedDict()
        #OFFdict['PRE'] = self.dd.products['MED_PRE'].copy()
        #OFFdict['IMG'] = self.dd.products['MED_IMG'].copy()
        #OFFdict['OVE'] = self.dd.products['MED_OVE'].copy()


        #OFF_cdp.ingest_inputs(data=OFFdict,
        #    header=OFF_hdr,
        #    meta=dict(units='ADU'))
        #OFF_cdp.savehardcopy()



        # Matrix of RONs (over-scan)

        ron_tb_cdp = cdp.Tables_CDP()
        roncdpdict = dict(
            caption='RON in the over-scan region.',
            valformat='%.2f')

        def _getRON_OVEval(self, Ckey, Q):
            return self.dd.products['STD_OVE'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getRON_OVEval,
                                    cdp=ron_tb_cdp,
                                    cdpdict=roncdpdict)

        # Save RONs to a json

        #OFF_hdr = OrderedDict()
        #OFF_hdr['title'] = 'OFFSETS'
        #OFF_hdr.update(CDP_header)

        #OFF_cdp = self.CDP_lib['OFF_JSON']
        #OFF_cdp.path=self.inputs['subpaths']['products']

        #OFFdict = OrderedDict()
        #OFFdict['PRE'] = self.dd.products['MED_PRE'].copy()
        #OFFdict['IMG'] = self.dd.products['MED_IMG'].copy()
        #OFFdict['OVE'] = self.dd.products['MED_OVE'].copy()


        #OFF_cdp.ingest_inputs(data=OFFdict,
        #    header=OFF_hdr,
        #    meta=dict(units='ADU'))
        #OFF_cdp.savehardcopy()


    def _get_XYdict_VPROFS(self):
        """ """
        
        data = dict(pre=self.dd.products['profiles_V_1D_pre'],
                    img=self.dd.products['profiles_V_1D_img'],
                    ove=self.dd.products['profiles_V_1D_ove'])

        regions = ['pre','img','ove']

        XYdict_VPROFS = dict(x=dict(),
                            y=dict(),
                            labelkeys=regions)
        for reg in regions:
            XYdict_VPROFS['x'][reg] = []
            XYdict_VPROFS['y'][reg] = []


        def assigner(VPROFS, data, Ckey):

            for Q in self.Quads:

                cenY = []

                for reg in regions:
                
                    _x = data[reg][Ckey][Q]['x'].copy()
                    _y = data[reg][Ckey][Q]['y'].copy()
                    VPROFS['x'][reg].append(_x)
                    VPROFS['y'][reg].append(_y)
                    cenY.append(np.median(_y))
                    #VPROFS['x'][reg] = [_x]
                    #VPROFS['y'][reg] =[_y]

                cenY = np.mean(cenY)

                for reg in regions: # recentering profiles
                    VPROFS['y'][reg][-1] -= cenY

            return VPROFS

        XYdict_VPROFS = self.iter_overCCDs(data, assigner, RetDict=XYdict_VPROFS)

        for reg in regions:
            XYdict_VPROFS['x'][reg] = np.array(XYdict_VPROFS['x'][reg]).T.tolist()
            XYdict_VPROFS['y'][reg] = np.array(XYdict_VPROFS['y'][reg]).T.tolist()


        return XYdict_VPROFS

    def _get_XYdict_HPROFS(self):
        """ """

        regions = ['pre','ove']

        XYdict_HPROFS = dict(x=dict(),
                            y=dict())
        for reg in regions:
            XYdict_HPROFS['x'][reg] = []
            XYdict_HPROFS['y'][reg] = []

        xgap = 3

        def assigner(HPROFS, data, Ckey):

            

            for Q in self.Quads:

                #_x = data[Ckey][Q]['x'].copy()
                _y = data[Ckey][Q]['y'].copy()

                _y -= np.nanmean(_y[-29+5:-5+1]) # subtract ove-scan offset
                
                _ypre = _y[0:51+20]
                _yove = _y[-(29+20):]

                HPROFS['x']['pre'].append(np.arange(len(_ypre)))
                HPROFS['y']['pre'].append(_ypre)

                HPROFS['x']['ove'].append(np.arange(len(_yove))+len(_ypre)+xgap)
                HPROFS['y']['ove'].append(_yove)
                

            return HPROFS

        XYdict_HPROFS = self.iter_overCCDs(self.dd.products['profiles_H_1D'], 
                        assigner, RetDict=XYdict_HPROFS)

        for reg in regions:
            XYdict_HPROFS['x'][reg] = np.array(XYdict_HPROFS['x'][reg]).T.tolist()
            XYdict_HPROFS['y'][reg] = np.array(XYdict_HPROFS['y'][reg]).T.tolist()


        XYdict_HPROFS['x']['pre2img'] = np.array([1,1])*51
        XYdict_HPROFS['y']['pre2img'] = np.array([-1.E3,1.E3])

        XYdict_HPROFS['x']['img2ove'] = np.array([1,1]) * (51+20+20+xgap)
        XYdict_HPROFS['y']['img2ove'] = np.array([-1.E3,1.E3])

        XYdict_HPROFS['labelkeys'] = ['pre','pre2img','img2ove','ove']

        return XYdict_HPROFS

    def _get_Histos_OFFSETS(self):
        """ """

        data = dict(pre=self.dd.products['MED_PRE'],
                    img=self.dd.products['MED_IMG'],
                    ove=self.dd.products['MED_OVE'])

        regions = ['pre','img','ove']

        HistosDict = dict(x=dict(),
                            y=dict(),
                            labelkeys=regions)

        for reg in regions:
            HistosDict['y'][reg] = []

        def assigner(Histos, data, Ckey):
            for Q in self.Quads:
                for reg in regions:
                    Histos['y'][reg].append(data[reg][Ckey][Q])
            return Histos

        HistosDict = self.iter_overCCDs(data, assigner, RetDict=HistosDict)

        minval = 1e6
        maxval = -1e6

        for reg in regions:

            regmin = int(np.nanmin(HistosDict['y'][reg]))
            regmax =  int(np.nanmax(HistosDict['y'][reg]))
            minval = min((minval,regmin))
            maxval = max((maxval,regmax))

        bin_edges = np.linspace(minval,maxval,min((20,int(maxval-minval))))
        
        for reg in regions:
            HistosDict['x'][reg] = bin_edges

        return HistosDict



    def meta_analysis(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)


        # Offsets in pre-/img-/ove-: histograms

        OFFSETS_HISTOS = self._get_Histos_OFFSETS()

        self.figdict['OFFSETS_HISTO'][1]['data'] = OFFSETS_HISTOS
        self.figdict['OFFSETS_HISTO'][1]['meta']['plotter'] = self.metacal.plot_HISTO


        # map difference in offset with reference

        def _get_GenValue_MAP(inData, Ckey, Q):
            return inData[Ckey][Q]

        OffsetsMap = self.get_FPAMAP(self.dd.products['MED_OVE'],
                                     extractor=_get_GenValue_MAP)

        RefOffsetsMap = self.dd.products['REF_OFFs'].copy()

        def _get_Diff_MAP(inData, Ckey, Q):
            """ """
            return inData[0][Ckey][Q]-inData[1][Ckey][Q]

        DiffOffsetsMap = self.get_FPAMAP((OffsetsMap,RefOffsetsMap),
                                        extractor=_get_Diff_MAP)
        
        self.figdict['DIFFOFFSETSMAP'][1]['data'] = DiffOffsetsMap
        self.figdict['DIFFOFFSETSMAP'][1]['meta']['plotter'] = self.metacal.plot_SimpleMAP

        DiffOffsetsList = []
        for Ckey in list(DiffOffsetsMap.keys()):
            for Q in self.Quads:
                DiffOffsetsList.append(DiffOffsetsMap[Ckey][Q])

        avgDiffOffset = np.mean(DiffOffsetsList)

        if self.report is not None:
            self.report.add_Text('Average Differerence in Offset with reference: %.2f ADU' % avgDiffOffset)

        # map RATIO of RON with reference


        RonsMap = self.get_FPAMAP(self.dd.products['STD_OVE'],
                                     extractor=_get_GenValue_MAP)

        RefRonsMap = self.dd.products['REF_RONs'].copy()

        def _get_Ratio_MAP(inData, Ckey, Q):
            """ """
            return inData[0][Ckey][Q]/inData[1][Ckey][Q]

        RatioRonsMap = self.get_FPAMAP((RonsMap,RefRonsMap),
                                        extractor=_get_Ratio_MAP)
        
        self.figdict['RATIORONSMAP'][1]['data'] = RatioRonsMap
        self.figdict['RATIORONSMAP'][1]['meta']['plotter'] = self.metacal.plot_SimpleMAP

        RatioRonsList = []
        for Ckey in list(RatioRonsMap.keys()):
            for Q in self.Quads:
                RatioRonsList.append(RatioRonsMap[Ckey][Q])

        avgRatioRon = np.mean(RatioRonsList)

        if self.report is not None:
            self.report.add_Text('Average Ratio of RON with reference: %.2f (adim.)' % avgRatioRon)


        # MATRIX of RON (electrons)

        ronele_tb_cdp = cdp.Tables_CDP()
        ronelecdpdict = dict(
            caption='RON [electrons] in the over-scan region.',
            valformat='%.2f')

        def _getRONele_OVEval(self, Ckey, Q):
            BLOCK, CCDk, _, _ = self.fpa.FPA_MAP[Ckey]
            gain = self.incdps['GAIN'][BLOCK][CCDk][Q][0]
            return self.dd.products['STD_OVE'][Ckey][Q] * gain


        TBronele, cdpronele = self.add_StandardQuadsTable(extractor=_getRONele_OVEval,
                                    cdp=ronele_tb_cdp,
                                    cdpdict=ronelecdpdict)
        RonEleList = []
        for Q in self.Quads:
            RonEleList += TBronele[Q].tolist()
        RonEleList = np.array(RonEleList)
        if self.report is not None:
            minronele = np.nanmin(RonEleList)
            maxronele = np.nanmax(RonEleList)
            Naboveronreq = len(np.where(RonEleList > 4.5)[0])
            
            self.report.add_Text('Min Ron [electrons]: %.2f e-' % minronele)
            self.report.add_Text('Max Ron [electrons]: %.2f e-' % maxronele)
            self.report.add_Text('N(Ron$>$4.5 e- rms): %i' % Naboveronreq)


        XYdict_HPROFS = self._get_XYdict_HPROFS()

        self.figdict['HPROFS'][1]['data'] = XYdict_HPROFS
        self.figdict['HPROFS'][1]['meta']['plotter'] = self.metacal.plot_XY        

        # show average profiles along rows (single plot)

        # show average profiles along columns (single plot): pre-/img-/ove-
        #       Is there stray-light?

        XYdict_VPROFS = self._get_XYdict_VPROFS()

        self.figdict['VPROFS'][1]['data'] = XYdict_VPROFS
        self.figdict['VPROFS'][1]['meta']['plotter'] = self.metacal.plot_XY


        

        # PLOTTING
        

        if self.report is not None:
            self.addFigures_ST(figkeys=['OFFSETS_HISTO','DIFFOFFSETSMAP', 
                                        'RATIORONSMAP','VPROFS','HPROFS'],
                               dobuilddata=False)


    def appendix(self):
        """Adds Appendices to Report."""


        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level = 0)


        # TABLE: reference values of OFFSETS
        
        def _getRefOffs(self, Ckey, Q):
            return self.dd.products['REF_OFFs'][Ckey][Q]
        
        cdpdictoff = dict(
            caption = 'Reference OFFSETs [ADU] in Over-scan. (from GRCALCAMP).',
            valformat = '%.1f')
        
        self.add_StandardQuadsTable(extractor=_getRefOffs,
                                    cdp=None,
                                    cdpdict=cdpdictoff)


        # TABLE: reference values of RONS
        
        
        def _getRefRons(self, Ckey, Q):
            return self.dd.products['REF_RONs'][Ckey][Q]
        
        cdpdictron = dict(
            caption = 'Reference RONs (rms, [ADU]) (from GRCALCAMP).',
            valformat = '%.2f')
        
        self.add_StandardQuadsTable(extractor=_getRefRons,
                                    cdp=None,
                                    cdpdict=cdpdictron)
        
