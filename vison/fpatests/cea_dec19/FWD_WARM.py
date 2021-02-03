# -*- coding: utf-8 -*-
"""

TEST: FPA-FWD (WARM)

Created on Thu Sep 26 16:11:14 2019

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

from vison.support import utils
from vison.image import bits
from vison.other import MOT_FFaux
from vison.fpatests import fpatask
from vison.datamodel import inputs
from vison.datamodel import cdp as cdpmod
from vison.fpatests.cea_dec19 import FW_aux

# END IMPORT


class FWD_WARM_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()


class FWD_WARM(fpatask.FpaTask):

    inputsclass = FWD_WARM_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.subtasks = [('check', self.check_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis),
                         ('appendix', self.appendix)]

        super(FWD_WARM, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                                       cleanafter=cleanafter)

        self.name = 'FWD_WARM'
        self.type = 'Simple'

        self.figdict = FW_aux.get_FWfigs()
        self.CDP_lib = FW_aux.get_CDP_lib()

        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')


    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(preprocessing=dict(
                                            offsetkwargs=dict(
                                                ignore_pover=True,
                                                trimscan=[25, 5],
                                                method='median',
                                                extension=-1,
                                                scan='pre'
                                            )))

    def load_references(self):
        """ """

        # Reference Ramps

        reframps_incdp = cdpmod.Json_CDP()
        reframps_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references']['warm_ramps'])

        def _get_ref_Rslopes_MAP(inData, Ckey, Q):
            return inData[Ckey][Q]['slope']


        RefRslopesMap = self.get_FPAMAP(reframps_incdp.data.copy(),
                                        extractor=_get_ref_Rslopes_MAP)

        self.dd.products['REFRAMPslopes'] = RefRslopesMap.copy()

        # Reference offsets

        refoffkey, offreg = 'offsets_fwd_cold', 'pre'

        refOFF_incdp = cdpmod.Json_CDP()
        refOFF_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refoffkey])


        def _get_ref_OFFs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][offreg]


        RefOFFsMap = self.get_FPAMAP(refOFF_incdp.data.copy(),
                                        extractor=_get_ref_OFFs_MAP)

        self.dd.products['REF_OFFs'] = RefOFFsMap.copy()


        # Reference RONs

        refronkey, ronreg = 'rons', 'pre'

        refRON_incdp = cdpmod.Json_CDP()
        refRON_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refronkey])


        def _get_ref_RONs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][ronreg]


        RefRONsMap = self.get_FPAMAP(refRON_incdp.data.copy(),
                                        extractor=_get_ref_RONs_MAP)

        self.dd.products['REF_RONs'] = RefRONsMap.copy()




    def _get_RAMPslope(self, vQ):
        """ """

        rawy = vQ.data['y'].copy()
        rownum = np.arange(len(rawy))
        try:
            rowsat = np.where(rawy == 2.**16 - 1)[0][0]
        except IndexError: # no saturation for whatever reason
            rowsat = len(rawy)
        ixsel = np.where((rawy < 5.E4) & (rownum < rowsat))

        try:
            pol = np.polyfit(rownum[ixsel], rawy[ixsel], 1)
            slope, intercept = pol[0], pol[1]
        except:
            slope, intercept = np.nan, np.nan

        return slope, intercept

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

        HERprofs = MOT_FFaux.extract_transcan_profiles(kccdobj,
                                                       [1.E3, 4.E4],
                                                       direction='serial')

        for Q in self.Quads:

            # HERprofile

            #avgQ = kccdobj.get_stats(Q,sector='img',statkeys=['mean'])[0]
            #print('%s%s %.2f' % (CCDID,Q,avgQ))

            # save HERprof

            self.dd.products['HER'][CCDID][Q] = HERprofs[Q].copy()

            # save HER value

            ixjump = HERprofs['ixjump']
            self.dd.products['HERval'][CCDID][Q] = HERprofs[Q]['y'][ixjump]

            if not debug:

                # Stats on pre-scan

                stats = kccdobj.get_stats(
                    Q,sector='pre',
                    statkeys=['median','std'],
                    trimscan=trimscans['pre'],
                    ignore_pover=True,
                    extension=-1,
                    VSTART=vstart,
                    VEND=vend)

                self.dd.products['MED_PRE'][CCDID][Q] = stats[0]
                self.dd.products['STD_PRE'][CCDID][Q] = stats[1]

                # vertical profile: RAMP

                vQ = kccdobj.get_1Dprofile(Q=Q, orient='ver', area='img', stacker='median',
                                           vstart=vstart, vend=vend)

                vQdata = vQ.data.copy()
                vx = vQdata['x'].copy()
                vy = vQdata['y'].copy()
                vx -= vx.min()
                if Q in ['E','F']:
                    vx = kccdobj.NAXIS2//2-vx

                ixsort = np.argsort(vx)
                vx = vx[ixsort]
                vy = vy[ixsort]

                vQdata = dict(x=vx.copy(),y=vy.copy())

                self.dd.products['profiles1D'][CCDID][Q] = vQdata.copy()

                # extract and save slope+intercept of profile for comparison

                rampslopeQ, rampinterQ = self._get_RAMPslope(vQ)

                self.dd.products['RAMPfits'][CCDID][Q] = (rampslopeQ, rampinterQ)

                # RAMP; bit-histograms extraction
                try:
                    ixsat = np.where(vy==2.**16-1)[0][0]
                except IndexError:
                    ixsat = len(vy)


                bitsmean = bits.get_histo_bits(kccdobj, Q, vstart=vstart, vend=ixsat)
                bitsbin = np.arange(0, 16) + 0.5

                # save bit histograms

                self.dd.products['BITS'][CCDID][Q] = dict(bins=bitsbin, H=bitsmean)

    def basic_analysis(self):
        """
        To-Do:
            Extract average profiles along columns (image area)
            Extract HER profiles and metrics
            Extract bits histograms
        """

        if self.report is not None:
            self.report.add_Section(keyword='extract', Title='Extraction', level=0)

        iObs = 0
        vstart = 0
        vend = 2086

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        LE1fits = self.dd.mx['File_name'][iObs]

        fullLE1fits = os.path.join(self.inputs['datapath'], LE1fits)

        LE1 = self.load_LE1(fullLE1fits)

        # DISPLAY IMAGE

        FWImgDict = self.get_ImgDictfromLE1(LE1, doequalise=True)

        self.figdict['FW_img'][1]['data'] = FWImgDict
        self.figdict['FW_img'][1]['meta']['plotter'] = self.metacal.plot_ImgFPA

        if self.report is not None:
            self.addFigures_ST(figkeys=['FW_img'], dobuilddata=False)

        prodkeys = ['profiles1D', 'RAMPfits',
                    'MED_PRE', 'STD_PRE',
                    'HER', 'HERval', 'BITS']

        for prodkey in prodkeys:

            self.dd.products[prodkey] = dict()

            for jY in range(1, self.NSLICES_FPA + 1):
                for iX in range(1, self.NCOLS_FPA + 1):
                    CCDID = 'C_%i%i' % (jY, iX)
                    self.dd.products[prodkey][CCDID] = dict()

        Bkwargs = dict(vstart=vstart,
                       vend=vend,
                       debug=False)

        self.iterate_over_CCDs_inLE1(LE1, FWD_WARM._basic_onLE1, **Bkwargs)

        if self.report is not None:

            for prodkey in prodkeys:
                self.report.add_Text('product: %s, all extracted!' % prodkey)

        # Matrices: HERvalues, RAMPslopes, Offsets, RONs

        # HERvalues matrix

        her_tb_cdp = self.CDP_lib['HER_TB']
        hercdpdict = dict(
            TBkey='HER_TB',
            meta=dict(),
            CDP_header=CDP_header,
            header_title='FWD WARM: HER TABLE',
            CDP_KEY='HER_TB_CDP',
            caption='HER value at first pixel.',
            valformat='%.2e'
        )

        def _getHERval(self, Ckey, Q):
            return self.dd.products['HERval'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getHERval,
                                    cdp=her_tb_cdp,
                                    cdpdict=hercdpdict)

        # RAMP Slopes

        rslope_tb_cdp = self.CDP_lib['RSLOPE_TB']
        rslopecdpdict = dict(
            TBkey='RSLOPE_TB',
            meta=dict(),
            CDP_header=CDP_header,
            header_title='FWD WARM: RAMP SLOPES TABLE',
            CDP_KEY='RSLOPE_TB_CDP',
            caption='RAMP Slope in ADU/pixel row',
            valformat='%.1f'
        )

        def _getRSlope(self, Ckey, Q):
            return self.dd.products['RAMPfits'][Ckey][Q][0]

        self.add_StandardQuadsTable(extractor=_getRSlope,
                                    cdp=rslope_tb_cdp,
                                    cdpdict=rslopecdpdict)

        # Matrix of offsets (pre-scan)

        off_tb_cdp = cdp.Tables_CDP()
        offcdpdict = dict(
            caption='Offsets in the pre-scan region.',
            valformat='%.2f')

        def _getOFF_PREval(self, Ckey, Q):
            return self.dd.products['MED_PRE'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getOFF_PREval,
                                    cdp=off_tb_cdp,
                                    cdpdict=offcdpdict)

        # Matrix of RONs (pre-scan)

        ron_tb_cdp = cdp.Tables_CDP()
        roncdpdict = dict(
            caption='RON in the pre-scan region.',
            valformat='%.2f')

        def _getRON_PREval(self, Ckey, Q):
            return self.dd.products['STD_PRE'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getRON_PREval,
                                    cdp=ron_tb_cdp,
                                    cdpdict=roncdpdict)




    def _get_XYdict_HER(self):
        """ """

        data = self.dd.products['HER']

        def assigner(HERDict, data, Ckey):

            for Q in self.Quads:
                labelkey = '%s_%s' % (Ckey, Q)
                _x = data[Ckey][Q]['x'].copy()
                _x -= _x.min()
                _y = data[Ckey][Q]['y'].copy()

                HERDict['x'][labelkey] = _x.copy()
                HERDict['y'][labelkey] = _y.copy()
                HERDict['labelkeys'].append(labelkey)

            return HERDict

        XYdict_HER = dict(x=dict(),
                          y=dict(),
                          labelkeys=[])

        XYdict_HER = self.iter_overCCDs(data, assigner, RetDict=XYdict_HER)

        return XYdict_HER

    def _get_XYdict_BITS(self):

        data = self.dd.products['BITS']

        def assigner(BITDict, data, Ckey):

            for Q in self.Quads:

                labelkey = '%s_%s' % (Ckey, Q)

                BITDict['x'][labelkey] = data[Ckey][Q]['bins'].copy()
                BITDict['y'][labelkey] = data[Ckey][Q]['H'].copy()
                BITDict['labelkeys'].append(labelkey)

            return BITDict

        XYdict_BIT = dict(x=dict(),
                          y=dict(),
                          labelkeys=[])

        XYdict_BIT = self.iter_overCCDs(data, assigner, RetDict=XYdict_BIT)

        return XYdict_BIT

    def meta_analysis(self):
        """ """

        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)


        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()


        # Display FPA-Map of ramp profiles

        def _assignProf1D(FWRampsDict, profdict, Ckey):
            """ """
            Cprofs = profdict[Ckey]

            res = dict(x=dict(), y=dict())

            for Q in self.Quads:
                res['x'][Q] = Cprofs[Q]['x']
                res['y'][Q] = Cprofs[Q]['y']

            FWRampsDict[Ckey] = res

            return FWRampsDict

        FWRampsDict = self.iter_overCCDs(self.dd.products['profiles1D'], _assignProf1D)
        FWRampsDict['labelkeys'] = self.Quads

        self.figdict['FW_RAMPS'][1]['data'] = FWRampsDict
        self.figdict['FW_RAMPS'][1]['meta']['plotter'] = self.metacal.plot_XYMAP

        # Display FPA-Map of ramp slopes

        def _get_Rslopes_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][0]

        RslopesMap = self.get_FPAMAP(self.dd.products['RAMPfits'],
                                     extractor=_get_Rslopes_MAP)

        self.figdict['SLOPESMAP'][1]['data'] = RslopesMap
        self.figdict['SLOPESMAP'][1]['meta']['plotter'] = self.metacal.plot_SimpleMAP

        # Display FPA-Map of ramp-slope differences with CALCAMP

        def _get_Diff_Rslopes_MAP(inData, Ckey, Q):
            return inData[0][Ckey][Q]-inData[1][Ckey][Q]

        RefRslopesMap = self.dd.products['REFRAMPslopes'].copy()


        DiffRslopesMap = self.get_FPAMAP((RslopesMap,RefRslopesMap),
                                        extractor=_get_Diff_Rslopes_MAP)

        self.figdict['DIFFSLOPESMAP'][1]['data'] = DiffRslopesMap
        self.figdict['DIFFSLOPESMAP'][1]['meta']['plotter'] = self.metacal.plot_SimpleMAP

        # map difference in offset with reference

        def _get_GenValue_MAP(inData, Ckey, Q):
            return inData[Ckey][Q]

        OffsetsMap = self.get_FPAMAP(self.dd.products['MED_PRE'],
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


        RonsMap = self.get_FPAMAP(self.dd.products['STD_PRE'],
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


        # Dispay HER profiles (single plot)

        XYdict_HER = self._get_XYdict_HER()

        self.figdict['HERPROFS'][1]['data'] = XYdict_HER
        self.figdict['HERPROFS'][1]['meta']['plotter'] = self.metacal.plot_XY

        # Display FPA-Map of HER-values

        def _get_HERvals_MAP(inData, Ckey, Q):
            return inData[Ckey][Q]

        HERvalsMap = self.get_FPAMAP(self.dd.products['HERval'],
                                     extractor=_get_HERvals_MAP)

        self.figdict['HERVALSMAP'][1]['data'] = HERvalsMap
        self.figdict['HERVALSMAP'][1]['meta']['plotter'] = self.metacal.plot_SimpleMAP

        # Display FPA-Map of differences with CALCAMP HER-values
        #   On second thoughts, probably better to skip this comparison, it'll 
        #   probably just add "noise"
        # PENDING?

        # Average BIT HISTOGRAM values Matrix

        avbithist_tb_cdp = cdpmod.Tables_CDP()
        avbithistdict = dict(
            caption='Average Value of Bit Histogram for each quadrant in the FPA.',
            valformat='%.2f'
        )

        def _getAvBitHisto(self, Ckey, Q):
            return np.nanmean(self.dd.products['BITS'][Ckey][Q]['H'])

        self.add_StandardQuadsTable(extractor=_getAvBitHisto,
                                    cdp=avbithist_tb_cdp,
                                    cdpdict=avbithistdict)

        # Display all bit-histogram maps together

        XYdict_BITS = self._get_XYdict_BITS()

        self.figdict['BITHISTOS'][1]['data'] = XYdict_BITS
        self.figdict['BITHISTOS'][1]['meta']['plotter'] = self.metacal.plot_XY

        if self.report is not None:
            self.addFigures_ST(figkeys=['FW_RAMPS', 'SLOPESMAP', 'DIFFSLOPESMAP','HERPROFS',
                                        'HERVALSMAP', 'BITHISTOS'],
                               dobuilddata=False)

    def appendix(self):
        """Adds Appendices to Report."""

        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level=0)

        # TABLE: reference values of OFFSETS


        def _getRefOffs(self, Ckey, Q):
            return self.dd.products['REF_OFFs'][Ckey][Q]

        cdpdictoff = dict(
            caption = 'Reference OFFSETs [ADU] in Over-scan. (from GRCALCAMP).',
            valformat = '%.1f')

        self.add_StandardQuadsTable(extractor=_getRefOffs,
                                    cdp=None,
                                    cdpdict=cdpdictoff)

        # TABLE: reference values of SLOPES

        def _getRefRSlope(self, Ckey, Q):
            return self.dd.products['REFRAMPslopes'][Ckey][Q]

        cdpdict = dict(
            caption = 'Reference Slopes [ADU/row] of Dark Current Ramps in FWD-WARM test. Taken from MOT-WARM test'+\
                    ' results at det. chain level (GRCALCAMP).',
            valformat = '%.1f')

        self.add_StandardQuadsTable(extractor=_getRefRSlope,
                                    cdp=None,
                                    cdpdict=cdpdict)

        # TABLE: reference values of HER
