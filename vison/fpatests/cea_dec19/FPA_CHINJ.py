# -*- coding: utf-8 -*-
"""

TEST: FPA-CHINJ (COLD)

Created on Thu Nov 14 13:05:00 2019

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
from vison.fpatests.cea_dec19 import CI_aux
from vison.inject import lib as ilib
# END IMPORT


class CHINJ_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('non', ([int], 'Number of lines ON')),
        ('noff', ([int], 'Number of lines OFF'))
    ])))

class CHINJ(fpatask.FpaTask):

    inputsclass = CHINJ_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.subtasks = [('check', self.check_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis),
                         ('appendix', self.appendix)]

        super(CHINJ, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                                       cleanafter=cleanafter)

        self.name = 'CHINJ'
        self.type = 'Simple'


        self.figdict = CI_aux.get_CIfigs()
        self.CDP_lib = CI_aux.get_CDP_lib()

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

        # REFERENCE OFFSETS

        refoffkey, offreg = ('offsets_fwd_cold','ove')

        refOFF_incdp = cdpmod.Json_CDP()
        refOFF_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refoffkey])

        def _get_ref_OFFs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][offreg]

        RefOFFsMap = self.get_FPAMAP(refOFF_incdp.data.copy(),
                                        extractor=_get_ref_OFFs_MAP)

        self.dd.products['REF_OFFs'] = RefOFFsMap.copy()


        # Reference Injection Values

        refINJ_incdp = cdpmod.Json_CDP()
        refINJ_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references']['injection_levels'])

        def _get_ref_INJs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q]

        RefINJsMap = self.get_FPAMAP(refINJ_incdp.data.copy(),
                                        extractor=_get_ref_INJs_MAP)

        self.dd.products['REF_INJs'] = RefINJsMap.copy()


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

        for Q in self.Quads:

            # Extract Charge Injection

            non = self.inputs['non']
            noff = self.inputs['noff']
            nrep = (vend - vstart) // (non + noff) + 1

            pattern = (non, noff, nrep)

            if not debug:

                ext_res = ilib.extract_injection_lines(
                                    kccdobj, Q, pattern, VSTART=vstart, VEND=vend, 
                                    suboffmean=True)

                # profile along columns

                yalcols = ext_res['avprof_alcol'].copy()
                xalcols = np.arange(len(yalcols), dtype='float32')

                self.dd.products['profiles1D_V'][CCDID][Q] = \
                        dict(x=xalcols,y=yalcols)

                # profile along rows

                yalrows = ext_res['avprof_alrow'].copy()
                xalrows = np.arange(len(yalrows), dtype='float32')

                self.dd.products['profiles1D_H'][CCDID][Q] = \
                        dict(x=xalrows,y=yalrows)                


                # save INJECTION value

                self.dd.products['INJECTION'][CCDID][Q] = \
                        ext_res['stats_injection']['p50']

                # Extract Overscan Stats


                reg = 'ove'
                ovestats = kccdobj.get_stats(
                                    Q,sector=reg,statkeys=['median','std'],
                                    trimscan=trimscans[reg],
                                    ignore_pover=True,extension=-1,
                                    VSTART=vstart,VEND=vend)

                self.dd.products['OFF_OVE'][CCDID][Q] = ovestats[0]
                self.dd.products['STD_OVE'][CCDID][Q] = ovestats[1]



    def basic_analysis(self):
        """
        To-Do:
            Extract average profiles along columns (image area)
            Extract Injection Values
            Extract Over-scan stats (offset + ron)
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

        CIImgDict = self.get_ImgDictfromLE1(LE1, doequalise=True)

        self.figdict['CI_img'][1]['data'] = CIImgDict
        self.figdict['CI_img'][1]['meta']['plotter'] = self.metacal.plot_ImgFPA

        if self.report is not None:
            self.addFigures_ST(figkeys=['CI_img'], dobuilddata=False)

        prodkeys = ['profiles1D_H', 'profiles1D_V',
                    'INJECTION','STD_OVE','OFF_OVE']

        for prodkey in prodkeys:

            self.dd.products[prodkey] = dict()

            for jY in range(1, self.NSLICES_FPA + 1):
                for iX in range(1, self.NCOLS_FPA + 1):
                    CCDID = 'C_%i%i' % (jY, iX)
                    self.dd.products[prodkey][CCDID] = dict()


        CIkwargs = dict(vstart=vstart,
                       vend=vend,
                       debug=False)


        self.iterate_over_CCDs_inLE1(LE1, CHINJ._basic_onLE1, **CIkwargs)

        if self.report is not None:
            for prodkey in prodkeys:
                nprodkey = prodkey.replace('_','\_')
                self.report.add_Text('product: %s, all extracted!' % nprodkey)


        # Matrix of Injection Levels

        inj_tb_cdp = self.CDP_lib['INJ_TB']
        injcdpdict = dict(
            TBkey='INJ_TB',
            meta=dict(),
            CDP_header=CDP_header,
            header_title='CHINJ: INJECTION TABLE',
            CDP_KEY='INJ_TB_CDP',
            caption='Median Injection values in ADU.',
            valformat='%.1f'
        )

        def _getINJval(self, Ckey, Q):
            return self.dd.products['INJECTION'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getINJval,
                                    cdp=inj_tb_cdp,
                                    cdpdict=injcdpdict)

        # FPA MAP with extracted Injection profiles


        def _assignProf1D(InjProfsMapDict, profdict, Ckey):
            """ """
            Cprofs = profdict[Ckey]

            res = dict(x=dict(), y=dict())

            for Q in self.Quads:
                res['x'][Q] = Cprofs[Q]['x']
                res['y'][Q] = Cprofs[Q]['y']

            InjProfsMapDict[Ckey] = res

            return InjProfsMapDict

        InjProfsMapDict = self.iter_overCCDs(self.dd.products['profiles1D_V'], _assignProf1D)
        InjProfsMapDict['labelkeys'] = self.Quads

        self.figdict['INJ_PROFS_MAP'][1]['data'] = InjProfsMapDict
        self.figdict['INJ_PROFS_MAP'][1]['meta']['plotter'] = self.metacal.plot_XYMAP


        if self.report is not None:
            self.addFigures_ST(figkeys=['INJ_PROFS_MAP'],
                               dobuilddata=False)

    def _get_XYdict_PROFS(self, kind):
        """ """

        Quads = self.Quads

        XYdict_PROFS = dict(x=dict(),
                            y=dict())
        for Q in Quads:
            XYdict_PROFS['x'][Q] = []
            XYdict_PROFS['y'][Q] = []


        def assigner(PROFS, data, Ckey):


            for Q in self.Quads:

                _x = data[Ckey][Q]['x'].copy()
                _y = data[Ckey][Q]['y'].copy()


                PROFS['x'][Q].append(_x)
                PROFS['y'][Q].append(_y)


            return PROFS

        if kind == 'HOR':
            data = self.dd.products['profiles1D_H']
        elif kind == 'VER':
            data = self.dd.products['profiles1D_V']

        XYdict_PROFS = self.iter_overCCDs(data, 
                        assigner, RetDict=XYdict_PROFS)

        for Q in Quads:
            XYdict_PROFS['x'][Q] = np.array(XYdict_PROFS['x'][Q]).T.tolist()
            XYdict_PROFS['y'][Q] = np.array(XYdict_PROFS['y'][Q]).T.tolist()

        XYdict_PROFS['labelkeys'] = Quads

        return XYdict_PROFS



    def meta_analysis(self):
        """ """

        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)



        # Display map of injection levels minus reference

        INJsMap = self.dd.products['INJECTION'].copy()

        def _get_Diff_Inj_MAP(inData, Ckey, Q):
            return inData[0][Ckey][Q]-inData[1][Ckey][Q]

        RefINJsMap = self.dd.products['REF_INJs']

        DiffInjMap = self.get_FPAMAP((INJsMap,RefINJsMap),
                                        extractor=_get_Diff_Inj_MAP)

        self.figdict['DIFF_INJMAP'][1]['data'] = DiffInjMap
        self.figdict['DIFF_INJMAP'][1]['meta']['plotter'] = self.metacal.plot_SimpleMAP

        # Show average horizontal injection profiles

        XYdict_HPROFS = self._get_XYdict_PROFS(kind='HOR')

        self.figdict['INJ_PROFS_HOR'][1]['data'] = XYdict_HPROFS
        self.figdict['INJ_PROFS_HOR'][1]['meta']['plotter'] = self.metacal.plot_XY        

        # Show average vertical injection profiles

        XYdict_VPROFS = self._get_XYdict_PROFS(kind='VER')

        self.figdict['INJ_PROFS_VER'][1]['data'] = XYdict_VPROFS
        self.figdict['INJ_PROFS_VER'][1]['meta']['plotter'] = self.metacal.plot_XY


        if self.report is not None:
            self.addFigures_ST(figkeys=['DIFF_INJMAP','INJ_PROFS_VER',
                            'INJ_PROFS_HOR'],
                               dobuilddata=False) 



    def appendix(self):
        """Adds Appendices to Report."""

        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level=0)


        # TABLE: reference values of OFFSETS

        def _getRefOffs(self, Ckey, Q):
            return self.dd.products['REF_OFFs'][Ckey][Q]

        cdpdictoff = dict(
            caption = 'Reference OFFSETS [ADU] (GRCALCAMP).',
            valformat = '%.1f')

        self.add_StandardQuadsTable(extractor=_getRefOffs,
                                    cdp=None,
                                    cdpdict=cdpdictoff)


        # TABLE: reference values of INJECTION

        def _getRefINJs(self, Ckey, Q):
            return self.dd.products['REF_INJs'][Ckey][Q]

        cdpdictinj = dict(
            caption = 'Reference charge injection levels [ADU] (GRCALCAMP).',
            valformat = '%.1f')

        self.add_StandardQuadsTable(extractor=_getRefINJs,
                                    cdp=None,
                                    cdpdict=cdpdictinj)
