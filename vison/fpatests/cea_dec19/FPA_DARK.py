# -*- coding: utf-8 -*-
"""

TEST: FPA-DARK (COLD)

Created on Thu Nov 14 16:27:00 2019

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
from vison.fpatests.cea_dec19 import DK_aux
from vison.inject import lib as ilib
# END IMPORT


class DARK_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptime', ([float], 'Exposure Time')),
    ])))

class DARK(fpatask.FpaTask):

    inputsclass = DARK_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.subtasks = [('check', self.check_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis),
                         ('appendix', self.appendix)]

        super(DARK, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                                       cleanafter=cleanafter)

        self.name = 'DARK'
        self.type = 'Simple'


        self.figdict = DK_aux.get_DKfigs()
        self.CDP_lib = DK_aux.get_CDP_lib()

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

        # Load Reference Offsets
        
        refoffkey, offreg = ('offsets_fwd_cold','ove')

        refOFF_incdp = cdpmod.Json_CDP()
        refOFF_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refoffkey])
        

        def _get_ref_OFFs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][offreg]
        
        
        RefOFFsMap = self.get_FPAMAP(refOFF_incdp.data.copy(),
                                        extractor=_get_ref_OFFs_MAP)

        self.dd.products['REF_OFFs'] = RefOFFsMap.copy()


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

                self.dd.products['profiles_H_1D'][CCDID][Q] = hor1Dprof


                for reg in ['pre','img','ove']:

                    # produce average profile along cols
                    _ver1Dprof = kccdobj.get_1Dprofile(
                                Q=Q, orient='ver', area=reg, stacker='mean', 
                                vstart=vstart, vend=vend)
            
                    self.dd.products['profiles_V_1D_%s' % reg][CCDID][Q] = _ver1Dprof


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

                # DARK Measurement: PENDING
                

    def basic_analysis(self):
        """
        To-Do:
            Extract STDs and MEDs in pre/img/ove

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

        prodkeys = []
        for reg in ['pre','img','ove']:
            prodkeys.append('profiles_V_1D_%s' % reg)

        prodkeys += ['profiles_H_1D', 
            'STD_PRE', 'STD_IMG','STD_OVE',
            'MED_PRE', 'MED_IMG', 'MED_OVE',
            'DARK', 'EDARK']

        for prodkey in prodkeys:

            self.dd.products[prodkey] = dict()

            for jY in range(1, self.NSLICES_FPA + 1):
                for iX in range(1, self.NCOLS_FPA + 1):
                    CCDID = 'C_%i%i' % (jY, iX)
                    self.dd.products[prodkey][CCDID] = dict()

        Bkwargs = dict(vstart=vstart,
                       vend=vend,
                       debug=False)

        self.iterate_over_CCDs_inLE1(LE1, DARK._basic_onLE1, **Bkwargs)

        if self.report is not None:
            for prodkey in prodkeys:
                nprodkey = prodkey.replace('_','\_')
                self.report.add_Text('product: %s, all extracted!' % nprodkey)



        # Subtract OFFSET for display

        def offset_subtractor(ccdobj, **kwargs):
            """ """
            for Q in ccdobj.Quads:
                ccdobj.sub_offset(Q, **kwargs)
            return ccdobj
        
        suboffkwargs = dict(method='row', scan='ove', trimscan=[5, 5],
                   ignore_pover=True, extension=-1)

        LE1.apply_function_to_ccds(offset_subtractor, **suboffkwargs)

        # DISPLAY IMAGE

        SubOffImgDict = self.get_ImgDictfromLE1(LE1, doequalise=True)

        self.figdict['DK_img'][1]['data'] = SubOffImgDict
        self.figdict['DK_img'][1]['meta']['plotter'] = self.metacal.plot_ImgFPA

        if self.report is not None:
            self.addFigures_ST(figkeys=['DK_img'], dobuilddata=False)

    def meta_analysis(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)




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
        