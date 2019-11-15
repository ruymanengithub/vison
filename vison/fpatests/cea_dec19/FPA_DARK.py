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

            

            if not debug:

                pass





    def basic_analysis(self):
        """
        To-Do:
            PENDING
        """

        if self.report is not None:
            self.report.add_Section(keyword='extract', Title='Extraction', level=0)

 


    def meta_analysis(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)

        
        refoffkey, offreg = ('offsets_fwd_cold','ove')

        refOFF_incdp = cdpmod.Json_CDP()
        refOFF_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refoffkey])
        

        def _get_ref_OFFs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][offreg]
        
        
        RefOFFsMap = self.get_FPAMAP(refOFF_incdp.data.copy(),
                                        extractor=_get_ref_OFFs_MAP)

        self.dd.products['REF_OFFs'] = RefOFFsMap.copy()



    def appendix(self):
        """Adds Appendices to Report."""

        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level=0)

        # TABLE: reference values of OFFSETS
        
        def _getRefOffs(self, Ckey, Q):
            return self.dd.products['REF_OFFs'][Ckey][Q]
        
        cdpdictoff = dict(
            caption = 'Reference OFFSETS (GRCALCAMP).',
            valformat = '%.1f')
        
        self.add_StandardQuadsTable(extractor=_getRefOffs,
                                    cdp=None,
                                    cdpdict=cdpdictoff)
        