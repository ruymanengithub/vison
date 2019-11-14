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

from vison.support import utils
from vison.image import bits
from vison.other import MOT_FFaux
from vison.fpatests import fpatask
from vison.datamodel import inputs
from vison.datamodel import cdp as cdpmod
from vison.fpatests.cea_dec19 import BIAS_aux

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

        self.figdict = BIAS_aux.get_Bfigs(self.inputs['readmode'],self.inputs['temperature'])
        self.CDP_lib = dict()

        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')


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

            	# produce average profile along cols
            	ver1Dprof = kccdobj.get_1Dprofile(
                                Q=Q, orient='ver', area='all', stacker='mean', vstart=vstart, vend=vend)

            	
            	self.dd.products['profiles_V_1D'][CCDID][Q] = ver1Dprof


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

        prodkeys = ['profiles_V_1D', 'profiles_H_1D', 
        			'STD_PRE', 'STD_IMG','STD_OVE',
        			'MED_PRE', 'MED_IMG', 'MED_OVE']

        for prodkey in prodkeys:

            self.dd.products[prodkey] = dict()

            for jY in range(1, self.NSLICES_FPA + 1):
                for iX in range(1, self.NCOLS_FPA + 1):
                    CCDID = 'C_%i%i' % (jY, iX)
                    self.dd.products[prodkey][CCDID] = dict()

        Bkwargs = dict(vstart=vstart,
                       vend=vend,
                       debug=False)

        self.iterate_over_CCDs(LE1, FPA_BIAS._basic_onLE1, **Bkwargs)

        if self.report is not None:
            for prodkey in prodkeys:
                nprodkey = prodkey.replace('_','\_')
                self.report.add_Text('product: %s, all extracted!' % nprodkey)

        # # Matrices: HERvalues, RAMPslopes

        # # HERvalues matrix

        # her_tb_cdp = self.CDP_lib['HER_TB']
        # hercdpdict = dict(
        #     TBkey='HER_TB',
        #     meta=dict(),
        #     CDP_header=CDP_header,
        #     header_title='FWD WARM: HER TABLE',
        #     CDP_KEY='HER_TB_CDP',
        #     caption='HER value at first pixel.',
        #     valformat='%.2e'
        # )

        # def _getHERval(self, Ckey, Q):
        #     return self.dd.products['HERval'][Ckey][Q]

        # self.add_StandardQuadsTable(extractor=_getHERval,
        #                             cdp=her_tb_cdp,
        #                             cdpdict=hercdpdict)

        # # RAMP Slopes

        # rslope_tb_cdp = self.CDP_lib['RSLOPE_TB']
        # rslopecdpdict = dict(
        #     TBkey='RSLOPE_TB',
        #     meta=dict(),
        #     CDP_header=CDP_header,
        #     header_title='FWD WARM: RAMP SLOPES TABLE',
        #     CDP_KEY='RSLOPE_TB_CDP',
        #     caption='RAMP Slope in ADU/pixel row',
        #     valformat='%.1f'
        # )

        # def _getRSlope(self, Ckey, Q):
        #     return self.dd.products['RAMPfits'][Ckey][Q][0]

        # self.add_StandardQuadsTable(extractor=_getRSlope,
        #                             cdp=rslope_tb_cdp,
        #                             cdpdict=rslopecdpdict)



    def meta_analysis(self):
        """ """
        stop()


    def appendix(self):
        """Adds Appendices to Report."""



        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level=0)

        return

        # TABLE: reference values of OFFSETS
        
        def _getRefVal1(self, Ckey, Q):
            return self.dd.products['REFsomething1'][Ckey][Q]
        
        cdpdict1 = dict(
            caption = 'Reference OFFSETS (GRCALCAMP).',
            valformat = '%.1f')
        
        self.add_StandardQuadsTable(extractor=_getRefVal,
                                    cdp=None,
                                    cdpdict=cdpdict1)



        # TABLE: reference values of RONS
        
        
        def _getRefVal2(self, Ckey, Q):
            return self.dd.products['REFsomething1'][Ckey][Q]
        
        cdpdict2 = dict(
            caption = 'Reference OFFSETS (GRCALCAMP).',
            valformat = '%.1f')
        
        self.add_StandardQuadsTable(extractor=_getRefVal2,
                                    cdp=None,
                                    cdpdict=cdpdict2)
        
