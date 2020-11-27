#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Common Use Task for Trap-Pumping Analysis.

Created on Tue Jan  2 17:44:04 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import copy
import os
import pandas as pd
from collections import OrderedDict

from vison.support import utils
from vison.pipe.task import Task
from vison.datamodel import core, ccd
from vison.pipe import lib as pilib
#from vison.pipe.task import Task
from vison.inject.InjTask import InjTask
from vison.support.files import cPickleRead
from vison.pump import tptools
# END IMPORT


class PumpTask(InjTask):

    def __init__(self, *args, **kwargs):
        super(PumpTask, self).__init__(*args, **kwargs)

    def check_data(self, **kwargs):
        """ """
        test = self.inputs['test']
        if test in ['TP01', 'TP11']:
            kwargs = dict(pattern=(2066, 0, 1),
                          figkeys=['TP01checks_offsets', 'TP01checks_deltaoff',
                                   'TP01checks_stds',
                                   'TP01checks_injlevel', 'TP01checks_injstd'])
        elif test in ['TP02', 'TP21']:
            kwargs = dict(pattern=(2066, 0, 1),
                          figkeys=['TP02checks_offsets', 'TP02checks_deltaoff',
                                   'TP02checks_stds',
                                   'TP02checks_injlevel', 'TP02checks_injstd'])

        InjTask.check_data(self, **kwargs)

    def check_metrics_ST(self, **kwargs):
        """

        """
        super(PumpTask, self).check_metrics_ST(**kwargs)

    def charact_injection(self):
        """Characterises Charge Injection."""

        testname = self.inputs['test']

        if testname in ['TP01', 'TP11']:
            testtype = 'v'
        elif testname in ['TP02', 'TP21']:
            testtype = 's'

        if self.report is not None:
            self.report.add_Section(
                keyword='charact', Title='%s: Charge Injection Characterisation' % testname,
                level=0)

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
        chinj_dd = OrderedDict()
        chinjnoise_dd = OrderedDict()

        if not self.drill:

            for i, id_dly in enumerate(id_dlys):

                for jCCD, CCDk in enumerate(CCDs):

                    chinjnoise_dd[CCDk] = OrderedDict()
                    chinj_dd[CCDk] = OrderedDict()

                    ixsel = np.where((self.dd.mx['id_dly'][:] == id_dly) &
                                     (self.dd.mx['%s_tpump' % testtype][:] == 0) &
                                     (self.dd.mx['CCD'][:] == CCDk))

                    iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ixsel][0]
                    iccdobj = cPickleRead(os.path.join(ccdpicklespath, iccdobj_f))

                    char_res_dict[CCDk] = tptools.charact_injection(iccdobj)

                    vstart = int(iccdobj.extensions[-1].header['vstart'])
                    vend = min(int(iccdobj.extensions[-1].header['vend']), iccdobj.NrowsCCD)
                    Nrows = vend - vstart

                    if (i == 0) and jCCD == 0:
                        pdeg = char_res_dict[CCDk]['pdeg']
                    char_res_dict[CCDk].pop('pdeg')

                    for Q in Quads:
                        chinjnoise_dd[CCDk][Q] = np.nanmedian(char_res_dict[CCDk][Q]['injnoise'])

                        pco = char_res_dict[CCDk][Q]['polycoeffs']
                        X = np.arange(Nrows)

                        injmap = np.outer(pco[:, 0], X**2.) + np.outer(pco[:, 1], X) +\
                            np.outer(pco[:, 2], np.ones((Nrows)))

                        chinj_dd[CCDk][Q] = np.nanmedian(injmap)

        # Saving the charge injection characterisation results

        self.dd.products['chinj_mx'] = copy.deepcopy(chinj_dd)
        self.dd.products['chinjnoise_mx'] = copy.deepcopy(chinjnoise_dd)

        chchar_cdp = self.CDP_lib['CHINJCHARACT']

        chchar_cdp.header = CDP_header.copy()
        chchar_cdp.meta = dict(pdeg=pdeg)
        chchar_cdp.path = productspath
        chchar_cdp.data = char_res_dict.copy()

        self.save_CDP(chchar_cdp)
        self.pack_CDP_to_dd(chchar_cdp, 'CHINJCHARACT')

        # REPORTS

        # Table with injection noise per CCD/Quadrant

        chinj_cdp = self.CDP_lib['CHINJ']
        chinj_cdp.path = self.inputs['subpaths']['products']
        chinj_dddf = OrderedDict(
            CHINJ=pd.DataFrame.from_dict(chinj_dd),
            CHINJNOISE=pd.DataFrame.from_dict(chinjnoise_dd))
        chinj_cdp.ingest_inputs(data=chinj_dddf.copy(),
                                meta=dict(),
                                header=CDP_header.copy())

        chinj_cdp.init_wb_and_fillAll(header_title='%s: CHARGE INJECTION' % self.inputs['test'])
        self.save_CDP(chinj_cdp)
        self.pack_CDP_to_dd(chinj_cdp, 'CHINJ_CDP')

        # Table with injection levels per CCD/Quadrant

        if self.report is not None:

            def ff(x): return '%.1f' % x

            formatters = [ff, ff, ff]
            
            CHINJtex = chinj_cdp.get_textable(sheet='CHINJ',
                                              caption='%s: Charge Injection Levels (ADU).' %
                                              self.inputs['test'],
                                              longtable=False,
                                              fitwidth=True,
                                              index=True,
                                              formatters=formatters)
            self.report.add_Text(CHINJtex)

        if self.report is not None:

            def ff(x): return '%.2f' % x

            formatters = [ff, ff, ff]

            CHINJNOISEtex = chinj_cdp.get_textable(
                sheet='CHINJNOISE',
                caption='%s: Charge Injection Noise (ADU, rms).' %
                self.inputs['test'],
                longtable=False,
                fitwidth=True,
                index=True,
                formatters=formatters)
            self.report.add_Text(CHINJNOISEtex)
