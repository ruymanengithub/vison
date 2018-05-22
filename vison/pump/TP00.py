#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP00 - used to test charge injection functionality using ELVIS.
  NOT intended for performance evaluation.


Created on Tue Jan 13 16:12:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from copy import deepcopy
from collections import OrderedDict

from vison.pipe import lib as pilib
from vison.support import context
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import scriptic as sc
from vison.pipe.task import Task
from PumpTask import PumpTask
from vison.image import performance
from vison.datamodel import inputs
# END IMPORT

isthere = os.path.exists

HKKeys = []

IG1 = 6
IG2 = 5

TP00_commvalues = dict(program='CALCAMP', test='TP00',
                       flushes=7, exptime=0., shuttr=0,
                       e_shuttr=0, vstart=0, vend=20,
                       siflsh=1, siflsh_p=500,
                       IDL=11, IDH=18,
                       IG1_1_T=IG1, IG1_2_T=IG1, IG1_3_T=IG1,
                       IG1_1_B=IG1, IG1_2_B=IG1, IG1_3_B=IG1,
                       IG2_T=IG2, IG2_B=IG2,
                       chinj=1, chinj_on=2066, chinj_of=0,
                       chin_dly=0,
                       id_wid=60,
                       v_tpump=0, s_tpump=0,
                       # v_tp_cnt=1000,
                       dwell_v=0,
                       # dwell_s=0,
                       motr_on=0,
                       comments='')


class TP00_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('Nshuffles_V', ([int], 'Number of Shuffles, Vertical Pumping.')),
        ('toi_tpv', ([list], 'Vector of TOI TP-V values.')),
        ('Nshuffles_S', ([int], 'Number of Shuffles Serial Pumping.')),
        ('dwell_tpsv',
         ([list], 'Vector of dwelling times Serial Pumping.'))
    ])))


class TP00(PumpTask):
    """ """

    inputsclass = TP00_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(TP00, self).__init__(inputs, log, drill, debug)
        self.name = 'TP00'
        self.type = 'Simple'
        self.subtasks = [('check', self.check_data)]
        self.HKKeys = HKKeys
        self.figdict = dict()
        self.inputs['subpaths'] = dict(figs='figs', ccdpickles='ccdpickles')

    def set_inpdefaults(self, **kwargs):
        """ """

        self.inpdefaults = dict(
            Nshuffles_V=500,
            toi_tpv=[200, 1000, 4000],
            Nshuffles_S=100,
            dwell_tpsv=[0, 16, 32],
        )


    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ """

        toi_chinj = 500
        id_delays = np.array([2.5, 1.5]) * toi_chinj

        vpumpmodes = [123, 234, 341, 412]
        Nshuffles_V = self.inputs['Nshuffles_V']
        toi_tpv = self.inputs['toi_tpv']

        Nshuffles_S = self.inputs['Nshuffles_S']
        dwell_tpsv = self.inputs['dwell_tpsv']
        spumpmodes = [23, 31]

        TP00_sdict = dict()

        TP00_commvalues['vstart'] = 0
        TP00_commvalues['vend'] = 20
        TP00_commvalues['toi_ch'] = toi_chinj

        # VERTICAL PUMPING

        colcounter = 1

        for i, id_delay in enumerate(id_delays):

            for m, vpumpmode in enumerate(vpumpmodes):

                for t, toi_tp in enumerate(toi_tpv):

                    colkey = 'col%i' % colcounter

                    TP00_sdict[colkey] = dict(frames=1, toi_tp=toi_tp,
                                              id_dly=id_delay, v_tpump=1, v_tpmod=vpumpmode,
                                              v_tp_cnt=Nshuffles_V,
                                              comments='V%i_%i' % (vpumpmode, toi_tp))

                    colcounter += 1

        # SERIAL PUMPING

        for i, id_delay in enumerate(id_delays):

            for s, spumpmode in enumerate(spumpmodes):

                for d, dwell_tps in enumerate(dwell_tpsv):

                    colkey = 'col%i' % colcounter

                    TP00_sdict[colkey] = dict(frames=1,
                                              id_dly=id_delay, s_tpump=0, s_tpmod=spumpmode,
                                              s_tp_cnt=Nshuffles_S, dwell_s=dwell_tps,
                                              comments='S%i_%i' % (spumpmode, dwell_tps))

                    colcounter += 1

        Ncols = len(TP00_sdict.keys())
        TP00_sdict['Ncols'] = Ncols

        commvalues = deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP00_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        TP00_sdict = sc.update_structdict(TP00_sdict, commvalues, diffvalues)

        return TP00_sdict

    def filterexposures(self, structure, explogf, datapath, OBSID_lims):
        """ """
        wavedkeys = []
        return super(TP00, self).filterexposures(structure, explogf, datapath, OBSID_lims, colorblind=True,
                                                 wavedkeys=wavedkeys)

    def check_data(self):
        """ 

        TP01: Checks quality of ingested data.


        **METACODE**

        ::

            check common HK values are within safe / nominal margins
            check voltages in HK match commanded voltages, within margins

            f.e.ObsID:
                f.e.CCD:
                    f.e.Q.:
                        measure offsets in pre-, over-
                        measure std in pre-, over-
                        measure mean in img-

            assess std in pre- (~RON) is within allocated margins
            assess offsets in pre-, and over- are equal, within allocated margins
            assess offsets are within allocated margins
            assess injection level is within expected margins

            plot histogram of injected levels for each Q
            [plot std vs. time]

            issue any warnings to log
            issue update to report          


        """
        raise NotImplementedError
