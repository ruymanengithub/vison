#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (vertical)

Created on Tue Aug 29 17:37:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
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
from vison.support.files import cPickleRead
import TP01aux
import tptools
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

IDL=11.
IDH=18.
IG1 = 5.
IG2 = 6.5

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
                       exptime=0., shuttr=0,e_shuttr=0, 
                       mirr_on=0,
                       wave=4,
                       motr_on=0,
                       source='flat',
                       comments='')


class TP01_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
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

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(TP01, self).__init__(inputs, log, drill, debug)
        self.name = 'TP01'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = TP01aux.TP01figs.copy()
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

        TP01_commvalues['ver_shuffles'] = Nshuffles_V

        # First Injection Drain Delay

        TP01_sdict['col1'] = dict(frames=1, v_tpump=0, comments='BGD',
                                  id_dly=id_delays[0], toi_ch=toi_chinj)

        colcounter = 2
        for i, toi_tp in enumerate(toi_tpv):

            for k, vpumpmode in enumerate(vpumpmodes):
                colkey = 'col%i' % colcounter
                TP01_sdict[colkey] = dict(frames=1, toi_tp=toi_tp,
                                          id_dly=id_delays[0], v_tpmod=vpumpmode,
                                          toi_ch=toi_chinj)

                colcounter += 1

        # Second Injection Drain Delay

        TP01_sdict['col%i' % colcounter] = dict(frames=1, v_tpump=0, comments='BGD',
                                                id_dly=id_delays[1], toi_ch=toi_chinj)
        colcounter += 1

        for j, toi_tp in enumerate(toi_tpv):

            for k, vpumpmode in enumerate(vpumpmodes):

                colkey = 'col%i' % colcounter
                #print colkey
                TP01_sdict[colkey] = dict(frames=1, toi_tp=toi_tp,
                                          id_dly=id_delays[1], v_tpmod=vpumpmode,
                                          toi_ch=toi_chinj)

                colcounter += 1

        Ncols = len(TP01_sdict.keys())
        TP01_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP01_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
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
             doMask=False, # ON TESTS!
             doOffset=True, doBias=False, doFF=False)

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
        
        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='TP01 Extraction', level=0)
        
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
                        vend = self.dd.mx['vend'][ix,jCCD]


                        ioutf = 'TP01_rawmap_%i_IDDLY_%i_ROE1_%s' % (
                            ObsID, id_dly, CCDk)

                        iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ix, jCCD]

                        try:
                            iccdobj = cPickleRead(
                                os.path.join(ccdpicklespath, iccdobj_f))
                            irawmap = copy.deepcopy(iccdobj)

                            irawmap = tptools.gen_raw_dpmap_vtpump(
                                irawmap, Navgrows=-1, vstart=vstart, vend=vend)

                            irawmap.writeto(os.path.join(productspath, \
                                                         '%s.fits' % ioutf))

                            self.dd.mx['dipoles_raw'][ix, jCCD] = ioutf

                        except:  # TESTS
                            pass

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

        return # TESTS
        raise NotImplementedError

    def meta_analysis(self):
        """

        Meta-analysis of data:

            Try to identify tau and pixel-phase location for each trap.
            Need to associate dipoles across TOI_TPs and TP-patterns


        **METACODE**

        ::

            across TOI_TP, patterns:
                build catalog of traps: x,y,I-phase, Amp
                from Amp(TOI) -> tau, Pc

            Report on :
                Histogram of Taus
                Histogram of Pc (capture probability)
                Histogram of I-phases (larger phases should have more traps, 
                                  statistically) -> check

                Total Count of Traps

        """
        return # TESTS
        raise NotImplementedError
