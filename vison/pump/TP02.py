#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: TP01

Trap-Pumping calibration (serial)

Created on Tue Aug 29 17:38:00 2017

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
#from vison.pipe.task import Task
from PumpTask import PumpTask
from vison.image import performance
from vison.datamodel import inputs
from vison.support.files import cPickleRead
import TP02aux
import tptools
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

IG1 = 6
IG2 = 5

TP02_commvalues = dict(program='CALCAMP', test='TP02',
                       exptime=0., shuttr=0,
                       vstart=0, vend=100,
                       siflsh=1, siflsh_p=500,
                       IDL=11, IDH=18,
                       IG1_1_T=IG1, IG1_2_T=IG1, IG1_3_T=IG1,
                       IG1_1_B=IG1, IG1_2_B=IG1, IG1_3_B=IG1,
                       IG2_T=IG2, IG2_B=IG2,
                       chinj=1, chinj_on=2066, chinj_of=0,
                       chin_dly=0,
                       id_wid=60,
                       toi_chinj=500,
                       s_tpump=1, s_tp_cnt=5000,
                       v_tp_cnt=0, dwell_v=0, dwell_s=0,
                       comments='')


class TP02_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
        ('Nshuffles_H',
         ([int], 'Number of Shuffles, Horizontal/Serial Pumping.')),
        ('dwell_sv', ([list], 'Dwell Times list [serial].')),
        ('id_delays',
         ([list], 'Injection Drain Delays [2, one per CCDs section].')),
        ('toi_chinj', ([int], 'TOI Charge Injection.')),
        ('spumpmodes',
         ([list], 'Horizontal/Serial Pumping Starting points.'))
    ])))


class TP02(PumpTask):
    """ """

    inputsclass = TP02_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data),
                         ('prep', self.prepare_images),
                         ('extract', self.extract),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(TP02, self).__init__(inputs, log, drill, debug)
        self.name = 'TP02'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = TP02aux.TP02figs.copy()
        self.inputs['subpaths'] = dict(figs='figs', 
                   ccdpickles='ccdpickles',
                   products='products')

    def set_inpdefaults(self, **kwargs):
        """ """
        toi_chinj = 500
        self.inpdefaults = dict(toi_chinj=toi_chinj,
                                Nshuffles_H=5000,
                                dwell_sv=[0., 4.75, 14.3, 28.6],
                                id_delays=np.array([2.5, 1.5])*toi_chinj,
                                spumpmodes=[23, 31])

    def set_perfdefaults(self, **kwargs):
        super(TP02, self).set_perfdefaults(**kwargs)

        Flu_lims, FluGrad_lims = self.get_FluenceAndGradient_limits()

        self.perfdefaults['Flu_lims'] = Flu_lims.copy()
        self.perfdefaults['FluGrad_lims'] = FluGrad_lims.copy()

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """ """

        Nshuffles_H = self.inputs['Nshuffles_H']
        dwell_sv = self.inputs['dwell_sv']
        id_delays = self.inputs['id_delays']
        toi_chinj = self.inputs['toi_chinj']
        spumpmodes = self.inputs['spumpmodes']

        assert len(id_delays) == 2

        TP02_sdict = dict()

        TP02_commvalues['ser_shuffles'] = Nshuffles_H

        # First Injection Drain Delay

        TP02_sdict['col1'] = dict(frames=1, v_tpump=0, s_tpump=0,
                                  comments='BGD', id_dly=id_delays[0], toi_ch=toi_chinj)

        colcounter = 2
        for i, dwell_s in enumerate(dwell_sv):

            for k, sermode in enumerate(spumpmodes):
                colkey = 'col%i' % colcounter
                TP02_sdict[colkey] = dict(frames=1, dwell_s=dwell_s,
                                          id_dly=id_delays[0], s_tpmod=sermode, toi_ch=toi_chinj)

                colcounter += 1

        # Second Injection Drain Delay

        TP02_sdict['col%i' % colcounter] = dict(frames=1, v_tpump=0, s_tpump=0,
                                                comments='BGD', id_dly=id_delays[1], toi_ch=toi_chinj)
        colcounter += 1

        for j, dwell_s in enumerate(dwell_sv):

            for k, sermode in enumerate(spumpmodes):

                colkey = 'col%i' % colcounter
                #print colkey
                TP02_sdict[colkey] = dict(frames=1, dwell_s=dwell_s,
                                          id_dly=id_delays[1], s_tpmod=sermode, toi_ch=toi_chinj)

                colcounter += 1

        Ncols = len(TP02_sdict.keys())
        TP02_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(TP02_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        TP02_sdict = sc.update_structdict(TP02_sdict, commvalues, diffvalues)

        return TP02_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(TP02, self).filterexposures(structure, explog, OBSID_lims, colorblind=True,
                                                 wavedkeys=wavedkeys)
    
    def prepare_images(self):
        super(TP02, self).prepare_images(doExtract=True, 
             doMask=False, # ON TESTS!
             doOffset=True, doBias=False, doFF=False)
    
    def extract(self):
        """ 
        Obtain Maps of Serial Dipoles.
        
        """
        
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
                    
                    ixsel = np.where((self.dd.mx['id_dly'][:,jCCD] == id_dly) & (
                        self.dd.mx['s_tpump'][:,jCCD] != 0))
                    ixref = np.where((self.dd.mx['id_dly'][:,jCCD] == id_dly) & (
                        self.dd.mx['s_tpump'][:,jCCD] == 0))
                    
                    Rvstart = self.dd.mx['vstart'][ixref[0][0], jCCD]
                    Rvend = self.dd.mx['vend'][ixref[0][0], jCCD]
                    
                    ccdref_f = '%s.pick' % self.dd.mx['ccdobj_name'][ixref[0][0], jCCD]
                    ccdref = cPickleRead( os.path.join(ccdpicklespath, ccdref_f))
                    
                    InjProfiles = OrderedDict()
                    
                    for Q in ccdref.Quads:
                        InjProfiles[Q] = tptools.get_InjProfile(ccdref, 
                                   Q, Navgrows=-1, vstart=Rvstart, vend=Rvend, 
                                   extension=-1)

                    for ix in ixsel[0]:
                        
                        ObsID = self.dd.mx['ObsID'][ix]
                        vstart = self.dd.mx['vstart'][ix, jCCD]
                        vend = self.dd.mx['vend'][ix,jCCD]

                        ioutf = 'TP02_rawmap_%i_IDDLY_%i_ROE1_%s' % (
                            ObsID, id_dly, CCDk)

                        iccdobj_f = '%s.pick' % self.dd.mx['ccdobj_name'][ix, jCCD]

                        try:
                            iccdobj = cPickleRead(
                                os.path.join(ccdpicklespath, iccdobj_f))
                            irawmap = copy.deepcopy(iccdobj)


                            irawmap = tptools.gen_raw_dpmap_stpump(
                                irawmap, InjProfiles[Q], vstart=vstart, vend=vend)

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
                        load raw 1D map of relative pumping (from extract_data)
                        identify dipoles:
                            x, rel-amplitude, orientation (E or W)

            produce & report:  
                map location of dipoles
                PDF of dipole amplitudes (for E and W)
                Counts of dipoles (and E vs. W)

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
                build catalog of traps: x,y,R-phase, amp(dwell)
                from Amp(dwell) -> tau, Pc

            Report on :
               Histogram of Taus
               Histogram of Pc (capture probability)
               Histogram of R-phases

               Total Count of Traps



        """
        return #TESTS
        raise NotImplementedError
