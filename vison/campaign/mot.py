#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration Campaign.

M.O.T. Test Set: used to inspect PERFORMANCE at the lowest level possible

Created on Mon Jul 30 17:07:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from collections import OrderedDict
from pdb import set_trace as stop
import numpy as np
import copy

from vison.dark import BIAS0X
from vison.flat import PTC0X, BF01
from vison.inject import CHINJ01
from vison.pump import TP01, TP02
from vison.other import MOT_FF
from vison.ogse import ogse as ogsemod
from vison.support import context
# END IMPORT


def generate_mot_sequence(diffvalues, toGen, elvis=context.elvis, 
                          CHAMBER=None, purpose='scripts'):
    """ """

    ogse = ogsemod.Ogse(CHAMBER)

    test_sequence = OrderedDict()

    _toGen = dict(
            BIAS01 = False,
            CHINJ01 = False,
            TP01 = False,
            TP02 = False,
            PTC01 = False,
            MOT_FF = False
            ) 
    
    _toGen.update(toGen)

    # BIAS

    if _toGen['BIAS01']:

        print 'BIAS01...'

        Nbias01 = 5
        diffBIAS01 = dict(mirr_on=0)
        diffBIAS01.update(diffvalues)

        bias01 = BIAS0X.BIAS0X(inputs=dict(test='BIAS01',
                                           N=Nbias01,
                                           diffvalues=diffBIAS01,
                                           elvis=elvis,
                                           CHAMBER=CHAMBER))

        test_sequence['BIAS01'] = copy.deepcopy(bias01)


    # CHARGE INJECTION

    # CHINJ01

    if _toGen['CHINJ01']:

        print 'CHINJ01...'

        # CHINJ01

        IDL = 11.
        IDH = 18.
        IG2 = 7.5
        IG1s = [3., 8.]
        dIG1 = 0.75
        toi_chinj01 = 500
        id_delays = [toi_chinj01*2.5, toi_chinj01*1.5]

        diffCHINJ01 = dict(mirr_on=0)
        diffCHINJ01.update(diffvalues)

        chinj01 = CHINJ01.CHINJ01(inputs=dict(elvis=elvis,
                                              CHAMBER=CHAMBER,
                                              test='CHINJ01',
                                              IDL=IDL, IDH=IDH, 
                                              IG2=IG2,IG1s=IG1s,
                                              dIG1=dIG1,
                                              toi_chinj=toi_chinj01,
                                              id_delays=id_delays,
                                              diffvalues=diffCHINJ01))
        #structCHINJ01 = chinj01.build_scriptdict(elvis=elvis)

        test_sequence['CHINJ01'] = copy.deepcopy(chinj01)


    # TRAP-PUMPING
    
    
    # TP01

    if _toGen['TP01']:

        print 'TP01...'

        TOI_TPv = [200]
        Nshuffles_V = 5000
        toi_chinjTP01 = 250  # quick injection
        id_delays_TP01 = (np.array([2.5, 1.5]) * toi_chinjTP01).tolist()
        vpumpmodes = [123,234,341,412]

        diffTP01 = dict()
        diffTP01.update(diffvalues)

        tp01 = TP01.TP01(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP01',
                                     toi_tpv=TOI_TPv, 
                                     toi_chinj=toi_chinjTP01,
                                     Nshuffles_V=Nshuffles_V,
                                     id_delays=id_delays_TP01,
                                     vpumpmodes=vpumpmodes,
                                     diffvalues=diffTP01))
        #structTP01 = tp01.build_scriptdict(elvis=elvis)

        test_sequence['TP01'] = copy.deepcopy(tp01)

    # TP02   
    

    if _toGen['TP02']:

        print 'TP02...'

        Nshuffles_H = 5000
        dwell_sv = [0., 4.75, 14.3, 28.6]  # us
        toi_chinjTP02 = 250  # quick injection
        id_delays_TP02 = (np.array([2.5, 1.5])*toi_chinjTP02).tolist()

        diffTP02 = dict(mirr_on=0)
        diffTP02.update(diffvalues)

        tp02 = TP02.TP02(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP02',
                                     Nshuffles_H=Nshuffles_H,
                                     dwell_sv=dwell_sv, toi_chinj=toi_chinjTP02,
                                     id_delays=id_delays_TP02,
                                     diffvalues=diffTP02))

        #structTP02 = tp02.build_scriptdict(elvis=elvis)

        test_sequence['TP02'] = copy.deepcopy(tp02)

    
     # PTC
    
    if _toGen['PTC01']:

        print 'PTC01...'

        diffPTC01 = dict(mirr_on=0,
                         vstart=0,
                         vend=2086)

        diffPTC01.update(diffvalues)
        tFWC_flat800 = ogse.profile['tFWC_flat']['nm800']
        # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
        exptsPTC01 = (np.array([10., 30.,50., 70., 90.])/100.*tFWC_flat800).tolist()  # ms
        frsPTC01 = [2, 2, 2, 2, 2]

        ptc01 = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                        CHAMBER=CHAMBER,
                                        test='PTC01', exptimes=exptsPTC01,
                                        frames=frsPTC01, wavelength=800,
                                        diffvalues=diffPTC01))
        #structPTC01 = ptc01.build_scriptdict(diffvalues=diffPTC01,elvis=elvis)

        test_sequence['PTC01'] = copy.deepcopy(ptc01)

    
    if _toGen['MOT_FF']:
        
        print 'MOT_FF...'

        diffMOT_FF = dict(mirr_on=0,
                         vstart=0,
                         vend=2086)

        diffMOT_FF.update(diffvalues)
        tFWC_flat800 = ogse.profile['tFWC_flat']['nm800']
        # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
        exptsMOT_FF = (np.array([10., 30.,50., 70., 90.])/100.*tFWC_flat800).tolist()  # ms
        frsMOT_FF = [2, 2, 2, 2, 2]

        mot_ff = MOT_FF.MOT_FF(inputs=dict(elvis=elvis,
                                        CHAMBER=CHAMBER,
                                        test='MOT_FF', 
                                        surrogate='PTC01',
                                        exptimes=exptsMOT_FF,
                                        frames=frsMOT_FF, 
                                        wavelength=800,
                                        Npix=5,
                                        diffvalues=diffMOT_FF))
        #structPTC01 = ptc01.build_scriptdict(diffvalues=diffPTC01,elvis=elvis)

        test_sequence['MOT_FF'] = copy.deepcopy(mot_ff)

    
    return test_sequence
