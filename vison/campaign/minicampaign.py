#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration Campaign.

Mini-Campaign Test Set: used to verify functionality of ROE + OGSE 
ahead of running the campaign itself.

Created on Mon Jan 15 15:56:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from collections import OrderedDict
from pdb import set_trace as stop
import numpy as np
import copy

from vison.point import FOCUS00, PSF0X
from vison.dark import BIAS0X, DARK01
from vison.flat import FLAT0X
from vison.inject import CHINJ00
from vison.pump import TP00
#from vison.other import PERSIST01 as PER01
#from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
#from vison.pipe import lib as pilib
from vison.support import context
# END IMPORT


def generate_reduced_test_sequence(equipment, toGen, elvis=context.elvis, 
                                   CHAMBER=None, purpose='scripts'):
    """ """

    ogse = ogsemod.Ogse(CHAMBER)

    operator = equipment['operator']
    sn_ccd1 = equipment['sn_ccd1']
    sn_ccd2 = equipment['sn_ccd2']
    sn_ccd3 = equipment['sn_ccd3']
    sn_roe = equipment['sn_roe']
    sn_rpsu = equipment['sn_rpsu']

    test_sequence = OrderedDict()

    # BIAS

    if toGen['BIAS01']:

        print 'BIAS01...'

        Nbias01 = 3
        diffBIAS01 = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        bias01 = BIAS0X.BIAS0X(inputs=dict(
                                            test='BIAS01',
                                            N=Nbias01,
                                           diffvalues=diffBIAS01,
                                           elvis=elvis,
                                           CHAMBER=CHAMBER))

        #structBIAS01 = bias01.build_scriptdict(elvis=elvis)
        test_sequence['BIAS01'] = copy.deepcopy(bias01)

    # DARKS

    if toGen['DARK01']:

        print 'DARK01...'

        Ndark01 = 1
        exptime_dark01 = 565.  # s
        diffDARK01 = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        dark01 = DARK01.DARK01(inputs=dict(N=Ndark01, exptime=exptime_dark01,
                                           diffvalues=diffDARK01,
                                           elvis=elvis,
                                           CHAMBER=CHAMBER))

        #structDARK01 = dark01.build_scriptdict(elvis=elvis)
        test_sequence['DARK01'] = copy.deepcopy(dark01)

    # CHARGE INJECTION

    # CHINJ00

    if toGen['CHINJ00']:

        print 'CHINJ00...'

        IDH = 18.
        toi_chinj00 = 500
        chinj_on = 30
        chinj_of = 300

        diffCHINJ00 = dict(mirr_on=0, sn_ccd1=sn_ccd1,
                           sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                           sn_rpsu=sn_rpsu, operator=operator)

        chinj00 = CHINJ00.CHINJ00(inputs=dict(IDH=IDH,
                                              toi_chinj=toi_chinj00,
                                              chinj_on=chinj_on,
                                              chinj_of=chinj_of,
                                              diffvalues=diffCHINJ00,
                                              elvis=elvis,
                                              CHAMBER=CHAMBER))
        #structCHINJ00 = chinj00.build_scriptdict(elvis=elvis)
        test_sequence['CHINJ00'] = copy.deepcopy(chinj00)

    # TRAP-PUMPING

    # TP00

    if toGen['TP00']:

        print 'TP00...'

        Nshuffles_V = 500
        TOI_TPv = [200, 1000, 4000]
        Nshuffles_S = 500
        dwell_tpsv = [0, 16, 32]

        diffTP00 = dict(sn_ccd1=sn_ccd1,
                        sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                        sn_rpsu=sn_rpsu, operator=operator)

        tp00 = TP00.TP00(inputs=dict(toi_tpv=TOI_TPv,
                                     Nshuffles_V=Nshuffles_V,
                                     Nshuffles_S=Nshuffles_S,
                                     dwell_tpsv=dwell_tpsv,
                                     diffvalues=diffTP00,
                                     elvis=elvis,
                                     CHAMBER=CHAMBER))

        #structTP00 = tp00.build_scriptdict(elvis=elvis)
        test_sequence['TP00'] = copy.deepcopy(tp00)

    # FLATS

    exptimes_FLAT0X = dict(nm0=ogse.profile['tFWC_flat']['nm0'],
                           nm590=ogse.profile['tFWC_flat']['nm590'],
                           nm640=ogse.profile['tFWC_flat']['nm640'],
                           nm730=ogse.profile['tFWC_flat']['nm730'],
                           nm800=ogse.profile['tFWC_flat']['nm800'],
                           nm880=ogse.profile['tFWC_flat']['nm880'])
    # FLAT-01

    if toGen['FLAT01']:

        print 'FLAT01...'

        t_dummy_F01 = np.array([25., 50., 75])/100.
        exptimesF01 = (exptimes_FLAT0X['nm800'] * t_dummy_F01).tolist()  # s
        framesF01 = [1, 1, 1]

        diffFLAT01 = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        inpF01 = dict(exptimes=exptimesF01,
                      frames=framesF01,
                      wavelength=800,
                      test='FLAT01_800',
                      diffvalues=diffFLAT01,
                      elvis=elvis,
                      CHAMBER=CHAMBER)

        flat01 = FLAT0X.FLAT0X(inputs=inpF01)

        # structFLAT01 = flat01.build_scriptdict(
        #    diffvalues=diffFLAT01, elvis=elvis)
        test_sequence['FLAT01'] = copy.deepcopy(flat01)

    # FLAT-02

    if toGen['FLAT02']:

        wavesFLAT02 = [0, 590, 640, 730, 880]
        t_dummy_F02 = np.array([50.])/100.
        framesF02 = [2]

        diffFLAT02 = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        for iw, wave in enumerate(wavesFLAT02):

            itestkey = 'FLAT02_%i' % wave
            print '%s...' % itestkey

            iexptimesF02 = (
                exptimes_FLAT0X['nm%i' % wave] * t_dummy_F02).tolist()

            inpF02 = dict(exptimes=iexptimesF02,
                          frames=framesF02,
                          wavelength=wave,
                          test=itestkey,
                          diffvalues=diffFLAT02,
                          elvis=elvis,
                          CHAMBER=CHAMBER)

            flat02 = FLAT0X.FLAT0X(inputs=inpF02)

            # istructFLAT02 = flat02.build_scriptdict(
            #    diffvalues=diffFLAT02, elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(flat02)

    # FOCUS

    if toGen['FOCUS00']:

        print 'FOCUS00...'

        wavesFOCUS00w = [590, 800, 880]

        diffFOCUS00w = dict(sn_ccd1=sn_ccd1,
                            sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                            sn_rpsu=sn_rpsu, operator=operator)

        for iw, wave in enumerate(wavesFOCUS00w):

            tFWC_pointw = ogse.profile['tFWC_point']['nm%i' % wave]
            iexptimeF00 = 60./100. * tFWC_pointw

            itestkey = 'FOCUS00_%i' % wave

            diffFOCUS00w['test'] = itestkey

            print '%s...' % itestkey

            focus00 = FOCUS00.FOCUS00(inputs=dict(wavelength=wave,
                                                  exptime=iexptimeF00,
                                                  diffvalues=diffFOCUS00w,
                                                  elvis=elvis,
                                                  CHAMBER=CHAMBER))
            # istructFOCUS00w = focus00.build_scriptdict(diffvalues=diffFOCUS00w,
            #                                           elvis=elvis)
            test_sequence[itestkey] = copy.deepcopy(focus00)

    # PSF

    if toGen['PSF01']:

        print 'PSF01...'

        wavePSF01w = 800

        diffPSF01w = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        tFWC_pointPSF01 = ogse.profile['tFWC_point']['nm%i' % wavePSF01w]
        exptsPSF01w = (np.array([25., 50., 75.])/100. *
                       tFWC_pointPSF01).tolist()
        frsPSF01w = [4, 4, 4]

        itestkey = 'PSF01_%i' % wavePSF01w
        diffPSF01w['test'] = itestkey

        print '%s...' % itestkey

        psf01w = PSF0X.PSF0X(inputs=dict(wavelength=wave,
                                         exptimes=exptsPSF01w,
                                         frames=frsPSF01w,
                                         test=itestkey,
                                         diffvalues=diffPSF01w,
                                         elvis=elvis,
                                         CHAMBER=CHAMBER))

        # istructPSF01w = psf01w.build_scriptdict(diffvalues=diffPSF01w,
        #                                        elvis=elvis)

        test_sequence[itestkey] = copy.deepcopy(psf01w)

    return test_sequence
