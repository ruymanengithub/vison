#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Description of the Ground Calibration Campaign.


Created on Wed Oct 11 11:43:54 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
from collections import OrderedDict
from pdb import set_trace as stop
import numpy as np
import copy

from vison.point import FOCUS00, PSF0X
from vison.dark import BIAS01, DARK01
from vison.flat import NL01, PTC0X, FLAT0X
from vison.inject import CHINJ01, CHINJ02
from vison.pump import TP01, TP02
from vison.other import PERSIST01 as PER01
from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
#from vison.pipe import lib as pilib
from vison.support import context
# END IMPORT


def generate_test_sequence(diffvalues, toGen, elvis=context.elvis, 
                           CHAMBER=None):
    """ """
    
    ogse = ogsemod.Ogse(CHAMBER)

    print 'GENERATING TEST SEQUENCE...'

    test_sequence = OrderedDict()

    # BIAS

    if toGen['BIAS01']:

        print 'BIAS01...'

        Nbias01 = 25
        diffBIAS01 = dict(mirr_on=0)
        diffBIAS01.update(diffvalues)

        bias01 = BIAS01.BIAS01(inputs=dict(elvis=elvis,
                                           CHAMBER=CHAMBER,
                                           test='BIAS01',
                                           N=Nbias01, 
                                           diffvalues=diffBIAS01))

        #structBIAS01 = bias01.build_scriptdict(elvis=elvis)
        test_sequence['BIAS01'] = copy.deepcopy(bias01)

    # DARKS

    if toGen['DARK01']:

        print 'DARK01...'

        Ndark01 = 4
        exptime_dark01 = 565.  # s
        diffDARK01 = dict(mirr_on=0)
        diffDARK01.update(diffvalues)

        dark01 = DARK01.DARK01(inputs=dict(elvis=elvis,
                                           CHAMBER=CHAMBER,
                                           test='DARK01',
                                           N=Ndark01, 
                                           exptime=exptime_dark01,
                                           diffvalues=diffDARK01))

        #structDARK01 = dark01.build_scriptdict(elvis=elvis)
        test_sequence['DARK01'] = copy.deepcopy(dark01)

    # CHARGE INJECTION

    # CHINJ01

    if toGen['CHINJ01']:

        print 'CHINJ01...'

        # CHINJ01

        IDL = 11.
        IDH = 18.
        IG1s = [2., 6.]
        toi_chinj01 = 500
        id_delays = [toi_chinj01*2.5, toi_chinj01*1.5]

        diffCHINJ01 = dict(mirr_on=0)
        diffCHINJ01.update(diffvalues)

        chinj01 = CHINJ01.CHINJ01(inputs=dict(elvis=elvis,
                                              CHAMBER=CHAMBER,
                                              test='CHINJ01',
                                              IDL=IDL, IDH=IDH, IG1s=IG1s,
                                              toi_chinj=toi_chinj01,
                                              id_delays=id_delays,
                                              diffvalues=diffCHINJ01))
        #structCHINJ01 = chinj01.build_scriptdict(elvis=elvis)

        test_sequence['CHINJ01'] = copy.deepcopy(chinj01)

    # CHINJ02

    if toGen['CHINJ02']:

        print 'CHINJ02...'

        IDLs = [10., 13.]
        IDH = 18.
        #IG1s = [2000,6000]
        toi_chinj02 = 500
        id_delays = [toi_chinj02*2.5, toi_chinj02*1.5]
        diffCHINJ02 = dict(mirr_on=0)
        diffCHINJ02.update(diffvalues)

        chinj02 = CHINJ02.CHINJ02(inputs=dict(elvis=elvis,
                                              CHAMBER=CHAMBER,
                                              test='CHINJ02',
                                              IDLs=IDLs, IDH=IDH, toi_chinj=toi_chinj02,
                                              id_delays=id_delays, diffvalues=diffCHINJ02))
        #structCHINJ02 = chinj02.build_scriptdict(elvis=elvis)
        test_sequence['CHINJ02'] = copy.deepcopy(chinj02)

    # TRAP-PUMPING

    # TP01

    if toGen['TP01']:

        print 'TP01...'

        TOI_TPv = [200, 1000, 2000, 4000, 8000]
        toi_chinjTP01 = 250  # quick injection
        id_delays_TP01 = (np.array([2.5, 1.5]) * toi_chinjTP01).tolist()

        diffTP01 = dict()
        diffTP01.update(diffvalues)

        tp01 = TP01.TP01(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP01',
                                     toi_tpv=TOI_TPv, toi_chinj=toi_chinjTP01,
                                     id_delays=id_delays_TP01,
                                     diffvalues=diffTP01))
        #structTP01 = tp01.build_scriptdict(elvis=elvis)

        test_sequence['TP01'] = copy.deepcopy(tp01)

    # TP02

    if toGen['TP02']:

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

    # FLATS

    exptimes_FLAT0X = dict(nm590=ogse.profile['tFWC_flat']['nm590'],
                           nm640=ogse.profile['tFWC_flat']['nm640'],
                           nm730=ogse.profile['tFWC_flat']['nm730'],
                           nm800=ogse.profile['tFWC_flat']['nm800'],
                           nm880=ogse.profile['tFWC_flat']['nm880'])
    # FLAT-01

    if toGen['FLAT01']:

        print 'FLAT01...'

        t_dummy_F01 = np.array([25., 50., 75])/100.
        exptimesF01 = (exptimes_FLAT0X['nm800'] * t_dummy_F01).tolist()  # s
        framesF01 = [80, 60, 30]

        diffFLAT01 = dict(mirr_on=0)
        diffFLAT01.update(diffvalues)

        inpF01 = dict(elvis=elvis,
                      CHAMBER=CHAMBER,
                      test='FLAT01',
                      exptimes=exptimesF01,
                      frames=framesF01,
                      wavelength=800,
                      diffvalues=diffFLAT01)

        flat01 = FLAT0X.FLAT0X(inputs=inpF01)

        #structFLAT01 = flat01.build_scriptdict(diffvalues=diffFLAT01,elvis=elvis)
        test_sequence['FLAT01'] = copy.deepcopy(flat01)

    # FLAT-02

    if toGen['FLAT02']:

        wavesFLAT02 = [590, 640, 880]
        t_dummy_F02 = np.array([25., 75])/100.
        framesF02 = [80, 30]

        diffFLAT02 = dict(mirr_on=0)
        diffFLAT02.update(diffvalues)

        for iw, wave in enumerate(wavesFLAT02):

            itestkey = 'FLAT02_%i' % wave
            print '%s...' % itestkey

            iexptimesF02 = (
                exptimes_FLAT0X['nm%i' % wave] * t_dummy_F02).tolist()

            inpF02 = dict(elvis=elvis,
                          CHAMBER=CHAMBER,
                          exptimes=iexptimesF02,
                          frames=framesF02,
                          wavelength=wave,
                          test=itestkey,
                          diffvalues=diffFLAT02)

            flat02 = FLAT0X.FLAT0X(inputs=inpF02)

            #istructFLAT02 = flat02.build_scriptdict(diffvalues=diffFLAT02,elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(flat02)

    # PTC

    # PTC-01

    if toGen['PTC01']:

        print 'PTC01...'

        diffPTC01 = dict(mirr_on=0,
                         vstart=0,
                         vend=2086)

        diffPTC01.update(diffvalues)
        tFWC_flat800 = ogse.profile['tFWC_flat']['nm800']
        # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
        exptsPTC01 = (np.array([5., 10., 20., 30., 50., 70., 80., 90.,
                                100., 110., 120.])/100.*tFWC_flat800).tolist()  # ms
        frsPTC01 = [10, 10, 10, 10, 10, 10, 10, 10, 4, 4, 4]

        ptc01 = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                        CHAMBER=CHAMBER,
                                        test='PTC01', exptimes=exptsPTC01,
                                        frames=frsPTC01, wavelength=800,
                                        diffvalues=diffPTC01))
        #structPTC01 = ptc01.build_scriptdict(diffvalues=diffPTC01,elvis=elvis)

        test_sequence['PTC01'] = copy.deepcopy(ptc01)

    # PTC-02 - wavelength

    if toGen['PTC02WAVE']:

        print 'PTC02WAVE...'

        wavesPTC02w = [590, 640, 730, 880]

        diffPTC02w = dict(mirr_on=0)
        diffPTC02w.update(diffvalues)

        for iw, wave in enumerate(wavesPTC02w):

            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
            
            tFWC_flatw = ogse.profile['tFWC_flat']['nm%i' % wave]

            exptsPTC02w = (np.array(
                [10., 30., 50., 70., 80., 90.])/100.*tFWC_flatw).tolist()
            frsPTC02w = [4, 4, 4, 4, 4, 4]

            itestkey = 'PTC02_%i' % wave
            diffPTC02w['test'] = itestkey

            print '%s...' % itestkey

            ptc02w = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                             CHAMBER=CHAMBER,
                                             test=itestkey,
                                             wavelength=wave,
                                             exptimes=exptsPTC02w,
                                             frames=frsPTC02w,
                                             diffvalues=diffPTC02w))

            #istructPTC02w = ptc02w.build_scriptdict(diffvalues=diffPTC02w,elvis=elvis)
            test_sequence[itestkey] = copy.deepcopy(ptc02w)

   # PTC-02 - Temp.

    if toGen['PTC02TEMP']:

        print 'PTC02TEMP...'

        wavePTC02T = 800
        TempsPTC02T = [150., 156.]

        diffPTC02T = dict(mirr_on=0)
        diffPTC02T.update(diffvalues)

        # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
        tFWC_flatw = ogse.profile['tFWC_flat']['nm%i' % wavePTC02T]
        exptsPTC02T = (np.array([10., 30., 50., 70., 80., 90.]) /
                       100.*tFWC_flatw).tolist()
        frsPTC02T = [4, 4, 4, 4, 4, 4]

        for it, T in enumerate(TempsPTC02T):

            itestkey = 'PTC02_%iK' % T
            print '%s...' % itestkey

            diffPTC02T['test'] = itestkey

            ptc02t = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                             CHAMBER=CHAMBER,
                                             test=itestkey,
                                             wavelength=wavePTC02T,
                                             frames=frsPTC02T,
                                             exptimes=exptsPTC02T,
                                             diffvalues=diffPTC02T))

            # istructPTC02T = ptc02t.build_scriptdict(diffvalues=diffPTC02T,
            #                                        elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(ptc02t)

    # NL-01

    if toGen['NL01']:

        print 'NL01...'
        waveNL01 = 880

        diffNL01 = dict(mirr_on=0)
        diffNL01.update(diffvalues)

        # 5 frames per fluence: 1%, 2%, 3%, 5%, 10%, 20%,30%, 50%,70%,80%,85%,90%,95%
        tFWC_flatNL01 = ogse.profile['tFWC_flat']['nm%i' % waveNL01]
        exptsNL01 = (np.array([5., 10., 20., 30., 50., 70., 80., 90., 100.,
                               110., 120.])/100. * tFWC_flatNL01).tolist()  # ms
        exptinterNL01 = 0.5 * tFWC_flatNL01
        frsNL01 = (np.ones(11, dtype='int32')*5).tolist()

        nl01 = NL01.NL01(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='NL01',
                                     exptimes=exptsNL01, exptinter=exptinterNL01,
                                     frames=frsNL01, wavelength=waveNL01,
                                     diffvalues=diffNL01))

        #structNL01 = nl01.build_scriptdict(diffvalues=diffNL01,elvis=elvis)
        test_sequence['NL01'] = copy.deepcopy(nl01)

    # FOCUS

    if toGen['FOCUS00']:

        print 'FOCUS00...'

        # wavesFOCUS00w = [590,640,730,800,880] # TESTS
        wavesFOCUS00w = [800]

        diffFOCUS00w = dict()
        diffFOCUS00w.update(diffvalues)

        for iw, wave in enumerate(wavesFOCUS00w):
            
            tFWC_pointw = ogse.profile['tFWC_point']['nm%i' % wave]

            iexptimeF00 = 60./100. * tFWC_pointw

            itestkey = 'FOCUS00_%i' % wave

            diffFOCUS00w['test'] = itestkey

            print '%s...' % itestkey

            focus00 = FOCUS00.FOCUS00(inputs=dict(elvis=elvis,
                                                  CHAMBER=CHAMBER,
                                                  test=itestkey,
                                                  wavelength=wave,
                                                  exptime=iexptimeF00,
                                                  diffvalues=diffFOCUS00w))
            # istructFOCUS00w = focus00.build_scriptdict(diffvalues=diffFOCUS00w,
            #                                                   elvis=elvis)
            test_sequence[itestkey] = copy.deepcopy(focus00)

    # PSF

    if toGen['PSF01']:

        print 'PSF01...'

        #wavesPSF01w = [590,640,800,880]
        PSF0Xtestdefaults = PSF0X.get_testdefaults(ogse)
        wavesPSF01w = PSF0Xtestdefaults['waves']

        diffPSF01w = dict()
        diffPSF01w.update(diffvalues)

        for iw, wave in enumerate(wavesPSF01w):

            exptsPSF01w = PSF0Xtestdefaults['exptimes']['nm%i' % wave]
            frsPSF01w = PSF0Xtestdefaults['frames']

            itestkey = 'PSF01_%i' % wave
            diffPSF01w['test'] = itestkey

            print '%s...' % itestkey

            psf01w = PSF0X.PSF0X(inputs=dict(elvis=elvis,
                                             CHAMBER=CHAMBER,
                                             wavelength=wave,
                                             exptimes=exptsPSF01w,
                                             frames=frsPSF01w,
                                             test=itestkey,
                                             diffvalues=diffPSF01w))

            # istructPSF01w = psf01w.build_scriptdict(diffvalues=diffPSF01w,
            #                                       elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(psf01w)

    # PSF02

    if toGen['PSF02']:

        print 'PSF02...'

        wPSF02 = 800
        temps_PSF02 = [150, 156]
        
        tFWC_pointPSF02  = ogse.profile['tFWC_point']['nm%i' % wPSF02]

        exptsPSF02 = np.array([5., 25., 50., 75., 90.]) / \
            100. * tFWC_pointPSF02
        frsPSF02 = [20, 14, 10, 4, 4]

        diffPSF02 = dict()
        diffPSF02.update(diffvalues)

        for it, temp in enumerate(temps_PSF02):

            itestkey = 'PSF02_%iK' % temp

            diffPSF02['test'] = itestkey

            print '%s...' % itestkey

            psf02k = PSF0X.PSF0X(inputs=dict(elvis=elvis,
                                             CHAMBER=CHAMBER,
                                             wavelength=800,
                                             exptimes=exptsPSF02.tolist(),
                                             frames=frsPSF02,
                                             test=itestkey,
                                             diffvalues=diffPSF02))

            # istructPSF02 = psf02k.build_scriptdict(diffvalues=diffPSF02.copy(),
            #                                       elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(psf02k)

    # OTHER

    # PERSIST

    if toGen['PERSIST01']:

        print 'PERSIST01...'
        
        tFWC_point590 = ogse.profile['tFWC_point']['nm590']

        exptPER01_SATUR = tFWC_point590*100.   # s
        exptPER01_LATEN = 565.  # s

        diffPER01 = dict()
        diffPER01.update(diffvalues)

        persist01 = PER01.PERSIST01(inputs=dict(elvis=elvis,
                                                CHAMBER=CHAMBER,
                                                test='PERSIST01',
                                                exptSATUR=exptPER01_SATUR,
                                                exptLATEN=exptPER01_LATEN,
                                                diffvalues=diffPER01))

        #structPER01 = persist01.build_scriptdict(diffvalues=diffPER01,elvis=elvis)

        test_sequence['PERSIST01'] = copy.deepcopy(persist01)

    return test_sequence
