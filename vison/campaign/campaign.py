#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Description of the Ground Calibration Campaign.

:History:
Created on Wed Oct 11 11:43:54 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from collections import OrderedDict
from pdb import set_trace as stop
import numpy as np
import copy

from vison.point import FOCUS00, PSF0X
from vison.dark import BIAS0X, DARK01
from vison.flat import NL01, NL02, PTC0X, FLAT0X, BF01
from vison.inject import CHINJ01, CHINJ02
from vison.pump import TP01, TP11, TP02, TP21
from vison.other import PERSIST01 as PER01
from vison.other import MOT_WARM, COSMETICS00
from vison.point import lib as polib
from vison.ogse import ogse as ogsemod
#from vison.pipe import lib as pilib
from vison.support import context
from vison.support import utils
# END IMPORT


def generate_test_sequence(diffvalues, toGen, elvis=context.elvis,
                           CHAMBER=None, purpose='scripts'):
    """
    | Function that generates a number of tests, as instances of their corresponding
    task classes. 

    | Aimed at the VGCC.

    | Now supporting test repetitions."""
    taskslist = toGen.keys()

    test_sequence = OrderedDict()
    for taskname in taskslist:
        if not toGen[taskname]:
            continue

        strip_taskname, iteration = utils.remove_iter_tag(taskname, Full=True)
        _toGen = OrderedDict()
        _toGen[strip_taskname] = True
        ans = _generate_test_sequence(diffvalues, _toGen,
                                      elvis, CHAMBER, purpose)

        if iteration is not None:
            for key in list(ans.keys()):
                test_sequence['%s.%i' % (key, iteration)] = copy.deepcopy(ans[key])
        else:
            for key in list(ans.keys()):
                test_sequence[key] = copy.deepcopy(ans[key])

    return test_sequence


def _generate_test_sequence(diffvalues, toGen, elvis=context.elvis,
                            CHAMBER=None, purpose='scripts'):
    """ """

    ogse = ogsemod.Ogse(CHAMBER)

    #print 'GENERATING TEST SEQUENCE...'

    test_sequence = OrderedDict()

    _toGen = dict(
        COSMETICS00=False,
        BIAS01=False,
        BIAS02=False,
        DARK01=False,
        CHINJ01=False,
        CHINJ02=False,
        TP01=False,
        TP11=False,
        TP02=False,
        TP21=False,
        FLATFLUX00=False,
        FLAT01=False,
        FLAT02=False,
        FLAT_STB=False,
        PTC01=False,
        PTC02WAVE=False,
        PTC02TEMP=False,
        PTC02RD=False,
        NL01=False,
        NL02=False,
        NL02RD=False,
        FOCUS00=False,
        PSF01=False,
        PSFLUX00=False,
        PSF02=False,
        PERSIST01=False,
        MOT_WARM=False,
        BF01=False,
        BF01WAVE=False,
        MOT_FF=False)

    _toGen.update(toGen)

    # BIAS

    if _toGen['BIAS01']:

        Nbias01 = 10
        diffBIAS01 = dict(mirr_on=0)
        diffBIAS01.update(diffvalues)

        bias01 = BIAS0X.BIAS0X(inputs=dict(elvis=elvis,
                                           CHAMBER=CHAMBER,
                                           test='BIAS01',
                                           N=Nbias01,
                                           diffvalues=diffBIAS01))

        #structBIAS01 = bias01.build_scriptdict(elvis=elvis)
        test_sequence['BIAS01'] = copy.deepcopy(bias01)

    # BIAS

    if _toGen['BIAS02']:

        Nbias02 = 10
        diffBIAS02 = dict(mirr_on=0,
                          rdmode='rwd_bas_v',
                          swellw=context.sumwell['rwd_bas_v'][0],
                          swelldly=context.sumwell['rwd_bas_v'][1])
        diffBIAS02.update(diffvalues)

        bias02 = BIAS0X.BIAS0X(inputs=dict(elvis=elvis,
                                           CHAMBER=CHAMBER,
                                           test='BIAS02',
                                           N=Nbias02,
                                           diffvalues=diffBIAS02))

        #structBIAS01 = bias01.build_scriptdict(elvis=elvis)
        test_sequence['BIAS02'] = copy.deepcopy(bias02)

    # DARKS

    if _toGen['DARK01']:

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

    if _toGen['CHINJ01']:

        # CHINJ01

        IDL = 11.
        IDH = 18.
        IG2 = 6.00
        IG1s = [2.5, 6.75]
        toi_chinj01 = 500
        id_delays = [toi_chinj01 * 2.5, toi_chinj01 * 1.5]  # bottom [EF], top [GH] CCD half

        diffCHINJ01 = dict(mirr_on=0)
        diffCHINJ01.update(diffvalues)

        chinj01 = CHINJ01.CHINJ01(inputs=dict(elvis=elvis,
                                              CHAMBER=CHAMBER,
                                              test='CHINJ01',
                                              IDL=IDL, IDH=IDH,
                                              IG2=IG2, IG1s=IG1s,
                                              toi_chinj=toi_chinj01,
                                              id_delays=id_delays,
                                              diffvalues=diffCHINJ01))
        #structCHINJ01 = chinj01.build_scriptdict(elvis=elvis)

        test_sequence['CHINJ01'] = copy.deepcopy(chinj01)

    # CHINJ02

    if _toGen['CHINJ02']:

        IDLs = [9.5, 12.5]
        IDH = 18.
        dIDL = 0.25  # V
        IG1 = 4.
        IG2 = 6.
        toi_chinj02 = 500
        id_delays = [toi_chinj02 * 2.5, toi_chinj02 * 1.5]  # bottom [EF], top [GH] CCD half
        diffCHINJ02 = dict(mirr_on=0)
        diffCHINJ02.update(diffvalues)

        chinj02 = CHINJ02.CHINJ02(inputs=dict(elvis=elvis,
                                              CHAMBER=CHAMBER,
                                              test='CHINJ02',
                                              IDLs=IDLs, dIDL=dIDL,
                                              IDH=IDH, IG1=IG1, IG2=IG2,
                                              toi_chinj=toi_chinj02,
                                              id_delays=id_delays,
                                              diffvalues=diffCHINJ02))
        #structCHINJ02 = chinj02.build_scriptdict(elvis=elvis)
        test_sequence['CHINJ02'] = copy.deepcopy(chinj02)

    # TRAP-PUMPING

    # TP01

    if _toGen['TP01']:

        toi_chinjTP01 = 250  # quick injection
        Nshuffles_V = 5000
        TOI_TPv = [50, 200, 1000, 2000, 4000, 8000]
        vpumpmodes = [123, 234, 341, 412]

        # top, bottom CCD half
        id_delays_TP01 = (np.array([2.5, 1.5]) * toi_chinjTP01).tolist()

        diffTP01 = dict()
        diffTP01.update(diffvalues)

        tp01 = TP01.TP01(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP01',
                                     toi_chinj=toi_chinjTP01,
                                     Nshuffles_V=Nshuffles_V,
                                     id_delays=id_delays_TP01,
                                     toi_tpv=TOI_TPv,
                                     vpumpmodes=vpumpmodes,
                                     diffvalues=diffTP01))
        #structTP01 = tp01.build_scriptdict(elvis=elvis)

        test_sequence['TP01'] = copy.deepcopy(tp01)

    # TP11

    if _toGen['TP11']:

        toi_chinjTP11 = 250  # quick injection
        Nshuffles_V = 5000
        TOI_TPv = [50, 200, 1000, 2000, 4000, 8000]
        vpumpmodes = [123, 234, 341, 412]

        # top, bottom CCD half
        id_delays_TP11 = (np.array([2.5, 1.5]) * toi_chinjTP11).tolist()

        diffTP11 = dict()
        diffTP11.update(diffvalues)

        tp11 = TP11.TP11(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP11',
                                     toi_chinj=toi_chinjTP11,
                                     Nshuffles_V=Nshuffles_V,
                                     id_delays=id_delays_TP11,
                                     toi_tpv=TOI_TPv,
                                     vpumpmodes=vpumpmodes,
                                     diffvalues=diffTP11))

        test_sequence['TP11'] = copy.deepcopy(tp11)

    # TP02

    if _toGen['TP02']:

        Nshuffles_H = 5000
        dwell_sv = [0., 4.75, 14.3, 28.6]  # us
        toi_chinjTP02 = 250  # quick injection
        # top, bottom CCD half
        id_delays_TP02 = (np.array([2.5, 1.5]) * toi_chinjTP02).tolist()
        spumpmodes = [23, 31]

        diffTP02 = dict(mirr_on=0)
        diffTP02.update(diffvalues)

        tp02 = TP02.TP02(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP02',
                                     toi_chinj=toi_chinjTP02,
                                     Nshuffles_H=Nshuffles_H,
                                     dwell_sv=dwell_sv,
                                     id_delays=id_delays_TP02,
                                     spumpmodes=spumpmodes,
                                     diffvalues=diffTP02))

        #structTP02 = tp02.build_scriptdict(elvis=elvis)

        test_sequence['TP02'] = copy.deepcopy(tp02)

    # TP21

    if _toGen['TP21']:

        Nshuffles_H = 5000
        dwell_sv = [0., 4.75, 14.3, 28.6]  # us
        toi_chinjTP21 = 250  # quick injection
        # top, bottom CCD half
        id_delays_TP21 = (np.array([2.5, 1.5]) * toi_chinjTP21).tolist()
        spumpmodes = [23, 31]

        diffTP21 = dict(mirr_on=0)
        diffTP21.update(diffvalues)

        tp21 = TP21.TP21(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='TP21',
                                     toi_chinj=toi_chinjTP21,
                                     Nshuffles_H=Nshuffles_H,
                                     dwell_sv=dwell_sv,
                                     id_delays=id_delays_TP21,
                                     spumpmodes=spumpmodes,
                                     diffvalues=diffTP21))

        #structTP02 = tp02.build_scriptdict(elvis=elvis)

        test_sequence['TP21'] = copy.deepcopy(tp21)

    # FLATS

    exptimes_FLAT0X = dict(nm590=ogse.profile['tFWC_flat']['nm590'],
                           nm2000=ogse.profile['tFWC_flat']['nm2000'],
                           nm730=ogse.profile['tFWC_flat']['nm730'],
                           nm800=ogse.profile['tFWC_flat']['nm800'],
                           nm880=ogse.profile['tFWC_flat']['nm880'],
                           nm0=ogse.profile['tFWC_flat']['nm0'])

    if _toGen['FLATFLUX00']:

        wavesFLATFLUX00 = [590, 2000, 730, 800, 880, 0]

        diffFLATFLUX00w = dict()
        diffFLATFLUX00w.update(diffvalues)

        for iw, wave in enumerate(wavesFLATFLUX00):

            tFWC_flatw = exptimes_FLAT0X['nm%i' % wave]

            diffFLATFLUX00w = dict(mirr_on=0)
            diffFLATFLUX00w.update(diffvalues)

            exptsFLATFLUX00w = (np.array(
                [5., 20., 50., 80.]) / 100. * tFWC_flatw).tolist()
            frsFLATFLUX00w = [1, 1, 1, 1]

            itestkey = 'FLATFLUX00_%i' % wave
            diffFLATFLUX00w['test'] = itestkey

            flatflux00w = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                                  CHAMBER=CHAMBER,
                                                  test=itestkey,
                                                  exptimes=exptsFLATFLUX00w,
                                                  frames=frsFLATFLUX00w,
                                                  wavelength=wave,
                                                  diffvalues=diffFLATFLUX00w))

            test_sequence[itestkey] = copy.deepcopy(flatflux00w)

    # FLAT-01

    if _toGen['FLAT01']:

        t_dummy_F01 = np.array([25., 50., 75]) / 100.
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

    if _toGen['FLAT02']:

        wavesFLAT02 = [590, 730, 880]
        t_dummy_F02 = np.array([25., 75]) / 100.
        framesF02 = [80, 30]

        diffFLAT02 = dict(mirr_on=0)
        diffFLAT02.update(diffvalues)

        for iw, wave in enumerate(wavesFLAT02):

            itestkey = 'FLAT02_%i' % wave

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

    # FLAT-STABILITY

    if _toGen['FLAT_STB']:

        waveFLATSTB = 800
        t_dummy_FSTB = np.array([30.]) / 100.
        framesFST = [15]

        diffFLATSTB = dict(mirr_on=0)
        diffFLATSTB.update(diffvalues)

        testkey = 'FLAT_STB'

        exptimesFSTB = (
            exptimes_FLAT0X['nm%i' % waveFLATSTB] * t_dummy_FSTB).tolist()

        inpFSTB = dict(elvis=elvis,
                       CHAMBER=CHAMBER,
                       exptimes=exptimesFSTB,
                       frames=framesFST,
                       wavelength=waveFLATSTB,
                       test=testkey,
                       diffvalues=diffFLATSTB)

        flatstb = FLAT0X.FLAT0X(inputs=inpFSTB)

        #istructFLAT02 = flat02.build_scriptdict(diffvalues=diffFLAT02,elvis=elvis)

        test_sequence[testkey] = copy.deepcopy(flatstb)

    # PTC

    # PTC-01

    if _toGen['PTC01']:

        diffPTC01 = dict(mirr_on=0,
                         vstart=0,
                         vend=2086)

        diffPTC01.update(diffvalues)
        tFWC_flat800 = ogse.profile['tFWC_flat']['nm800']
        # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
        exptsPTC01 = (np.array([5., 10., 20., 30., 50., 70., 80., 90.,
                                100., 110., 120.]) / 100. * tFWC_flat800).tolist()  # ms
        frsPTC01 = [10, 10, 10, 10, 10, 10, 10, 10, 4, 4, 4]

        ptc01 = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                        CHAMBER=CHAMBER,
                                        test='PTC01',
                                        exptimes=exptsPTC01,
                                        frames=frsPTC01,
                                        wavelength=800,
                                        diffvalues=diffPTC01))
        #structPTC01 = ptc01.build_scriptdict(diffvalues=diffPTC01,elvis=elvis)

        test_sequence['PTC01'] = copy.deepcopy(ptc01)

    # PTC-02 - wavelength

    if _toGen['PTC02WAVE']:

        wavesPTC02w = [590, 730, 880]

        diffPTC02w = dict(mirr_on=0)
        diffPTC02w.update(diffvalues)

        for iw, wave in enumerate(wavesPTC02w):

            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.

            tFWC_flatw = ogse.profile['tFWC_flat']['nm%i' % wave]

            exptsPTC02w = (np.array(
                [10., 30., 50., 70., 80., 90.]) / 100. * tFWC_flatw).tolist()
            frsPTC02w = [4, 4, 4, 4, 4, 4]

            itestkey = 'PTC02_%i' % wave
            diffPTC02w['test'] = itestkey

            ptc02w = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                             CHAMBER=CHAMBER,
                                             test=itestkey,
                                             exptimes=exptsPTC02w,
                                             frames=frsPTC02w,
                                             wavelength=wave,
                                             diffvalues=diffPTC02w))

            #istructPTC02w = ptc02w.build_scriptdict(diffvalues=diffPTC02w,elvis=elvis)
            test_sequence[itestkey] = copy.deepcopy(ptc02w)

    if _toGen['BF01']:

        diffBF01 = dict(mirr_on=0,
                        vstart=0,
                        vend=2086)

        diffBF01.update(diffvalues)
        tFWC_flat800 = ogse.profile['tFWC_flat']['nm800']
        # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
        exptsBF01 = (np.array([5., 10., 20., 30., 50., 70., 80., 90.,
                               100., 110., 120.]) / 100. * tFWC_flat800).tolist()  # ms
        frsBF01 = [10, 10, 10, 10, 10, 10, 10, 10, 4, 4, 4]

        bf01 = BF01.BF01(inputs=dict(
            elvis=elvis,
            CHAMBER=CHAMBER,
            test='BF01',
            exptimes=exptsBF01,
            frames=frsBF01,
            wavelength=800,
            Npix=5,
            clipsigma=4.,
            covfunc='ver2',
            doBiasCorr=True,
            central='mean',
            surrogate='PTC01',
            diffvalues=diffBF01))

        test_sequence['BF01'] = copy.deepcopy(bf01)

    if _toGen['BF01WAVE']:

        wavesBF01w = [590, 730, 880]

        diffBF01w = dict(mirr_on=0)
        diffBF01w.update(diffvalues)

        for iw, wave in enumerate(wavesBF01w):

            # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.

            tFWC_flatw = ogse.profile['tFWC_flat']['nm%i' % wave]

            exptsBF01w = (np.array(
                [10., 30., 50., 70., 80., 90.]) / 100. * tFWC_flatw).tolist()
            frsBF01w = [4, 4, 4, 4, 4, 4]

            itestkey = 'BF01_%i' % wave
            diffBF01w['test'] = itestkey
            isurrogate = 'PTC02_%i' % wave

            bf01w = BF01.BF01(inputs=dict(elvis=elvis,
                                          CHAMBER=CHAMBER,
                                          test=itestkey,
                                          exptimes=exptsBF01w,
                                          frames=frsBF01w,
                                          wavelength=wave,
                                          Npix=5,
                                          clipsigma=4.,
                                          covfunc='ver2',
                                          doBiasCorr=True,
                                          central='mean',
                                          surrogate=isurrogate,
                                          diffvalues=diffBF01w))

            test_sequence[itestkey] = copy.deepcopy(bf01w)

    if _toGen['FLATFLUX00']:

        wavesFLATFLUX00 = [590, 730, 800, 880, 0]

        diffFLATFLUX00w = dict()
        diffFLATFLUX00w.update(diffvalues)

        for iw, wave in enumerate(wavesFLATFLUX00):

            tFWC_flatw = ogse.profile['tFWC_point']['nm%i' % wave]
            exptsFLATFLUX00w = (np.array([0.1, 0.3, 0.6, 0.8]) * tFWC_flatw).tolist()
            frsFLATFLUX00 = [1, 1, 1, 1]

            itestkey = 'FLATFLUX00_%i' % wave
            diffFLATFLUX00w['test'] = itestkey

            flatflux00w = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                                  CHAMBER=CHAMBER,
                                                  wavelength=wave,
                                                  exptimes=exptsFLATFLUX00w,
                                                  frames=frsFLATFLUX00,
                                                  test=itestkey,
                                                  diffvalues=diffFLATFLUX00w))

            test_sequence[itestkey] = copy.deepcopy(flatflux00w)

    # PTC-02 - Temp.

    if _toGen['PTC02TEMP']:

        wavePTC02T = 800
        TempsPTC02T = [148., 158.]

        diffPTC02T = dict(mirr_on=0)
        diffPTC02T.update(diffvalues)

        # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
        tFWC_flatw = ogse.profile['tFWC_flat']['nm%i' % wavePTC02T]
        exptsPTC02T = (np.array([10., 30., 50., 70., 80., 90.]) /
                       100. * tFWC_flatw).tolist()
        frsPTC02T = [4, 4, 4, 4, 4, 4]

        for it, T in enumerate(TempsPTC02T):

            itestkey = 'PTC02_%iK' % T

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

    # PTC-02 - RD.

    if _toGen['PTC02RD']:

        wavePTC02RD = 800
        PTC02RDs = [15.5, 16.0, 16.5]

        # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
        tFWC_flatw = ogse.profile['tFWC_flat']['nm%i' % wavePTC02RD]
        exptsPTC02RD = (np.array([10., 30., 50., 70., 80., 90.]) /
                        100. * tFWC_flatw).tolist()
        frsPTC02RD = [4, 4, 4, 4, 4, 4]

        for ir, RD in enumerate(PTC02RDs):

            diffPTC02RD = dict(mirr_on=0, RD_T=RD,
                               RD_B=RD)
            diffPTC02RD.update(diffvalues)

            itestkey = 'PTC02_%iR' % (RD * 10.,)

            diffPTC02RD['test'] = itestkey

            ptc02rd = PTC0X.PTC0X(inputs=dict(elvis=elvis,
                                              CHAMBER=CHAMBER,
                                              test=itestkey,
                                              wavelength=wavePTC02RD,
                                              frames=frsPTC02RD,
                                              exptimes=exptsPTC02RD,
                                              diffvalues=diffPTC02RD))

            # istructPTC02T = ptc02t.build_scriptdict(diffvalues=diffPTC02T,
            #                                        elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(ptc02rd)

    # NL-01

    if _toGen['NL01']:

        waveNL01 = 0

        diffNL01 = dict(mirr_on=0)
        diffNL01.update(diffvalues)

        # 5 frames per fluence: 1%, 2%, 3%, 5%, 10%, 20%,30%, 50%,70%,80%,85%,90%,95%
        tFWC_flatNL01 = ogse.profile['tFWC_flat']['nm%i' % waveNL01]
        exptsNL01 = (np.array([0.5, 0.7, 1., 2., 3., 5., 10., 20., 30., 55., 70., 80., 85.,
                               90., 95., 100., 110.]) / 100. * tFWC_flatNL01).tolist()  # ms
        exptinterNL01 = 0.5 * tFWC_flatNL01
        frsNL01 = (np.ones(len(exptsNL01), dtype='int32') * 4).tolist()

        nl01 = NL01.NL01(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='NL01',
                                     exptimes=exptsNL01,
                                     exptinter=exptinterNL01,
                                     frames=frsNL01,
                                     wavelength=waveNL01,
                                     diffvalues=diffNL01))

        #structNL01 = nl01.build_scriptdict(diffvalues=diffNL01,elvis=elvis)
        test_sequence['NL01'] = copy.deepcopy(nl01)

    # NL-02

    if _toGen['NL02']:

        diffNL02 = dict(mirr_on=0)
        diffNL02.update(diffvalues)

        FLUDIVIDE = 20.

        relfluencesNL02 = np.array([0.5, 0.7, 1., 2., 3., 5., 10., 15.,
                                    20., 25., 50., 70., 80., 85., 90., 95., 100., 110.])

        ixFLUDIVIDE = np.where(relfluencesNL02 == FLUDIVIDE)[0][0]

        waveNL02A = 0
        tFWCwA = ogse.profile['tFWC_flat']['nm%i' % waveNL02A]

        ixLOWFLU = (np.arange(ixFLUDIVIDE + 1 + 1),)
        exptsNL02A = (relfluencesNL02[ixLOWFLU] / 100. * tFWCwA).tolist()  # ms
        framesNL02A = (np.ones(len(exptsNL02A), dtype='int32') * 4).tolist()

        waveNL02B = 880
        tFWCwB = ogse.profile['tFWC_flat']['nm%i' % waveNL02B]

        ixHIFLU = (np.arange(ixFLUDIVIDE - 1, len(relfluencesNL02)),)
        exptsNL02B = (relfluencesNL02[ixHIFLU] / 100. * tFWCwB).tolist()  # ms
        framesNL02B = (np.ones(len(exptsNL02B), dtype='int32') * 4).tolist()

        nl02 = NL02.NL02(inputs=dict(elvis=elvis,
                                     CHAMBER=CHAMBER,
                                     test='NL02',
                                     wavelengthA=waveNL02A,
                                     exptimesA=exptsNL02A,
                                     framesA=framesNL02A,
                                     wavelengthB=waveNL02B,
                                     exptimesB=exptsNL02B,
                                     framesB=framesNL02B,
                                     exptinter=0.5 * tFWCwB,
                                     diffvalues=diffNL02))

        test_sequence['NL02'] = copy.deepcopy(nl02)

    # NL-02

    if _toGen['NL02RD']:

        RDs = [15.5, 16.0, 16.5]

        FLUDIVIDE = 20.

        relfluencesNL02RD = np.array([0.5, 0.7, 1., 2., 3., 5., 10., 20., 30., 55., 70., 80., 85.,
                                      90., 95., 100., 110.])

        waveNL02A = 0
        tFWCwA = ogse.profile['tFWC_flat']['nm%i' % waveNL02A]

        ixLOWFLU = np.where(relfluencesNL02RD < FLUDIVIDE)
        exptsNL02A = (relfluencesNL02RD[ixLOWFLU] / 100. * tFWCwA).tolist()  # ms
        framesNL02A = (np.ones(len(exptsNL02A), dtype='int32') * 4).tolist()

        waveNL02B = 880
        tFWCwB = ogse.profile['tFWC_flat']['nm%i' % waveNL02B]

        ixHIFLU = np.where(relfluencesNL02RD >= FLUDIVIDE)
        exptsNL02B = (relfluencesNL02RD[ixHIFLU] / 100. * tFWCwB).tolist()  # ms
        framesNL02B = (np.ones(len(exptsNL02B), dtype='int32') * 4).tolist()

        for ir, RD in enumerate(RDs):

            diffNL02RD = dict(mirr_on=0, RD_T=RD,
                              RD_B=RD)
            diffNL02RD.update(diffvalues)

            itestkey = 'NL02_%iR' % (RD * 10.,)

            diffNL02RD['test'] = itestkey

            nl02rd = NL02.NL02(inputs=dict(elvis=elvis,
                                           CHAMBER=CHAMBER,
                                           test=itestkey,
                                           wavelengthA=waveNL02A,
                                           exptimesA=exptsNL02A,
                                           framesA=framesNL02A,
                                           wavelengthB=waveNL02B,
                                           exptimesB=exptsNL02B,
                                           framesB=framesNL02B,
                                           exptinter=0.5 * tFWCwB,
                                           diffvalues=diffNL02RD))

            test_sequence[itestkey] = copy.deepcopy(nl02rd)

    # FOCUS

    if _toGen['FOCUS00']:

        # wavesFOCUS00w = [590,640,730,800,880] # TESTS
        wavesFOCUS00w = [590, 730, 800, 880]

        diffFOCUS00w = dict()
        diffFOCUS00w.update(diffvalues)

        for iw, wave in enumerate(wavesFOCUS00w):

            tFWC_pointw = ogse.profile['tFWC_point']['nm%i' % wave]

            iexptimeF00 = 60. / 100. * tFWC_pointw

            deltafocusF00 = 0.2  # mm

            itestkey = 'FOCUS00_%i' % wave

            diffFOCUS00w['test'] = itestkey

            focus00 = FOCUS00.FOCUS00(inputs=dict(elvis=elvis,
                                                  CHAMBER=CHAMBER,
                                                  test=itestkey,
                                                  exptime=iexptimeF00,
                                                  wavelength=wave,
                                                  deltafocus=deltafocusF00,
                                                  diffvalues=diffFOCUS00w))
            # istructFOCUS00w = focus00.build_scriptdict(diffvalues=diffFOCUS00w,
            #                                                   elvis=elvis)
            test_sequence[itestkey] = copy.deepcopy(focus00)

    # PSF

    if _toGen['PSF01']:

        #
        PSF0Xtestdefaults = PSF0X.get_testdefaults(ogse)
        #wavesPSF01w = PSF0Xtestdefaults['waves']
        wavesPSF01w = [590, 730, 800, 880]

        diffPSF01w = dict()
        diffPSF01w.update(diffvalues)

        for iw, wave in enumerate(wavesPSF01w):

            exptsPSF01w = PSF0Xtestdefaults['exptimes']['nm%i' % wave]
            frsPSF01w = PSF0Xtestdefaults['frames']

            itestkey = 'PSF01_%i' % wave
            diffPSF01w['test'] = itestkey

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

    if _toGen['PSFLUX00']:

        wavesPSFLUX00 = [590, 730, 800, 880, 0]

        diffPSFLUX00w = dict()
        diffPSFLUX00w.update(diffvalues)

        for iw, wave in enumerate(wavesPSFLUX00):

            tFWC_pointw = ogse.profile['tFWC_point']['nm%i' % wave]
            exptsPSFLUX00w = (
                np.array([0.1, 0.3, 0.6, 0.8]) * tFWC_pointw).tolist()
            frsPSFLUX00 = [1, 1, 1, 1]

            itestkey = 'PSFLUX00_%i' % wave
            diffPSFLUX00w['test'] = itestkey

            psflux00w = PSF0X.PSF0X(inputs=dict(elvis=elvis,
                                                CHAMBER=CHAMBER,
                                                wavelength=wave,
                                                exptimes=exptsPSFLUX00w,
                                                frames=frsPSFLUX00,
                                                test=itestkey,
                                                diffvalues=diffPSFLUX00w))

            test_sequence[itestkey] = copy.deepcopy(psflux00w)

    # PSF02

    if _toGen['PSF02']:

        wPSF02 = 800
        temps_PSF02 = [148, 158]

        tFWC_pointPSF02 = ogse.profile['tFWC_point']['nm%i' % wPSF02]

        exptsPSF02 = np.array([5., 25., 50., 75., 90.]) / \
            100. * tFWC_pointPSF02
        frsPSF02 = [20, 14, 10, 4, 4]

        diffPSF02 = dict()
        diffPSF02.update(diffvalues)

        for it, temp in enumerate(temps_PSF02):

            itestkey = 'PSF02_%iK' % temp

            diffPSF02['test'] = itestkey

            psf02k = PSF0X.PSF0X(inputs=dict(elvis=elvis,
                                             CHAMBER=CHAMBER,
                                             wavelength=wPSF02,
                                             exptimes=exptsPSF02.tolist(),
                                             frames=frsPSF02,
                                             test=itestkey,
                                             diffvalues=diffPSF02))

            # istructPSF02 = psf02k.build_scriptdict(diffvalues=diffPSF02.copy(),
            #                                       elvis=elvis)

            test_sequence[itestkey] = copy.deepcopy(psf02k)

    # OTHER

    # PERSIST

    if _toGen['PERSIST01']:

        #wavePERS = 2000
        #tFWC_point_PERS = ogse.profile['tFWC_point']['nm%i' % wavePERS]

        # exptPER01_SATUR = tFWC_point_PERS*2500.   # s
        exptPER01_SATUR = 600.  # s., HARDWIRED
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

    if _toGen['MOT_WARM']:

        diffMOT_WM = dict()
        diffMOT_WM.update(diffvalues)

        mot_wm = MOT_WARM.MOT_WARM(inputs=dict(elvis=elvis,
                                               CHAMBER=CHAMBER,
                                               test='MOT_WARM',
                                               diffvalues=diffMOT_WM))

        test_sequence['MOT_WARM'] = copy.deepcopy(mot_wm)

    if _toGen['COSMETICS00']:

        diffCOS = dict()
        diffCOS.update(diffvalues)

        cos = COSMETICS00.COSMETICS00(inputs=dict(elvis=elvis,
                                                  CHAMBER=CHAMBER,
                                                  test='COSMETICS00',
                                                  diffvalues=diffCOS))

        test_sequence['COSMETICS00'] = copy.deepcopy(cos)

    return test_sequence
