#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

VIS Ground Calibration Campaign.

Mini-Campaign Test Set: used to verify functionality of ROE + OGSE 
ahead of running the campaign itself.

Created on Mon Jan 15 15:56:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from collections import OrderedDict
from pdb import set_trace as stop
import numpy as np

from vison.point import FOCUS00, PSF0X
from vison.dark import BIAS01, DARK01
from vison.flat import FLAT0X
from vison.inject import CHINJ00
from vison.pump import TP00
#from vison.other import PERSIST01 as PER01
#from vison.point import lib as polib
from vison.ogse import ogse
#from vison.pipe import lib as pilib
from vison.support import context
# END IMPORT


def generate_reduced_test_sequence(equipment, toGen, elvis=context.elvis):
    """ """

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

        bias01 = BIAS01.BIAS01(inputs=dict(N=Nbias01, diffvalues=diffBIAS01))

        structBIAS01 = bias01.build_scriptdict(elvis=elvis)

        test_sequence['BIAS01'] = structBIAS01

    # DARKS

    if toGen['DARK01']:

        print 'DARK01...'

        Ndark01 = 1
        exptime_dark01 = 565.  # s
        diffDARK01 = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        dark01 = DARK01.DARK01(inputs=dict(N=Ndark01, exptime=exptime_dark01,
                                           diffvalues=diffDARK01))

        structDARK01 = dark01.build_scriptdict(elvis=elvis)

        test_sequence['DARK01'] = structDARK01

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
                                              diffvalues=diffCHINJ00))
        structCHINJ00 = chinj00.build_scriptdict(elvis=elvis)

        test_sequence['CHINJ00'] = structCHINJ00


# ==============================================================================
#     if toGen['CHINJ01']:
#
#         print 'CHINJ01...'
#
#         # CHINJ01
#
#         IDL = 11.
#         IDH = 18.
#         IG1s = [2.,6.]
#         toi_chinj01 = 500
#         id_delays = [toi_chinj01*3,toi_chinj01*2]
#
#         diffCHINJ01 = dict(mirr_on=0,sn_ccd1=sn_ccd1,
#                       sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                       sn_rpsu=sn_rpsu,operator=operator)
#
#         chinj01 = CHINJ01.CHINJ01(inputs=dict(IDL=IDL,IDH=IDH,IG1s=IG1s,
#                                               toi_chinj=toi_chinj01,
#                                               id_delays=id_delays,
#                                               diffvalues=diffCHINJ01))
#         structCHINJ01 = chinj01.build_scriptdict(elvis=elvis)
#
#         test_sequence['CHINJ01'] = structCHINJ01
#
#
#
#     # CHINJ02
#
#     if toGen['CHINJ02']:
#
#         print 'CHINJ02...'
#
#         IDLs = [13.,16.]
#         IDH = 18.
#         #IG1s = [2000,6000]
#         toi_chinj02 = 500
#         id_delays = [toi_chinj02*3,toi_chinj02*2]
#         diffCHINJ02 = dict(pos_cal_mirror=polib.mirror_nom['F4'],sn_ccd1=sn_ccd1,
#                           sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                           sn_rpsu=sn_rpsu,operator=operator)
#
#         chinj02 = CHINJ02.CHINJ02(inputs=dict(IDLs=IDLs,IDH=IDH,toi_chinj=toi_chinj02,
#                                               id_delays=id_delays,diffvalues=diffCHINJ02))
#         structCHINJ02 = chinj02.build_scriptdict(elvis=elvis)
#
#         test_sequence['CHINJ02'] = structCHINJ02
#
# ==============================================================================

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
                                     diffvalues=diffTP00))

        structTP00 = tp00.build_scriptdict(elvis=elvis)

        test_sequence['TP00'] = structTP00


# ==============================================================================
#     # TP02
#
#     if toGen['TP02']:
#
#         print 'TP02...'
#
#         Nshuffles_H=5000
#         dwell_sv = [0.,4.75,14.3,28.6] # us
#         toi_chinjTP02 = 250 # quick injection
#         id_delays_TP02 = np.array([3.,2.])*toi_chinjTP02
#
#         diffTP02 = dict(sn_ccd1=sn_ccd1,
#                           sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                           sn_rpsu=sn_rpsu,operator=operator,toi_chinj=toi_chinjTP02)
#
#         tp02 = TP02.TP02(inputs=dict(Nshuffles_H=Nshuffles_H,
#                                      dwell_sv=dwell_sv,toi_chinj=toi_chinjTP02,
#                                      id_delays=id_delays_TP02,
#                                      diffvalues=diffTP02))
#
#         structTP02 = tp02.build_scriptdict(elvis=elvis)
#
#         test_sequence['TP02'] = structTP02
#
#
# ==============================================================================
    # FLATS

    exptimes_FLAT0X = dict(nm0=ogse.tFWC_flat['nm0'],
                           nm590=ogse.tFWC_flat['nm590'],
                           nm640=ogse.tFWC_flat['nm640'],
                           nm730=ogse.tFWC_flat['nm730'],
                           nm800=ogse.tFWC_flat['nm800'],
                           nm880=ogse.tFWC_flat['nm880'])
    # FLAT-01

    if toGen['FLAT01']:

        print 'FLAT01...'

        t_dummy_F01 = np.array([25., 50., 75])/100.
        exptimesF01 = (exptimes_FLAT0X['nm800'] * t_dummy_F01).tolist()  # s
        framesF01 = [1, 1, 1]

        inpF01 = dict(exptimes=exptimesF01,
                      frames=framesF01,
                      wavelength=800,
                      test='FLAT01_800')

        diffFLAT01 = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        flat01 = FLAT0X.FLAT0X(inputs=inpF01)

        structFLAT01 = flat01.build_scriptdict(
            diffvalues=diffFLAT01, elvis=elvis)
        test_sequence['FLAT01'] = structFLAT01

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
                          test=itestkey)

            flat02 = FLAT0X.FLAT0X(inputs=inpF02)

            istructFLAT02 = flat02.build_scriptdict(
                diffvalues=diffFLAT02, elvis=elvis)

            test_sequence[itestkey] = istructFLAT02


# ==============================================================================
#     # PTC
#
#     # PTC-01
#
#     if toGen['PTC01']:
#
#         print 'PTC01...'
#
#         diffPTC01 = dict(test='PTC01',
#                          vstart=0,
#                          vend=2086,
#                          sn_ccd1=sn_ccd1,
#                          sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                          sn_rpsu=sn_rpsu,operator=operator)
#
#         # 5%, 10%, 20%, 30%, 50%, 70%, 80%, 90%, 100%, 110%, 120%
#         exptsPTC01 = np.array([5.,30.,50.,80.,100.,120.])/100.*ogse.tFWC_flat['nm800'] # ms
#         frsPTC01 = [2,2,2,2,2,2]
#
#         ptc01 = PTC0X.PTC0X(inputs=dict(test='PTC01',exptimes=exptsPTC01,
#                                         frames=frsPTC01,wavelength=800))
#         structPTC01 = ptc01.build_scriptdict(diffvalues=diffPTC01,elvis=elvis)
#
#         test_sequence['PTC01'] = structPTC01
#
#
# ==============================================================================
# ==============================================================================
#     # PTC-02 - wavelength
#
#
#     if toGen['PTC02WAVE']:
#
#         print 'PTC02WAVE...'
#
#         wavesPTC02w = [590,640,730,880]
#
#         diffPTC02w = dict(sn_ccd1=sn_ccd1,
#                       sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                       sn_rpsu=sn_rpsu,operator=operator)
#
#         for iw, wave in enumerate(wavesPTC02w):
#
#             # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
#
#             exptsPTC02w = np.array([10.,30.,50.,70.,80.,90.])/100.*ogse.tFWC_flat['nm%i' % wave]
#             frsPTC02w = [4,4,4,4,4,4]
#
#             ptc02w = PTC0X.PTC0X(inputs=dict(test='PTC02_%i' % wave,
#                                              wavelength=wave,
#                                              exptimes=exptsPTC02w,
#                                              frames=frsPTC02w))
#
#             itestkey = 'PTC02_%i' % wave
#             diffPTC02w['test'] = itestkey
#
#             print '%s...' % itestkey
#
#             istructPTC02w = ptc02w.build_scriptdict(diffvalues=diffPTC02w,elvis=elvis)
#             test_sequence[itestkey] = istructPTC02w
#
#
#    # PTC-02 - Temp.
#
#
#     if toGen['PTC02TEMP']:
#
#         print 'PTC02TEMP...'
#
#         wavePTC02T = 800
#         TempsPTC02T = [150.,156.]
#
#         diffPTC02T = dict(sn_ccd1=sn_ccd1,
#                           sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                           sn_rpsu=sn_rpsu,operator=operator)
#
#         # 10%, 30%, 50%, 70%, 80%, 90% x FWC. 4 frames per fluence.
#         exptsPTC02T = np.array([10.,30.,50.,70.,80.,90.])/100.*ogse.tFWC_flat['nm%i' % wavePTC02T]
#         frsPTC02T = [4,4,4,4,4,4]
#
#         for it,T in enumerate(TempsPTC02T):
#
#             itestkey = 'PTC02_%iK' % T
#             print '%s...' % itestkey
#
#             diffPTC02T['test'] = itestkey
#
#             ptc02t = PTC0X.PTC0X(inputs=dict(test=itestkey,
#                                              wavelength=wavePTC02T,
#                                              frames=frsPTC02T,
#                                              exptimes=exptsPTC02T))
#             istructPTC02T = ptc02t.build_scriptdict(diffvalues=diffPTC02T,
#                                                     elvis=elvis)
#
#             test_sequence[itestkey] = istructPTC02T
#
#
#     # NL-01
#
#
#     if toGen['NL01']:
#
#         print 'NL01...'
#
#         diffNL01 = dict(test='NL01',sn_ccd1=sn_ccd1,
#                           sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                           sn_rpsu=sn_rpsu,operator=operator)
#
#         # 5 frames per fluence: 1%, 2%, 3%, 5%, 10%, 20%,30%, 50%,70%,80%,85%,90%,95%
#         exptsNL01 = np.array([5.,10.,20.,30.,50.,70.,80.,90.,100.,110.,120.])/100. * ogse.tFWC_flat['nm0'] # ms
#         exptinterNL01 = 0.5 * ogse.tFWC_flat['nm0']
#         frsNL01 = np.ones(11,dtype='int32')*5
#
#         nl01 = NL01.NL01(inputs=dict(exptimes=exptsNL01,exptinter=exptinterNL01,
#                                      frames=frsNL01,wavelength=0))
#         structNL01 = nl01.build_scriptdict(diffvalues=diffNL01,elvis=elvis)
#
#         test_sequence['NL01'] = structNL01
#
# ==============================================================================
    # FOCUS

    if toGen['FOCUS00']:

        print 'FOCUS00...'

        wavesFOCUS00w = [590, 800, 880]

        diffFOCUS00w = dict(sn_ccd1=sn_ccd1,
                            sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                            sn_rpsu=sn_rpsu, operator=operator)

        for iw, wave in enumerate(wavesFOCUS00w):

            iexptimeF00 = 60./100. * ogse.tFWC_point['nm%i' % wave]

            itestkey = 'FOCUS00_%i' % wave

            diffFOCUS00w['test'] = itestkey

            print '%s...' % itestkey

            focus00 = FOCUS00.FOCUS00(inputs=dict(wavelength=wave,
                                                  exptime=iexptimeF00))
            istructFOCUS00w = focus00.build_scriptdict(diffvalues=diffFOCUS00w,
                                                       elvis=elvis)
            test_sequence[itestkey] = istructFOCUS00w

    # PSF

    if toGen['PSF01']:

        print 'PSF01...'

        wavePSF01w = 800

        diffPSF01w = dict(sn_ccd1=sn_ccd1,
                          sn_ccd2=sn_ccd2, sn_ccd3=sn_ccd3, sn_roe=sn_roe,
                          sn_rpsu=sn_rpsu, operator=operator)

        exptsPSF01w = (np.array([25., 50., 75.])/100. *
                       ogse.tFWC_point['nm%i' % wavePSF01w]).tolist()
        frsPSF01w = [4, 4, 4]

        itestkey = 'PSF01_%i' % wavePSF01w
        diffPSF01w['test'] = itestkey

        print '%s...' % itestkey

        psf01w = PSF0X.PSF0X(inputs=dict(wavelength=wave,
                                         exptimes=exptsPSF01w,
                                         frames=frsPSF01w,
                                         test=itestkey))

        istructPSF01w = psf01w.build_scriptdict(diffvalues=diffPSF01w,
                                                elvis=elvis)

        test_sequence[itestkey] = istructPSF01w


# ==============================================================================
#     # PSF02
#
#
#     if toGen['PSF02']:
#
#
#         print 'PSF02...'
#
#         wPSF02 = 800
#         temps_PSF02 = [150,156]
#
#         exptsPSF02 = np.array([5.,25.,50.,75.,90.])/100. * ogse.tFWC_point['nm%i' % wPSF02]
#         frsPSF02 = [20,15,10,4,3]
#
#         diffPSF02 = dict(sn_ccd1=sn_ccd1,
#                           sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                           sn_rpsu=sn_rpsu,operator=operator)
#
#
#         for it,temp in enumerate(temps_PSF02):
#
#
#             itestkey = 'PSF02_%iK' % temp
#
#             diffPSF02['test'] = itestkey
#
#             print '%s...' % itestkey
#
#             psf02k = PSF0X.PSF0X(inputs=dict(wavelength=800,
#                                              exptimes=exptsPSF02,
#                                              frames=frsPSF02,
#                                              test=itestkey))
#
#             istructPSF02 = psf02k.build_scriptdict(diffvalues=diffPSF02,
#                                                    elvis=elvis)
#
#             test_sequence[itestkey] = istructPSF02
#
#     # OTHER
#
#     # PERSIST
#
#
#     if toGen['PERSIST01']:
#
#         print 'PERSIST01...'
#
#         exptPER01_SATUR = 15.   # s
#         exptPER01_LATEN = 565. # s
#
#         diffPER01 = dict(sn_ccd1=sn_ccd1,
#                           sn_ccd2=sn_ccd2,sn_ccd3=sn_ccd3,sn_roe=sn_roe,
#                           sn_rpsu=sn_rpsu,operator=operator)
#
#         persist01 = PER01.PERSIST01(inputs=dict(exptSATUR=exptPER01_SATUR,
#                                                 exptLATEN=exptPER01_LATEN))
#
#         structPER01 = persist01.build_scriptdict(diffvalues=diffPER01,elvis=elvis)
#
#         test_sequence['PERSIST01'] = structPER01
#
# ==============================================================================

    return test_sequence