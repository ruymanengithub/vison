#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in NL analysis.

Created on Mon Feb 5 15:51:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import datetime
from collections import OrderedDict

from scipy import interpolate
# END IMPORT

FullDynRange = 2.**16
NLdeg = 7


def get_exptime_atmiddynrange(flu1D, exp1D, method='spline', debug=False):
    """ """
    
    X = flu1D/FullDynRange
    Y = exp1D.copy()
    ixsort = np.argsort(X)
    X = X[ixsort].copy()
    Y = Y[ixsort].copy()
    
    splinefunc = interpolate.interp1d(X,Y,kind='linear')
    
    mod1d_fit = np.polyfit(X, Y, 2)  # a linear approx. is fine
    mod1d_pol = np.poly1d(mod1d_fit)
    t50poly = mod1d_pol(0.5)
    t50spline = splinefunc(0.5)
    
    if debug:
        try:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(X,Y,'b.-')
            ax.axvline(0.5,c='k',ls='--')
            ax.axhline(t50poly,c='r',ls='--')
            ax.axhline(t50spline,c='b',ls='-')
            plt.show()
        except:
            stop()
    
    if method == 'spline':
        return t50spline
    elif method == 'poly':
        return t50poly
    


def fitNL(fluencesNL, exptimes, display=False):
    """ """

    assert fluencesNL.shape[0] == exptimes.shape[0]
    assert fluencesNL.ndim <= 2
    assert exptimes.ndim == 1

    nomG = 3.5 # e/ADU
    minfitFl = 2000. # ADU
    maxfitFl = FullDynRange-10000. # ADU

    Nexp = len(exptimes)

    if fluencesNL.ndim == 2:
        
        Nsec = fluencesNL.shape[1]
        #_exptimes = np.repeat(exptimes.reshape(Nexp,1),Nsec,axis=1)

        t50 = np.zeros(Nsec, dtype='float32') + np.nan

        for i in range(Nsec):
            t50[i] = get_exptime_atmiddynrange(fluencesNL[:, i], exptimes, 
               method='spline', debug=False)

        t50 = np.repeat(t50.reshape(1, Nsec), Nexp, axis=0)

        exptimes_bc = np.repeat(exptimes.reshape(Nexp, 1), Nsec, axis=1)

    else:

        t50 = get_exptime_atmiddynrange(fluencesNL[:, i], exptimes,
                                        method='spline')

        exptimes_bc = exptimes.copy()

    YL = exptimes_bc/t50 * FullDynRange/2.
    Z = 100.*(fluencesNL/YL-1.) 
    
    efNL = np.sqrt(fluencesNL*nomG)/nomG
    
    W = 100.*(efNL/YL)
    
    X = fluencesNL.flatten().copy()
    Y = Z.flatten().copy()
    ixsort = np.argsort(X)
    X = X[ixsort].copy()
    Y = Y[ixsort].copy()
    W = W.flatten()[ixsort].copy()
    
    selix = np.where((X>minfitFl) & (X<maxfitFl))

    NLfit = np.polyfit(X[selix], Y[selix], w=W[selix], deg=NLdeg, full=False)

    NLpol = np.poly1d(NLfit)

    fkfluencesNL = np.arange(minfitFl, maxfitFl, dtype='float32')
    # array with NL fluences (1 to 2**16, in steps of ADU)
    Z_bestfit = NLpol(fkfluencesNL)
    
    ixmax = np.abs(Z_bestfit).argmax()
    maxNLpc = Z_bestfit[ixmax]
    flu_maxNLpc = fkfluencesNL[ixmax]
    
    if display:
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(X,Y,'k.')
        ax.plot(X[selix],Y[selix],'b.')
        ax.plot(fkfluencesNL,Z_bestfit,'r--')
        #ax.set_ylim([-10.,10.])
        plt.show()

    fitresults = OrderedDict(
                             coeffs=NLfit, 
                             NLdeg=NLdeg, 
                             maxNLpc=maxNLpc,
                             flu_maxNLpc=flu_maxNLpc,
                             inputcurve = OrderedDict(
                                     YL=fluencesNL.flatten().copy(),
                                     Z=Z.flatten().copy()),
                             outputcurve = OrderedDict(
                                     YL=fkfluencesNL.copy(),
                                     Z=Z_bestfit.copy()))
                             
    
    
    return fitresults


def wrap_fitNL(fluences, variances, exptimes, col_labels, times=np.array([]), TrackFlux=True, subBgd=True):
    """ """
    # col001 == BGD
    # colEVEN = STAB
    # colODD = Fluences != 0

    NObsIDs, Nsecs = fluences.shape

    dtimes = np.array(
        [(times[i]-times[0]).seconds for i in range(NObsIDs)], dtype='float32')

    col_numbers = np.array([int(item[3:]) for item in col_labels])

    ixboo_bgd = col_numbers == 1
    ixboo_stab = (col_numbers % 2 == 0) & (col_numbers > 1)  # EVEN
    ixboo_fluences = (col_numbers % 2 != 0)  & (col_numbers > 1)# ODD

    if subBgd:
        bgd = np.nanmean(fluences[ixboo_bgd, :], axis=0)
        bgd = np.repeat(bgd.reshape(1, Nsecs), NObsIDs, axis=0).copy()
        fluences -= bgd
    else:
        bgd = 0.

    if TrackFlux and len(times) == NObsIDs:

        st_dtimes = dtimes[ixboo_stab].copy()
        st_fluences = np.nanmean(fluences[ixboo_stab, :], axis=1).copy()

        track_fit = np.polyfit(st_dtimes, st_fluences,
                               2, full=False, cov=False)
        track_pol = np.poly1d(track_fit)

        track = track_pol(dtimes[ixboo_fluences])
        track /= np.median(track)
        #track_bc = np.repeat(track.reshape(len(track), 1), Nsecs, axis=1)
        fluences[ixboo_fluences, :] /= track.reshape(track.shape[0],-1)
    else:
        track = np.ones_like(fluences[ixboo_fluences,0])
    
    fitresults = fitNL(fluences[ixboo_fluences, :], exptimes[ixboo_fluences],
                       display=False)
    fitresults['bgd'] = bgd
    fitresults['stability_pc'] = np.std(track)*100.

    return fitresults


def test_wrap_fitNL():
    """ """

    col_labels = np.array(['col001', 'col001', 'col001', 'col001',
                           'col002', 'col002', 'col002', 'col002', 'col002',
                           'col003',
                           'col004', 'col004', 'col004', 'col004', 'col004',
                           'col005',
                           'col006', 'col006', 'col006', 'col006', 'col006',
                           'col007',
                           'col008', 'col008', 'col008', 'col008', 'col008',
                           'col009',
                           'col010', 'col010', 'col010', 'col010', 'col010',
                           'col013',
                           'col014', 'col014', 'col014', 'col014', 'col014',
                           'col015',
                           'col016', 'col016', 'col016', 'col016', 'col016',
                           'col017',
                           'col018', 'col018', 'col018', 'col018', 'col018',
                           'col019',
                           'col020', 'col020', 'col020', 'col020', 'col020',
                           'col021',
                           'col022', 'col022', 'col022', 'col022', 'col022',
                           'col023',
                           'col024', 'col024', 'col024', 'col024', 'col024',
                           'col025'],
                          dtype='|S12')

    exptimes = np.array([0.,   0.,   0.,   0.,   1.5,   1.5,   1.5,   1.5,   1.5,
                         15.,   3.,   3.,   3.,   3.,   3.,  15.,   6.,   6.,
                         6.,   6.,   6.,  15.,   9.,   9.,   9.,   9.,   9.,
                         15.,  15.,  15.,  15.,  15.,  15.,  15.,  21.,  21.,
                         21.,  21.,  21.,  15.,  24.,  24.,  24.,  24.,  24.,
                         15.,  27.,  27.,  27.,  27.,  27.,  15.,  30.,  30.,
                         30.,  30.,  30.,  15.,  33.,  33.,  33.,  33.,  33.,
                         15.,  36.,  36.,  36.,  36.,  36.,  15.])

    raw_data = np.repeat((exptimes/exptimes.max()*2. **
                          16).reshape((70, 1)), 49, axis=1)
    dtobjs = np.array([datetime.datetime(2018, 2, 8, 3, 59, 17),
                       datetime.datetime(2018, 2, 8, 4, 0, 55), datetime.datetime(
                           2018, 2, 8, 4, 2, 28),
                       datetime.datetime(2018, 2, 8, 4, 4, 7), datetime.datetime(
                           2018, 2, 8, 4, 5, 43),
                       datetime.datetime(2018, 2, 8, 4, 7, 19), datetime.datetime(
                           2018, 2, 8, 4, 8, 54),
                       datetime.datetime(2018, 2, 8, 4, 10, 29), datetime.datetime(
                           2018, 2, 8, 4, 12, 5),
                       datetime.datetime(2018, 2, 8, 4, 14, 2), datetime.datetime(
                           2018, 2, 8, 4, 15, 40),
                       datetime.datetime(2018, 2, 8, 4, 17, 18), datetime.datetime(
                           2018, 2, 8, 4, 18, 55),
                       datetime.datetime(2018, 2, 8, 4, 20, 38), datetime.datetime(
                           2018, 2, 8, 4, 22, 15),
                       datetime.datetime(2018, 2, 8, 4, 24, 5), datetime.datetime(
                           2018, 2, 8, 4, 25, 46),
                       datetime.datetime(2018, 2, 8, 4, 27, 27), datetime.datetime(
                           2018, 2, 8, 4, 29, 6),
                       datetime.datetime(2018, 2, 8, 4, 30, 47), datetime.datetime(
                           2018, 2, 8, 4, 32, 28),
                       datetime.datetime(2018, 2, 8, 4, 34, 18), datetime.datetime(
                           2018, 2, 8, 4, 36, 2),
                       datetime.datetime(2018, 2, 8, 4, 37, 46), datetime.datetime(
                           2018, 2, 8, 4, 39, 28),
                       datetime.datetime(2018, 2, 8, 4, 41, 19), datetime.datetime(
                           2018, 2, 8, 4, 43, 11),
                       datetime.datetime(2018, 2, 8, 4, 45), datetime.datetime(
                           2018, 2, 8, 4, 46, 57),
                       datetime.datetime(2018, 2, 8, 4, 48, 46), datetime.datetime(
                           2018, 2, 8, 4, 50, 34),
                       datetime.datetime(2018, 2, 8, 4, 52, 24), datetime.datetime(
                           2018, 2, 8, 4, 54, 13),
                       datetime.datetime(2018, 2, 8, 4, 56, 10), datetime.datetime(
                           2018, 2, 8, 4, 58, 6),
                       datetime.datetime(2018, 2, 8, 5, 0, 5), datetime.datetime(
                           2018, 2, 8, 5, 2),
                       datetime.datetime(2018, 2, 8, 5, 4, 2), datetime.datetime(
                           2018, 2, 8, 5, 5, 58),
                       datetime.datetime(2018, 2, 8, 5, 7, 48), datetime.datetime(
                           2018, 2, 8, 5, 9, 48),
                       datetime.datetime(2018, 2, 8, 5, 11, 46), datetime.datetime(
                           2018, 2, 8, 5, 13, 43),
                       datetime.datetime(2018, 2, 8, 5, 15, 42), datetime.datetime(
                           2018, 2, 8, 5, 17, 40),
                       datetime.datetime(2018, 2, 8, 5, 19, 30), datetime.datetime(
                           2018, 2, 8, 5, 21, 40),
                       datetime.datetime(2018, 2, 8, 5, 23, 46), datetime.datetime(
                           2018, 2, 8, 5, 25, 46),
                       datetime.datetime(2018, 2, 8, 5, 27, 47), datetime.datetime(
                           2018, 2, 8, 5, 29, 49),
                       datetime.datetime(2018, 2, 8, 5, 31, 39), datetime.datetime(
                           2018, 2, 8, 5, 33, 44),
                       datetime.datetime(2018, 2, 8, 5, 35, 49), datetime.datetime(
                           2018, 2, 8, 5, 37, 52),
                       datetime.datetime(2018, 2, 8, 5, 39, 56), datetime.datetime(
                           2018, 2, 8, 5, 42, 8),
                       datetime.datetime(2018, 2, 8, 5, 43, 58), datetime.datetime(
                           2018, 2, 8, 5, 46, 7),
                       datetime.datetime(2018, 2, 8, 5, 48, 14), datetime.datetime(
                           2018, 2, 8, 5, 50, 26),
                       datetime.datetime(2018, 2, 8, 5, 52, 34), datetime.datetime(
                           2018, 2, 8, 5, 54, 43),
                       datetime.datetime(2018, 2, 8, 5, 56, 42), datetime.datetime(
                           2018, 2, 8, 5, 58, 58),
                       datetime.datetime(2018, 2, 8, 6, 1, 8), datetime.datetime(
                           2018, 2, 8, 6, 3, 18),
                       datetime.datetime(2018, 2, 8, 6, 5, 29), datetime.datetime(
                           2018, 2, 8, 6, 7, 40),
                       datetime.datetime(2018, 2, 8, 6, 9, 33)], dtype=object)

    fitresults = wrap_fitNL(raw_data, exptimes, col_labels,
                            times=dtobjs, TrackFlux=True, subBgd=True)


if __name__ == '__main__':
    test_wrap_fitNL()
