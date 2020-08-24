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
from pylab import plot, show  # TESTS
from sklearn import linear_model
from astropy.io import ascii
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit

from scipy import interpolate
# END IMPORT

FullDynRange = 2.**16
NLdeg = 4


def get_RANSAC_linear_model(X, Y):
    ransac = linear_model.RANSACRegressor()
    ransac.fit(np.expand_dims(X, 1), np.expand_dims(Y, 1))
    rpredictor = ransac.predict

    def predictor(scalar):
        return rpredictor(np.array([scalar]).reshape((1,1)))

    return predictor


def get_POLY_linear_model(X, Y):
    ixnonan = np.where(~np.isnan(X) & ~np.isnan(Y))
    mod1d_fit = np.polyfit(X[ixnonan], Y[ixnonan], 1)  # a linear approx. is fine
    predictor = np.poly1d(mod1d_fit)

    return predictor


def get_exptime_atfracdynrange(flu1D, exp1D, frac=0.5, method='spline',
                               maxrelflu=None,
                               debug=False):
    """ """

    X = flu1D / FullDynRange
    Y = exp1D.copy()
    if maxrelflu is None:
        ixvalid = np.where(~np.isnan(X))
    else:
        ixvalid = np.where(~np.isnan(X) & (X < maxrelflu))

    X = X[ixvalid].copy()
    Y = Y[ixvalid].copy()

    ixsort = np.argsort(X)
    X = X[ixsort].copy()
    Y = Y[ixsort].copy()

    if method == 'spline':
        predictor = interpolate.interp1d(X, Y, kind='linear')
        tfrac = predictor(frac)
    elif method == 'poly':
        predictor = get_POLY_linear_model(X, Y)
        tfrac = predictor(frac)
    elif method == 'ransac':
        predictor = get_RANSAC_linear_model(X, Y)
        tfrac = predictor(frac)[0, 0]

    if debug:
        try:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(X, Y, 'bo')
            Xplot = np.array([0., 1.])
            if method == 'ransac':
                Yplot = predictor(np.expand_dims(Xplot))[:, 0]
            else:
                Yplot = predictor(Xplot)
            ax.plot(Xplot, Yplot, 'k--')
            ax.axvline(frac, c='k', ls='--')
            ax.axhline(tfrac, c='r', ls='--')
            ax.set_xlabel('Flu/DynRange')
            ax.set_ylabel('exptime')
            plt.show()
        except BaseException:
            stop()

    return tfrac


def getXYW_NL(fluencesNL, exptimes, nomG, pivotfrac=0.5, maxrelflu=None, method='spline'):
    """ """

    assert fluencesNL.shape[0] == exptimes.shape[0]
    assert fluencesNL.ndim <= 2
    assert exptimes.ndim == 1

    Nexp = len(exptimes)

    if fluencesNL.ndim == 2:

        Nsec = fluencesNL.shape[1]
        #_exptimes = np.repeat(exptimes.reshape(Nexp,1),Nsec,axis=1)

        tpivot = np.zeros(Nsec, dtype='float32') + np.nan

        for i in range(Nsec):
            tpivot[i] = get_exptime_atfracdynrange(
                fluencesNL[:, i], exptimes, frac=pivotfrac, method=method, maxrelflu=maxrelflu, debug=False)
        tpivot = np.repeat(tpivot.reshape(1, Nsec), Nexp, axis=0)

        exptimes_bc = np.repeat(exptimes.reshape(Nexp, 1), Nsec, axis=1)

    else:

        tpivot = get_exptime_atfracdynrange(fluencesNL, exptimes,
                                            frac=pivotfrac,
                                            maxrelflu=maxrelflu,
                                            method=method)

        exptimes_bc = exptimes.copy()

    YL = exptimes_bc / tpivot * FullDynRange * pivotfrac
    Z = 100. * (fluencesNL / YL - 1.)

    efNL = np.sqrt(fluencesNL * nomG) / nomG

    W = 100. * (efNL / YL)

    X = fluencesNL.flatten().copy()
    Y = Z.flatten().copy()
    ixsort = np.argsort(X)
    X = X[ixsort].copy()
    Y = Y[ixsort].copy()
    W = W.flatten()[ixsort].copy()

    # if len(np.where(np.abs(Y)>5.)[0])>10: stop()# TESTS

    return X, Y, W


def getXYW_NL02(fluencesNL, exptimes, nomG, minrelflu=None, maxrelflu=None):
    """ """

    assert fluencesNL.shape[0] == exptimes.shape[0]
    assert fluencesNL.ndim <= 2
    assert exptimes.ndim == 1

    if maxrelflu is None:
        maxrelflu = 0.7
    if minrelflu is None:
        minrelflu = 0.02

    YL = np.zeros_like(fluencesNL, dtype='float32')
    if fluencesNL.ndim == 2:
        expix, regix = np.meshgrid(
            np.arange(
                fluencesNL.shape[0]), np.arange(
                fluencesNL.shape[1]), indexing='ij')
    else:
        expix = np.arange(fluencesNL.shape[0])
        regix = np.ones(fluencesNL.shape[0])

    uexptimes = np.unique(exptimes)

    for iu, uexp in enumerate(uexptimes):
        ix = np.where(exptimes == uexp)

        if fluencesNL.ndim == 2:
            iflu = fluencesNL[ix[0], ...].mean(axis=1)
        else:
            iflu = fluencesNL[ix[0]].copy()
        ixdeviants = np.where(sigma_clip(iflu, sigma=3).mask)
        if len(ixdeviants[0]) > 0:
            ixNaN = (ix[0][ixdeviants],)
            fluencesNL[ixNaN, ...] = np.nan

    if fluencesNL.ndim == 2:

        Nsec = fluencesNL.shape[1]
        #_exptimes = np.repeat(exptimes.reshape(Nexp,1),Nsec,axis=1)

        for isec in range(Nsec):
            #predictor = get_RANSAC_linear_model(exptimes[:],fluencesNL[:,isec])
            #Ypred = np.squeeze(predictor(np.expand_dims(exptimes,1)))

            ixnonan = np.where(~np.isnan(fluencesNL[:, 0]))

            _ixsel = np.where((fluencesNL[ixnonan, isec] >= 2.**16 * minrelflu) &
                              (fluencesNL[ixnonan, isec] <= 2.**16 * maxrelflu))

            ixsel = (ixnonan[0][_ixsel[1]],)

            xp = exptimes[ixsel]
            yp = np.squeeze(fluencesNL[ixsel, isec])
            predictor = get_POLY_linear_model(xp, yp)
            #predictor.coef[1] = 0.
            intersect = predictor.coef[1]
            YpredL = predictor(exptimes)
            YL[:, isec] = YpredL.copy()

            #print('sec:%i' % isec)
            # plot(exptimes[:],fluencesNL[:,isec]-YpredL,'k.')
            # plot(exptimes[:],fluencesNL[:,isec]-YpredL,'r-')
            # show()

            # plot(YpredL[3:],fluencesNL[3:,isec]/YpredL[3:]-1.,marker='.',ls='')
        # show()

    elif fluencesNL.ndim == 1:

        #predictor = get_RANSAC_linear_model(exptimes,fluencesNL)
        #YL[:] = np.squeeze(predictor(np.expand_dims(exptimes,1)))

        ixnonan = np.where(~np.isnan(fluencesNL))

        _ixsel = np.where((fluencesNL[ixnonan] >= 2.**16 * minrelflu) &
                          (fluencesNL[ixnonan] <= 2.**16 * maxrelflu))

        ixsel = (ixnonan[0][_ixsel],)

        #ixsel=np.where((fluencesNL>=2.**16*minrelflu) & (fluencesNL<=2.**16*maxrelflu))
        xp = exptimes[ixsel].copy()
        yp = np.squeeze(fluencesNL[ixsel]).copy()
        predictor = get_POLY_linear_model(xp, yp)

        #predictor.coef[1] = 0.
        intersect = predictor.coef[1]
        YpredL = predictor(exptimes)
        YL[:] = YpredL.copy()

        # plot(yp,yp/predictor(exptimes[ixsel])-1.,'k.')
        # show()

    Z = 100. * (fluencesNL - YL) / (YL - intersect)

    efNL = np.sqrt((fluencesNL - intersect) * nomG) / nomG

    W = 100. * (efNL / YL)

    ixsel = np.where((exptimes > 0.))

    # X = fluencesNL[ixsel].flatten().copy()
    X = fluencesNL[ixsel].flatten().copy()
    Y = Z[ixsel].flatten().copy()
    W = W[ixsel].flatten().copy()
    expix = expix[ixsel].flatten().copy()
    regix = regix[ixsel].flatten().copy()

    ixsort = np.argsort(X)
    X = X[ixsort].copy()
    Y = Y[ixsort].copy()
    W = W[ixsort].copy()
    expix = expix[ixsort].copy()
    regix = regix[ixsort].copy()

    return X, Y, W, expix, regix


def fitNL_pol(X, Y, W, Exptimes, minfitFl, maxfitFl, display=False):
    """ """

    ixnonan = np.where(~np.isnan(X))
    _selix = np.where((X[ixnonan] > minfitFl) & (X[ixnonan] < maxfitFl))
    selix = (ixnonan[0][_selix],)

    # TESTS
    doBin = False

    if doBin:

        Nbins = 30
        Xnanmin = np.nanmin(X[selix])
        Xnanmax = np.nanmax(X[selix])
        xfit = np.linspace(Xnanmin, Xnanmax, Nbins)
        ixbins = np.digitize(X[selix], xfit)
        yfit = np.array([np.median(Y[selix][np.where(ixbins == ix)]) for ix in range(Nbins)])
        ixnonans = np.where(~np.isnan(yfit))
        xfit = xfit[ixnonans].copy()
        yfit = yfit[ixnonans].copy()
    else:
        ixnonans = np.where(~(np.isnan(X) | np.isnan(Y)))
        xfit = X[ixnonans].copy()
        yfit = Y[ixnonans].copy()

    model = make_pipeline(PolynomialFeatures(NLdeg),
                          linear_model.Ridge(alpha=0.5,
                                             solver='auto', max_iter=1000))

    model.fit(xfit[:, np.newaxis], yfit[:, np.newaxis])

    # model.fit(X[selix],Y[selix])

    #x_plot = np.linspace(nanmin,nanmax,1000)[:,np.newaxis]
    #y_plot = model.predict(x_plot)
    NLfit = model._final_estimator.coef_[0][::-1]
    NLfit[-1] = model._final_estimator.intercept_[0]

    # plot(X[selix],Y[selix],'k.')
    # plot(xbins,ybins,'go')
    # plot(x_plot,y_plot,'r--')
    # show()

    #NLfit = np.polyfit(X[selix], Y[selix], w=W[selix], deg=NLdeg, full=False)

    NLpol = np.poly1d(NLfit)

    # array with NL fluences (1 to 2**16, in steps of 20 ADU)
    fkfluencesNL = np.linspace(minfitFl, maxfitFl, 200) * 1.
    Y_bestfit = NLpol(fkfluencesNL)

    # Based on fit
    #ixmax = np.abs(Y_bestfit).argmax()
    #maxNLpc = Y_bestfit[ixmax]
    #flu_maxNLpc = fkfluencesNL[ixmax]

    # Direct measure

    ixmax = np.abs(yfit).argmax()
    maxNLpc = yfit[ixmax]
    flu_maxNLpc = xfit[ixmax]

    if display:
        from matplotlib import pyplot as plt
        import matplotlib.cm as cm

        uexptimes, ixuexptimes = np.unique(Exptimes, return_index=True)

        ecolors = cm.rainbow(np.linspace(0, 1, len(uexptimes)))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(X, Y, 'k.')

        for i in selix[0]:
            ixcolor = np.where(np.isclose(uexptimes, Exptimes[i]))[0][0]
            ax.plot(X[i], Y[i], '.', color=ecolors[ixcolor, :])

        ax.plot(fkfluencesNL, Y_bestfit, 'r--')
        ax.set_xlabel('NL Fluence')
        ax.set_ylabel('NLpc')
        # ax.set_ylim([-10.,10.])
        plt.show()

    fitresults = OrderedDict(
        coeffs=NLfit,
        NLdeg=NLdeg,
        maxNLpc=maxNLpc,
        flu_maxNLpc=flu_maxNLpc,
        inputcurve=OrderedDict(
            X=X.copy(),
            Y=Y.copy()),
        outputcurve=OrderedDict(
            X=fkfluencesNL.copy(),
            Y=Y_bestfit.copy()))

    return fitresults


def fNL_wExp(x, *p):
    """ """
    return p[0] * np.exp(-(x - p[1]) / p[2]) + np.poly1d(p[3:])(x)


def fNL(x, *p):
    """ """
    return np.poly1d(p)(x)


def fitNL_taylored(X, Y, W, Exptimes, minfitFl, maxfitFl, display=False,
                   addExp=False):
    """ """

    ixnonan = np.where(~np.isnan(X))
    _selix = np.where((X[ixnonan] > minfitFl) & (X[ixnonan] < maxfitFl))
    selix = (ixnonan[0][_selix],)

    # TESTS
    doBin = False

    if doBin:

        Nbins = 30
        Xnanmin = np.nanmin(X[selix])
        Xnanmax = np.nanmax(X[selix])
        xfit = np.linspace(Xnanmin, Xnanmax, Nbins)
        ixbins = np.digitize(X[selix], xfit)
        yfit = np.array([np.median(Y[selix][np.where(ixbins == ix)]) for ix in range(Nbins)])
        ixnonans = np.where(~np.isnan(yfit))
        xfit = xfit[ixnonans].copy() / 2**16
        yfit = yfit[ixnonans].copy()
    else:
        #ixnonans = np.where(~(np.isnan(X) | np.isnan(Y)))
        xfit = X[selix].copy() / 2**16
        yfit = Y[selix].copy()

    ixsort = np.argsort(xfit)
    xfit = xfit[ixsort]
    yfit = yfit[ixsort]

    if addExp:
        #p0 = np.concatenate((np.array([10.,0.15,1.]),np.zeros(NLdeg+1)))
        p0 = np.concatenate((np.array([1., 0.01, 0.05]), np.zeros(NLdeg + 1)))
        bounds = []
        bounds.append([0., 10. / 2**16, 1.e-3] + [-100.] * (NLdeg) + [-10.])
        bounds.append([10., 10000. / 2**16, 2.E-1] + [100.] * (NLdeg) + [10.])
        # bounds = [[0.,  10., -1.E-3,-1.E-2, -10.],
        #          [1.E3,1.E3, 1.E-3, 1.E-2,  10.]]
    else:
        #p0 = np.concatenate((np.array([10.,0.15,1.]),np.zeros(NLdeg+1)))
        p0 = np.zeros(NLdeg + 1)
        bounds = []
        bounds.append([-100.] * (NLdeg) + [-10.])
        bounds.append([100.] * (NLdeg) + [10.])
        # bounds = [[0.,  10., -1.E-3,-1.E-2, -10.],
        #          [1.E3,1.E3, 1.E-3, 1.E-2,  10.]]

    if addExp:
        ff = fNL_wExp
    else:
        ff = fNL

    popt, pcov = curve_fit(ff, xfit, yfit,
                           p0=p0,
                           method='trf',
                           bounds=bounds,
                           absolute_sigma=False,
                           maxfev=10000)

    # array with NL fluences (1 to 2**16, in steps of 20 ADU)
    fkfluencesNL = np.linspace(minfitFl, maxfitFl, 400) * 1.
    Y_bestfit = ff(fkfluencesNL / 2**16, *popt)

    # Direct measure

    #ixmax = np.abs(yfit).argmax()
    #maxNLpc = yfit[ixmax]
    #flu_maxNLpc = xfit[ixmax]*2.**16

    # From the best fit curve

    ixmax = np.abs(Y_bestfit).argmax()
    maxNLpc = Y_bestfit[ixmax]
    flu_maxNLpc = fkfluencesNL[ixmax]

    if display:

        from matplotlib import pyplot as plt
        import matplotlib.cm as cm

        uexptimes, ixuexptimes = np.unique(Exptimes, return_index=True)

        ecolors = cm.rainbow(np.linspace(0, 1, len(uexptimes)))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(X, Y, 'k.')

        for i in selix[0]:
            ixcolor = np.where(np.isclose(uexptimes, Exptimes[i]))[0][0]
            ax.plot(X[i], Y[i], '.', color=ecolors[ixcolor, :])

        ax.plot(fkfluencesNL, Y_bestfit, 'r--')
        ax.set_xlabel('NL Fluence')
        ax.set_ylabel('NLpc')
        # ax.set_ylim([-10.,10.])
        plt.show()

    fitresults = OrderedDict(
        coeffs=popt,
        maxNLpc=maxNLpc,
        flu_maxNLpc=flu_maxNLpc,
        inputcurve=OrderedDict(
            X=X.copy(),
            Y=Y.copy()),
        outputcurve=OrderedDict(
            X=fkfluencesNL.copy(),
            Y=Y_bestfit.copy()))

    return fitresults


def wrap_fitNL_SingleFilter(fluences, variances, exptimes, times=np.array([]),
                            TrackFlux=True, subBgd=True):
    """ """
    # col001 == BGD
    # colEVEN = STAB
    # colODD = Fluences != 0

    nomG = 3.5  # e/ADU, used for noise estimates
    minfitFl = 2000.  # ADU
    maxfitFl = FullDynRange - 10000.  # ADU
    maxrelflu = 0.7  # used in determination of t_pivot

    NObsIDs, Nsecs = fluences.shape

    if TrackFlux:
        assert len(times) == len(exptimes)

        dtimes = np.array(
            [(times[i] - times[0]).seconds for i in range(NObsIDs)], dtype='float32')

    #col_numbers = np.array([int(item[3:]) for item in col_labels])

    nonzeroexptimes = exptimes[exptimes > 0.]
    unonzeroexptimes = np.unique(nonzeroexptimes)

    # the stability exposure time repeats the most
    _ixstab = np.array([(exptimes == iexptime).sum() for iexptime in unonzeroexptimes]).argmax()
    exptimestab = unonzeroexptimes[_ixstab]

    # boolean indices of the different types of exposures

    ixboo_bgd = exptimes == 0.
    ixboo_stab = exptimes == exptimestab
    ixboo_fluences = ((exptimes > 0.) & (exptimes != exptimestab))

    # Not used any more
    #ixboo_bgd = col_numbers == 1
    # ixboo_stab = (col_numbers % 2 == 0) & (col_numbers > 1)  # EVEN
    # ixboo_fluences = (col_numbers % 2 != 0)  & (col_numbers > 1)# ODD

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
        fluences[ixboo_fluences, :] /= track.reshape(track.shape[0], -1)
    else:
        track = np.ones_like(fluences[ixboo_fluences, 0])

    X, Y, W = getXYW_NL(fluences[ixboo_fluences, :],
                        exptimes[ixboo_fluences], nomG,
                        maxrelflu=maxrelflu,
                        method='ransac')

    fitresults = fitNL_pol(X, Y, W, minfitFl, maxfitFl, display=False)
    fitresults['bgd'] = bgd
    fitresults['stability_pc'] = np.std(track) * 100.

    return fitresults


def wrap_fitNL_TwoFilters(fluences, variances, exptimes, wave, times=np.array([]),
                          TrackFlux=True, subBgd=True, debug=False):
    """ """

    nomG = 3.5  # e/ADU, used for noise estimates
    minfitFl = 1000.  # ADU
    maxfitFl = FullDynRange - 10000.  # ADU
    pivotfrac = 0.2
    maxrelflu = 0.8

    NObsIDs, Nsecs = fluences.shape

    if TrackFlux:
        assert len(times) == len(exptimes)

        dtimes = np.array(
            [(times[i] - times[0]).seconds for i in range(NObsIDs)], dtype='float32')

    #col_numbers = np.array([int(item[3:]) for item in col_labels])

    nonzeroexptimes = exptimes[exptimes > 0.]
    unonzeroexptimes = np.unique(nonzeroexptimes)

    # the stability exposure time repeats the most
    _ixstab = np.array([(exptimes == iexptime).sum() for iexptime in unonzeroexptimes]).argmax()
    exptimestab = unonzeroexptimes[_ixstab]

    # boolean indices of the different types of exposures

    uwaves = np.unique(wave)

    ixboo_bgd = exptimes == 0.
    ixboo_stab = exptimes == exptimestab
    ixboo_fluA = ((exptimes > 0.) & (exptimes != exptimestab) & (wave == uwaves[0]))
    ixboo_fluB = ((exptimes > 0.) & (exptimes != exptimestab) & (wave == uwaves[1]))

    #ixboo_bgd = col_numbers == 1
    # ixboo_stab = (col_numbers % 2 == 0) & (col_numbers > 1)  # EVEN
    # ixboo_fluences = (col_numbers % 2 != 0)  & (col_numbers > 1)# ODD

    # assuming here the background is not affected by wavelength selection
    # in FW.. could be WRONG
    if subBgd:
        bgd = np.nanmean(fluences[ixboo_bgd, :], axis=0)
        bgd = np.repeat(bgd.reshape(1, Nsecs), NObsIDs, axis=0).copy()
        fluences -= bgd
    else:
        bgd = 0.

    if TrackFlux and len(times) == NObsIDs:
        st_dtimes = dtimes[ixboo_stab].copy()
        st_fluences = np.nanmean(fluences[ixboo_stab, :], axis=1).copy()

        track_fit = np.polyfit(st_dtimes, st_fluences, 2, full=False, cov=False)

        track_pol = np.poly1d(track_fit)

        trackA = track_pol(dtimes[ixboo_fluA])
        trackA /= np.median(trackA)

        trackB = track_pol(dtimes[ixboo_fluB])
        trackB /= np.median(trackB)

        #track_bc = np.repeat(track.reshape(len(track), 1), Nsecs, axis=1)
        fluences[ixboo_fluA, :] /= trackA.reshape(trackA.shape[0], -1)
        fluences[ixboo_fluB, :] /= trackB.reshape(trackB.shape[0], -1)
        trackstab = np.mean([trackA.std(), trackB.std()]) * 100.
    else:
        trackstab = 0.

    X_A, Y_A, W_A = getXYW_NL(fluences[ixboo_fluA, :],
                              exptimes[ixboo_fluA], nomG,
                              pivotfrac=pivotfrac,
                              maxrelflu=maxrelflu,
                              method='ransac')

    X_B, Y_B, W_B = getXYW_NL(fluences[ixboo_fluB, :],
                              exptimes[ixboo_fluB], nomG,
                              pivotfrac=pivotfrac,
                              maxrelflu=maxrelflu,
                              method='ransac')

    X = np.concatenate((X_A, X_B))
    Y = np.concatenate((Y_A, Y_B))
    W = np.concatenate((W_A, W_B))

    fitresults = fitNL_pol(X, Y, W, minfitFl, maxfitFl, display=debug)
    fitresults['bgd'] = bgd
    fitresults['stability_pc'] = trackstab

    return fitresults


def wrap_fitNL_TwoFilters_Alt(fluences, variances, exptimes, wave, times=np.array([]),
                              TrackFlux=True, debug=False, ObsIDs=None, NLdeg=NLdeg,
                              offset=0., XX=None, YY=None):
    """ """

    nomG = 3.5  # e/ADU, used for noise estimates
    minfitFl = 250.  # ADU
    maxfitFl = FullDynRange - 10000.  # ADU
    #pivotfrac = 0.2
    minrelflu = 0.10
    maxrelflu = 0.30

    NObsIDs, Nsecs = fluences.shape

    if TrackFlux:

        assert len(times) == len(exptimes)

        dtimes = np.array(
            [(times[i] - times[0]).seconds for i in range(NObsIDs)], dtype='float32')

    #col_numbers = np.array([int(item[3:]) for item in col_labels])

    nonzeroexptimes = exptimes[exptimes > 0.]
    unonzeroexptimes = np.unique(nonzeroexptimes)

    # hacky but effective way to identify stability exposures: the stability
    # exposure time repeats the most
    _ixstab = np.array([(exptimes == iexptime).sum() for iexptime in unonzeroexptimes]).argmax()
    exptimestab = unonzeroexptimes[_ixstab]

    # boolean indices of the different types of exposures

    uwaves = np.unique(wave)

    ixboo_bgd = exptimes == 0.
    ixboo_stab = exptimes == exptimestab
    ixboo_fluA = (exptimes != 0.) & (exptimes != exptimestab) & (wave == uwaves[0])
    ixboo_fluB = (exptimes != 0.) & (exptimes != exptimestab) & (wave == uwaves[1])

    if TrackFlux and len(times) == NObsIDs:

        st_dtimes = dtimes[ixboo_stab].copy()
        st_fluences = np.nanmean(fluences[ixboo_stab, :], axis=1).copy()

        ans = sigma_clip(st_fluences, sigma=3)
        ixgood = np.where(~ans.mask)

        # tmodel = make_pipeline(PolynomialFeatures(2),
        #                  linear_model.Ridge(alpha=0.1,
        #                  solver='auto',max_iter=1000))

        # tmodel.fit(st_dtimes[ixgood,np.newaxis],st_fluences[ixgood,np.newaxis])

        #track_fit = tmodel._final_estimator.coef_[0][::-1]
        #track_fit[-1] = tmodel._final_estimator.intercept_[0]

        track_fit = np.polyfit(st_dtimes[ixgood], st_fluences[ixgood], 2, full=False, cov=False)

        track_pol = np.poly1d(track_fit)

        trackA = track_pol(dtimes[ixboo_fluA])
        trackA /= np.median(trackA)

        trackB = track_pol(dtimes[ixboo_fluB])
        trackB /= np.median(trackB)

        #track_bc = np.repeat(track.reshape(len(track), 1), Nsecs, axis=1)
        fluences[ixboo_fluA, :] /= trackA.reshape(trackA.shape[0], -1)
        fluences[ixboo_fluB, :] /= trackB.reshape(trackB.shape[0], -1)
        trackstab = np.mean([trackA.std(), trackB.std()]) * 100.
    else:
        trackstab = 0.

    ixfitA = ixboo_fluA | ixboo_bgd | ixboo_stab
    X_A, Y_A, W_A, e_A, r_A = getXYW_NL02(fluences[ixfitA, :],
                                          exptimes[ixfitA], nomG,
                                          minrelflu=minrelflu,
                                          maxrelflu=maxrelflu)
    ixfitB = ixboo_fluB | ixboo_bgd
    X_B, Y_B, W_B, e_B, r_B = getXYW_NL02(fluences[ixfitB, :],
                                          exptimes[ixfitB], nomG,
                                          minrelflu=minrelflu,
                                          maxrelflu=maxrelflu)

    #bgdnoff = np.median(fluences[ixboo_bgd,:])

    X = np.concatenate((X_A, X_B)) - offset  # notice the small back-ground is left in!
    Y = np.concatenate((Y_A, Y_B))
    W = np.concatenate((W_A, W_B))
    #exps = np.concatenate((e_A,e_B))
    #regs = np.concatenate((r_A,r_B))

    Exptimes = np.concatenate((exptimes[ixfitA][e_A], exptimes[ixfitB][e_B]))

    Xcoo = np.concatenate((XX[e_A, r_A], XX[e_B, r_B]))
    Ycoo = np.concatenate((YY[e_A, r_A], YY[e_B, r_B]))

    fitresults = fitNL_taylored(X, Y, W, Exptimes, minfitFl, maxfitFl, display=debug,
                                addExp=True)
    fitresults['Xcoo'] = Xcoo.copy()
    fitresults['Ycoo'] = Ycoo.copy()

    #fitresults['bgd'] = bgd
    fitresults['stability_pc'] = trackstab

    # TESTS
    #import matplotlib.cm as cm
    #uexps = np.unique(exps)
    #colors = cm.rainbow(np.linspace(0,1,len(uexps)))
    # for iexp,uexp in enumerate(uexps):
    #    ixsel = np.where(uexp == exps)
    #    plot(X[ixsel],Y[ixsel],'.',c=colors[iexp])
    # show()
    # stop()

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

    exptimes = np.array([0., 0., 0., 0., 1.5, 1.5, 1.5, 1.5, 1.5,
                         15., 3., 3., 3., 3., 3., 15., 6., 6.,
                         6., 6., 6., 15., 9., 9., 9., 9., 9.,
                         15., 15., 15., 15., 15., 15., 15., 21., 21.,
                         21., 21., 21., 15., 24., 24., 24., 24., 24.,
                         15., 27., 27., 27., 27., 27., 15., 30., 30.,
                         30., 30., 30., 15., 33., 33., 33., 33., 33.,
                         15., 36., 36., 36., 36., 36., 15.])

    raw_data = np.repeat((exptimes / exptimes.max() * 2. **
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


def recalibrate_exptimes(exptimes, calibrationfile):
    """ """

    data = ascii.read(calibrationfile)
    commexptime = data['COMMEXPTIME'].data.copy()
    actexptime = data['ACTEXPTIME'].data.copy()

    predictor = interpolate.interp1d(commexptime, actexptime, kind='cubic')

    newexptimes = np.zeros_like(exptimes)
    ixnozero = np.where(exptimes != 0.)
    newexptimes[ixnozero] = predictor(exptimes[ixnozero])

    return newexptimes


if __name__ == '__main__':
    test_wrap_fitNL()
