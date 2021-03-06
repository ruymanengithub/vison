#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Trap-pumping Analysis Tools.

Created on Fri Mar 16 14:38:51 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
#import os
import numpy as np
import copy
from scipy import ndimage as nd
import pandas as pd
from collections import OrderedDict
from pylab import plot, show
from scipy.optimize import curve_fit
import string as st

from vison.datamodel import ccd as ccdmod
from vison.image.ds9reg import save_spots_as_ds9regs
# END IMPORT


def _thin_down_cat(cat, N):
    """ """
    ckeys = list(cat.keys())
    catlen = len(cat[ckeys[0]])
    if N >= catlen:
        return cat.copy()
    else:
        thcat = OrderedDict()
        ixsel = (np.random.choice(np.arange(catlen), N),)

        for key in ckeys:
            thcat[key] = cat[key][ixsel].copy()
        return thcat


def _get_N(dipdict):
    return len(dipdict['S'])


def _get_Ratio(dipdict):

    NN = len(np.where(dipdict['S'] == 1)[0])
    NS = len(np.where(dipdict['S'] == 0)[0])
    if NS > 0:
        Ratio = float(NN) / float(NS)
    else:
        Ratio = np.nan
    return Ratio


def _get_A(dipdict):
    return np.nanmedian(dipdict['A'])


def _aggregate_CQMT(masterdict, CCDs, allQuads, modkeys, tkeys, tkeyname):
    """ """

    allfuncts = OrderedDict()
    extracts = ['N', 'R', 'A']
    allfuncts['N'] = _get_N
    allfuncts['R'] = _get_Ratio
    allfuncts['A'] = _get_A

    #funct = allfuncts[extract]

    outdict = copy.deepcopy(masterdict)

    for CCDk in CCDs:
        for Q in allQuads:
            for modkey in modkeys:
                for tkey in tkeys:
                    values = []
                    for extract in extracts:
                        values.append(allfuncts[extract](masterdict[CCDk][Q][modkey][tkey]))
                        #values[extract] =  allfuncts[extract](masterdict[CCDk][Q][modkey][tkey])
                    outdict[CCDk][Q][modkey][tkey] = values

    reform = {(level1_key, level2_key, level3_key, level4_key): value
              for level1_key, level2_dict in list(outdict.items())
              for level2_key, level3_dict in list(level2_dict.items())
              for level3_key, level4_dict in list(level3_dict.items())
              for level4_key, value in list(level4_dict.items())}

    outdf = pd.DataFrame(reform).T
    colnames = dict()
    for i, c in enumerate(extracts):
        colnames[i] = c
    outdf.rename(columns=colnames, inplace=True)
    names = ['CCD', 'Q', 'mod', tkeyname]
    outdf.index.set_names(names, inplace=True)

    return outdf


def _get_tex(df):

    coretex = df.to_latex(multicolumn=True, multirow=True,
                          longtable=True, index=True,
                          escape=True)

    return coretex


def _get_txt(df):
    fok = df.to_string().split('\n')
    txt = ['\\begin{verbatim}'] + fok + ['\\end{verbatim}']
    return txt


def get_injprofile_tpnorm(ccdobj, vstart, vend):
    """Produces a 2D Map of charge injection to be used in trap-pumping analysis,
    to obtain dipole maps."""

    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    NAXIS1 = ccdobj.NAXIS1
    NAXIS2 = ccdobj.NAXIS2

    ccdobj_dummy = copy.deepcopy(ccdobj)

    for Q in Quads:

        qimg = ccdobj.extract_region(
            Q, 'img', vstart=vstart, vend=vend, canonical=True).copy()
        qimgmod = qimg.mean(axis=1).reshape(
            qimg.shape[0], 1).repeat(qimg.shape[1], axis=1)

        qmod = np.ones((NAXIS1 / 2, NAXIS2 / 2), dtype='float32')

        qmod[prescan:-overscan, vstart:vend] = qimgmod.copy()

        ccdobj_dummy.set_quad(qmod, Q, canonical=True, extension=-1)

    injprofile_tpnorm = ccdobj_dummy.extensions[-1].data.copy()

    return injprofile_tpnorm


def gen_raw_dpmap_vtpump(ccdobj, Navgrows=-1, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """ """

    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    for Q in Quads:

        model = get_InjProfile(ccdobj, Q, Navgrows=Navgrows,
                               vstart=vstart, vend=vend,
                               extension=extension)

        qdata = ccdobj.get_quad(Q, canonical=True, extension=extension).copy()

        dipmap = np.ones_like(qdata)
        dipmap[prescan:-overscan, vstart:vend] = \
            qdata[prescan:-overscan, vstart:vend] / model

        ccdobj.set_quad(dipmap, Q, canonical=True, extension=extension)

    return ccdobj


def gen_raw_dpmap_stpump(ccdobj, injprofiles, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """ """

    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    for Q in Quads:
        qdata = ccdobj.get_quad(Q, canonical=True, extension=extension).copy()

        dipmap = np.ones_like(qdata)
        dipmap[prescan:-overscan, vstart:vend] = \
            qdata[prescan:-overscan, vstart:vend] / injprofiles[Q]

        ccdobj.set_quad(dipmap, Q, canonical=True, extension=extension)

    return ccdobj


def get_InjProfile(ccdobj, Q, Navgrows=-1, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """ """

    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    qdata = ccdobj.get_quad(Q, canonical=True, extension=extension).copy()

    if Navgrows == -1:
        avgrow = qdata[prescan:-overscan, vstart:vend].mean(axis=1)
        model = avgrow[:, None]
    else:
        model = nd.filters.median_filter(qdata[prescan:-overscan, vstart:vend],
                                         size=(1, Navgrows),
                                         mode='reflect')
    return model


def wrap_gen_InjProfiles(ccdobj, Navgrows=-1, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """ """

    Quads = ccdobj.Quads
    InjProfiles = OrderedDict()

    for Q in Quads:
        InjProfiles[Q] = get_InjProfile(ccdobj, Q, Navgrows, vstart, vend, extension)

    return InjProfiles


def _fit_one_line(X, Y, pdeg=2):
    """ """

    pfit = np.polyfit(X, Y, pdeg)
    fpol = np.poly1d(pfit)
    Yfit = fpol(X)

    std = np.std(Y - Yfit)

    res = dict(coeffs=pfit, std=std)

    return res


def charact_injection(ccdobj):
    """ """

    Quads = ccdobj.Quads

    res = OrderedDict()
    for Q in Quads:
        res[Q] = OrderedDict()

    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    vstart = int(ccdobj.extensions[-1].header['vstart'])
    vend = min(int(ccdobj.extensions[-1].header['vend']), ccdobj.NrowsCCD)
    Nrows = vend - vstart

    Ncols = ccdobj.NcolsCCD

    X = np.arange(Nrows)
    pdeg = 2

    res['pdeg'] = pdeg

    for Q in Quads:

        quad = ccdobj.get_quad(Q, canonical=True, extension=-1)
        quadimg = quad[prescan:-overscan, vstart:vend].copy()

        quadpoly = np.zeros((Ncols, pdeg + 1), dtype='float32')
        quadinjnoise = np.zeros((Ncols,), dtype='float32')

        for i in range(Ncols):
            ires = _fit_one_line(X, quadimg[i, :], pdeg=pdeg)
            quadpoly[i, :] = ires['coeffs']
            quadinjnoise[i] = ires['std']

        res[Q]['polycoeffs'] = quadpoly.copy()
        res[Q]['injnoise'] = quadinjnoise.copy()

    return res


# def find_dipoles_vtpump_old(ccdobj, threshold, Q, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
#    """ """
#
#    prescan = ccdobj.prescan
#    overscan = ccdobj.overscan
#
#    qrawmap = ccdobj.get_quad(Q, canonical=True, extension=extension)[
#        prescan:-overscan, vstart:vend].copy()
#
#    asymmetry = 0.2
#
#    qraw_p1 = qrawmap[:, 1:].copy()
#    qraw_m1 = qrawmap[:, 0:-1].copy()
#    deltamap = qraw_p1-qraw_m1
#
#    rows0, cols0 = np.where((np.abs(deltamap) > threshold*2.) &
#                            (np.abs(qraw_m1-1.) > threshold) &
#                            (np.abs(qraw_p1-1.) > threshold) &
#                            ((2.*np.abs(qraw_m1-1.)-np.abs(deltamap)) < asymmetry*threshold) &
#                            ((2.*np.abs(qraw_p1-1.)-np.abs(deltamap)) < asymmetry*threshold))
#
#    if len(rows0) == 0:
#        return pd.DataFrame(dict(X=[], Y=[], S=[], A=[]), columns=['X', 'Y', 'S', 'A'])
#
#    X = rows0 + prescan
#    Y = cols0 + vstart
#
#    S = np.zeros_like(X)
#    S[np.where(deltamap[rows0, cols0] > 0)] = 1
#    S[np.where(deltamap[rows0, cols0] < 0)] = 0
#    A = deltamap[rows0, cols0]/2.
#
#    indict = OrderedDict(X=X, Y=Y, S=S, A=A)
#    df = pd.DataFrame(indict, columns=['X', 'Y', 'S', 'A'])
#
#    return df

def find_dipoles_vtpump(ccdobj, threshold, Q, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """Using Jesper Skottfelt's algorithm, as described in
    Trap_Pumping_Analysis_GCALCAMP_17OCT17_Azzollini.pdf

    'South' dipole: brighter pixel closer to serial register than dimmer pixel.
    'North' dipole: brighter pixel farther to serial register than dimmer pixel.

    """

    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    qrawmap = ccdobj.get_quad(Q, canonical=True, extension=extension)[
        prescan:-overscan, vstart:vend].copy()

    qmod = nd.filters.uniform_filter(qrawmap, size=21, mode='reflect')
    qresmap = qrawmap - qmod

    #qresmap = qrawmap - np.nanmedian(qrawmap)

    above_thresh = qresmap > threshold
    below_thresh = qresmap < -threshold
    within_thresh = ~above_thresh & ~below_thresh

    def rollit(array, disp):
        if disp == 0:
            return array
        return np.roll(array, disp, axis=1)

    dipoles_S = rollit(within_thresh, 0) & rollit(above_thresh, -1) & \
        rollit(below_thresh, -2) & rollit(within_thresh, -3)
    dipoles_N = rollit(within_thresh, 0) & rollit(below_thresh, -1) & \
        rollit(above_thresh, -2) & rollit(within_thresh, -3)

    def get_metrics(dipoles_mask, qresmap):
        NaxisY = dipoles_mask.shape[1]
        cols0, rows0 = np.where(dipoles_mask)
        X = cols0 + prescan
        Y = rows0 + vstart + 1
        ixnonrolled = np.where(Y < NaxisY - 4)
        A = (np.abs(rollit(qresmap, -1)) + np.abs(rollit(qresmap, -2)))[(cols0, rows0)] / 2.
        X = X[ixnonrolled]
        Y = Y[ixnonrolled]
        A = A[ixnonrolled]
        return X, Y, A

    XS, YS, AS = get_metrics(dipoles_S, qresmap)
    XN, YN, AN = get_metrics(dipoles_N, qresmap)

    Ndip = len(XS) + len(XN)
    S = np.zeros((Ndip,)) + np.nan

    S[0:len(XS)] = 1
    S[len(XS):None] = 0

    X = np.concatenate((XS, XN))
    Y = np.concatenate((YS, YN))
    A = np.concatenate((AS, AN))

    # if len(X) == 0:
    #    return pd.DataFrame(dict(X=[], Y=[], S=[], A=[]), columns=['X', 'Y', 'S', 'A'])

    outdict = OrderedDict()
    outdict['X'] = X.copy()
    outdict['Y'] = Y.copy()
    outdict['S'] = S.copy()
    outdict['A'] = A.copy()

    #df = pd.DataFrame(outdict, columns=['X', 'Y', 'S', 'A'])

    return outdict


def find_dipoles_stpump(ccdobj, threshold, Q, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """Using Jesper Skottfelt's algorithm, as described in
    Trap_Pumping_Analysis_GCALCAMP_17OCT17_Azzollini.pdf

    'West' dipole: brighter pixel closer to serial register than dimmer pixel.
    'East' dipole: brighter pixel farther to serial register than dimmer pixel.

    """

    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    qrawmap = ccdobj.get_quad(Q, canonical=True, extension=extension)[
        prescan:-overscan, vstart:vend].copy()

    qrawmap = qrawmap.mean(axis=1)
    qmod = nd.filters.uniform_filter(qrawmap, size=21, mode='reflect')
    qresmap = qrawmap - qmod

    #qresmap = qrawmap - np.nanmedian(qrawmap)

    above_thresh = qresmap > threshold
    below_thresh = qresmap < -threshold
    within_thresh = ~above_thresh & ~below_thresh

    def rollit(array, disp):
        if disp == 0:
            return array
        return np.roll(array, disp, axis=0)

    dipoles_W = rollit(within_thresh, 0) & rollit(above_thresh, -1) & \
        rollit(below_thresh, -2) & rollit(within_thresh, -3)
    dipoles_E = rollit(within_thresh, 0) & rollit(below_thresh, -1) & \
        rollit(above_thresh, -2) & rollit(within_thresh, -3)

    def get_metrics(dipoles_mask, qresmap):
        NaxisX = dipoles_mask.shape[0]
        cols0 = np.where(dipoles_mask)[0]
        X = cols0 + prescan
        ixnonrolled = np.where(X < NaxisX - 4)
        A = (np.abs(rollit(qresmap, -1)) + np.abs(rollit(qresmap, -2)))[cols0] / 2.
        X = X[ixnonrolled]
        A = A[ixnonrolled]
        return X, A

    XW, AW = get_metrics(dipoles_W, qresmap)
    XE, AE = get_metrics(dipoles_E, qresmap)

    Ndip = len(XW) + len(XE)
    S = np.zeros((Ndip,)) + np.nan

    S[0:len(XW)] = 1
    S[len(XW):None] = 0

    X = np.concatenate((XW, XE))
    A = np.concatenate((AW, AE))

    # if len(X) == 0:
    #    return pd.DataFrame(dict(X=[], S=[], A=[]), columns=['X', 'S', 'A'])

    #outdict = OrderedDict(X=X, S=S, A=A)
    outdict = OrderedDict()
    outdict['X'] = X.copy()
    outdict['S'] = S.copy()
    outdict['A'] = A.copy()

    return outdict


def save_dipcat2D_as_ds9regs(df, regfilename, clobber=True):
    """ """

    X = df['X']
    Y = df['Y']
    R = np.ones_like(X, dtype='float32') * 6.0

    data = dict(X=X,
                Y=Y,
                R=R)

    save_spots_as_ds9regs(data, regfilename, regtype='circle',
                          clobber=True)


def fcomp_distamp_dipoles_2D(merged, mcat):
    """ """

    _X = mcat['X'].as_matrix()
    _Y = mcat['Y'].as_matrix()
    _S = mcat['S'].as_matrix()

    _mix = np.zeros_like(_X, dtype='int32') - 1

    for ix in range(len(_X)):

        disc = np.isclose(((merged['uX'] - _X[ix])**2. + (merged['uY'] -
                                                          _Y[ix])**2.)**0.5, 0.) & (_S[ix] == merged['uS'])
        try:
            _mix[ix] = np.where(disc)[0][0]
        except BaseException:
            pass

    mcat['mix'] = pd.Series(_mix, index=mcat.index)

    return mcat


def fcomp_distamp_dipoles_1D(merged, mcat):
    """ """

    _X = mcat['X'].as_matrix()
    _S = mcat['S'].as_matrix()

    _mix = np.zeros_like(_X, dtype='int32') - 1

    for ix in range(len(_X)):

        disc = np.isclose(
            np.abs(merged['uX'] - _X[ix]), 0.) & (_S[ix] == merged['uS'])
        try:
            _mix[ix] = np.where(disc)[0][0]
        except BaseException:
            pass

    mcat['mix'] = pd.Series(_mix, index=mcat.index)

    return mcat


def merge_2dcats_generic(catsdict, catkeys, parentkey, columns, opcolumns, fcomp, dropna=False):
    """ """

    merged = catsdict[parentkey].copy()

    for col in opcolumns:
        merged['u%s' % col] = pd.Series(
            merged[col].as_matrix(), index=merged.index)

    pix = merged.index.values.copy()

    merged['mix'] = pd.Series(pix, index=merged.index)

    cats2match = [key for key in catkeys if key != parentkey]

    for ic, mkey in enumerate(cats2match):

        mcat = catsdict[mkey]

        #print mkey, len(mcat)

        mcat = fcomp(merged, mcat)

        merged = pd.merge(merged, mcat, how='outer', on='mix', suffixes=('', '_%s' % mkey),
                          copy=False)

        ixnonmerged = np.where(merged['mix'] == -1)
        Nnonmerged = len(ixnonmerged[0])
        maxmix = merged['mix'].max()
        if not np.isnan(maxmix):
            maxmix = int(maxmix)

        for ocol in opcolumns:
            merged.loc[ixnonmerged[0], 'u%s' %
                       ocol] = merged.loc[ixnonmerged[0], '%s_%s' % (ocol, mkey)]

        if (Nnonmerged > 0) and ~np.isnan(maxmix):

            merged.loc[ixnonmerged[0], 'mix'] = np.linspace(
                maxmix + 1, maxmix + Nnonmerged + 1, Nnonmerged, dtype='int32')

    renamerdict = dict()
    for key in columns:
        renamerdict[key] = '%s_%s' % (key, parentkey)

    merged = merged.rename(index=str, columns=renamerdict)
    merged = merged.drop(['mix'], axis=1)

    merged.drop(['u%s' % key for key in opcolumns], axis=1)

    if dropna:
        merged = merged.dropna()

    return merged


def merge_vtp_dipole_cats_bypos(catsdict, catkeys, parentkey, dropna=False):
    """ """
    columns = ['X', 'Y', 'S', 'A']
    opcolumns = ['X', 'Y', 'S']
    return merge_2dcats_generic(
        catsdict,
        catkeys,
        parentkey,
        columns,
        opcolumns,
        fcomp_distamp_dipoles_2D,
        dropna=False)


def merge_stp_dipole_cats_bypos(catsdict, catkeys, parentkey, dropna=False):
    columns = ['X', 'S', 'A']
    opcolumns = ['X', 'S']
    return merge_2dcats_generic(
        catsdict,
        catkeys,
        parentkey,
        columns,
        opcolumns,
        fcomp_distamp_dipoles_1D,
        dropna=False)


def get_f_A_vtp(N):
    def f_A_vtp(tph, logPc, tau):
        return N * 10.**logPc * (np.exp(-tph / tau) - np.exp(-2 * tph / tau))
    return f_A_vtp


def fit_PcTau_vtp(A, tois, Nshuffles=5000):
    """ """

    fitfunc = get_f_A_vtp(Nshuffles)

    p0 = [-6., tois[0]]
    pbounds = [[-7., -2.], [tois[0] / 10., tois[-1] * 10.]]

    ixsel = np.where(~np.isnan(A))
    if len(ixsel[0]) < 3:
        return np.nan, np.nan
    else:

        try:
            popt, pcov = curve_fit(fitfunc, tois[ixsel], A[ixsel], bounds=pbounds, p0=p0)

            Pc = 10.**popt[0]
            tau = popt[1]
        except BaseException:
            Pc, tau = np.nan, np.nan

    return Pc, tau


def batch_fit_PcTau_vtp(Amplitudes, tois, Nshuffles=5000):
    """ """

    Amp_mx = Amplitudes.as_matrix()

    Np = Amp_mx.shape[0]

    if Np == 0:
        return np.array([], dtype='float32'), np.array([], dtype='float32')

    assert Amp_mx.shape[1] == len(tois)

    Pc = np.zeros(Np, dtype='float32') + np.nan
    tau = np.zeros(Np, dtype='float32') + np.nan

    for i in range(Np):
        Pc[i], tau[i] = fit_PcTau_vtp(Amp_mx[i, :], tois, Nshuffles)

    return Pc, tau


def get_f_A_stp(N):
    def f_A_stp(tph, logPc, tau):
        return N * 10.**logPc * (1. - np.exp(-tph / tau))
    return f_A_stp


def fit_PcTau_stp(A, dwells, stoi, Nshuffles=5000):
    """ """

    fitfunc = get_f_A_stp(Nshuffles)

    p0 = [-6., stoi]
    pbounds = [[-7., -1.], [stoi / 10., stoi * 10.]]

    ixsel = np.where(~np.isnan(A))
    if len(ixsel[0]) < 3:
        return np.nan, np.nan
    else:

        try:
            popt, pcov = curve_fit(fitfunc, dwells[ixsel], A[ixsel], bounds=pbounds, p0=p0)

            Pc = 10.**popt[0]
            tau = popt[1]
        except BaseException:
            Pc, tau = np.nan, np.nan

    return Pc, tau


def batch_fit_PcTau_stp(Amplitudes, dwells, Nshuffles=5000):
    """ """

    stoi = 4.75  # us, serial TOI

    Amp_mx = Amplitudes.as_matrix()

    Np = Amp_mx.shape[0]

    if Np == 0:
        return np.array([], dtype='float32'), np.array([], dtype='float32')

    assert Amp_mx.shape[1] == len(dwells)

    Pc = np.zeros(Np, dtype='float32') + np.nan
    tau = np.zeros(Np, dtype='float32') + np.nan

    for i in range(Np):
        Pc[i], tau[i] = fit_PcTau_stp(Amp_mx[i, :], dwells + stoi, stoi, Nshuffles)
        # print(tau[i])

    return Pc, tau
