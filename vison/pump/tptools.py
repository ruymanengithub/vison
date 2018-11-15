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

from vison.datamodel import ccd as ccdmod
from vison.image.ds9reg import save_spots_as_ds9regs
# END IMPORT


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

        qmod = np.ones((NAXIS1/2, NAXIS2/2), dtype='float32')

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
              qdata[prescan:-overscan, vstart:vend]/model

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


def find_dipoles_vtpump(ccdobj, threshold, Q, vstart=0, vend=ccdmod.NrowsCCD, extension=-1):
    """ """

    prescan = ccdobj.prescan
    overscan = ccdobj.overscan

    qrawmap = ccdobj.get_quad(Q, canonical=True, extension=extension)[
        prescan:-overscan, vstart:vend].copy()
    
    asymmetry = 0.2

    qraw_p1 = qrawmap[:, 1:].copy()
    qraw_m1 = qrawmap[:, 0:-1].copy()
    deltamap = qraw_p1-qraw_m1

    rows0, cols0 = np.where((np.abs(deltamap) > threshold*2.) &
                            (np.abs(qraw_m1-1.) > threshold) &
                            (np.abs(qraw_p1-1.) > threshold) &
                            ((2.*np.abs(qraw_m1-1.)-np.abs(deltamap)) < asymmetry*threshold) &
                            ((2.*np.abs(qraw_p1-1.)-np.abs(deltamap)) < asymmetry*threshold))

    if len(rows0) == 0:
        return pd.DataFrame(dict(X=[], Y=[], S=[], A=[]), columns=['X', 'Y', 'S', 'A'])

    X = rows0 + prescan
    Y = cols0 + vstart

    S = np.zeros_like(X)
    S[np.where(deltamap[rows0, cols0] > 0)] = 1
    S[np.where(deltamap[rows0, cols0] < 0)] = 0
    A = deltamap[rows0, cols0]/2.

    indict = OrderedDict(X=X, Y=Y, S=S, A=A)
    df = pd.DataFrame(indict, columns=['X', 'Y', 'S', 'A'])

    return df


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


def fcomp_distamp_dipoles(merged, mcat):
    """ """

    _X = mcat['X'].as_matrix()
    _Y = mcat['Y'].as_matrix()
    _S = mcat['S'].as_matrix()

    _mix = np.zeros_like(_X, dtype='int32') - 1

    for ix in range(len(_X)):
        disc = np.isclose(
            ((merged['uX']-_X[ix])**2.+(merged['uY']-_Y[ix])**2.)**0.5, 0.) & (_S[ix] == merged['uS'])
        try:
            _mix[ix] = np.where(disc)[0][0]
        except:
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
        maxmix = int(merged['mix'].max())

        for ocol in opcolumns:
            merged.loc[ixnonmerged[0], 'u%s' %
                       ocol] = merged.loc[ixnonmerged[0], '%s_%s' % (ocol, mkey)]

        merged.loc[ixnonmerged[0], 'mix'] = np.linspace(
            maxmix+1, maxmix+Nnonmerged+1, Nnonmerged, dtype='int32')

    renamerdict = dict()
    for key in columns:
        renamerdict[key] = '%s_%s' % (key, parentkey)

    merged = merged.rename(index=str, columns=renamerdict)
    merged = merged.drop(['mix'], axis=1)

    merged.drop([u'u%s' % key for key in opcolumns], axis=1)

    if dropna:
        merged = merged.dropna()

    return merged


def merge_vtp_dipole_cats_bypos(catsdict, catkeys, parentkey, dropna=False):
    """ """
    columns = ['X', 'Y', 'S', 'A']
    opcolumns = ['X', 'Y', 'S']
    return merge_2dcats_generic(catsdict, catkeys, parentkey, columns, opcolumns, fcomp_distamp_dipoles, dropna=False)
