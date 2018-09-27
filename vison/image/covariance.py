#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Tools to retrieve covariance matrices for (differences of) Flat-Field images.
Used in the context of Brighter-Fatter analysis, mainly.

Created on Wed Mar  7 11:54:54 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from scipy import stats
# END IMPORT


def f_get_covmap(sq1, sq2, N, debug=False):
    """ """

    difimg = sq1 - sq2

    difimg -= np.median(difimg)  # ?

    NAXIS1 = difimg.shape[0]
    NAXIS2 = difimg.shape[1]

    def clip(img): return stats.sigmaclip(img, 4, 4).clipped

    var = np.nanvar(clip(difimg))
    mu = np.nanmedian(sq1)

    covmap = np.zeros((N, N), dtype='float32')

    x = np.arange(0, NAXIS1-N)
    y = np.arange(0, NAXIS2-N)

    x2d, y2d = np.meshgrid(x, y, indexing='ij')

    for i in range(N):
        for j in range(N):

            EXY = clip(difimg[x2d, y2d]*difimg[x2d+i, y2d+j]).mean()
            EX = clip(difimg[x2d, y2d]).mean()
            EY = clip(difimg[x2d+i, y2d+j]).mean()

            covmap[i, j] = (EXY-EX*EY) / var
            # covmap[i,j] = j # test

    # covmap[0,0] = 0. # not terribly interesting to see correlation of a pixel with itself

    if debug:
        stop()

    return covmap, mu, var


def get_cov_maps(ccdobjList, Npix=4, doTest=False):
    """ """

    Quads = ccdobjList[0].Quads

    Nframes = len(ccdobjList)
    Npairs = Nframes / 2

    tcovmapv = np.zeros((Npix, Npix, Npairs), dtype='float32')
    tmuv = np.zeros(Npairs, dtype='float32')
    tvarv = np.zeros(Npairs, dtype='float32')

    covmapv = dict()
    for Q in Quads:
        covmapv[Q] = tcovmapv.copy()
    muv = dict()
    for Q in Quads:
        muv[Q] = tmuv.copy()
    varv = dict()
    for Q in Quads:
        varv[Q] = tvarv.copy()

    if not doTest:

        for iP in range(Npairs):

            print 'Processing pair %i/%i' % (iP+1, Npairs)

            ccd1 = ccdobjList[iP*2]
            ccd2 = ccdobjList[iP*2+1]

            for Q in Quads:

                sq1 = ccd1.extract_region(
                    Q, area='img', canonical=True, extension=-1)
                sq2 = ccd2.extract_region(
                    Q, area='img', canonical=True, extension=-1)

                covmap, mu, var = f_get_covmap(sq1, sq2, Npix)

                covmapv[Q][:, :, iP] = covmap.copy()
                muv[Q][iP] = mu
                varv[Q][iP] = var

    av_var = dict()
    av_mu = dict()
    av_covmap = dict()

    for Q in Quads:
        av_covmap[Q] = covmapv[Q].mean(axis=2)
        av_var[Q] = np.mean(varv[Q])
        av_mu[Q] = np.mean(muv[Q])

    cov_dict = dict(varv=varv, covmapv=covmapv, muv=muv,
                    av_mu=av_mu, av_var=av_var,
                    av_covmap=av_covmap)

    cov_dict['Npairs'] = Npairs
    #cov_dict['Files'] = Files

    return cov_dict
