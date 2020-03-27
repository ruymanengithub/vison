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
from vison.datamodel.ccd_aux import Model2D
from vison.datamodel import ccd as ccdmod
# END IMPORT


def get_model2d(img, pdegree=5, doFilter=False, doBin=True,
                filtsize=1, binsize=1, filtertype='mean'):

    regmodel = Model2D(img)

    if doFilter:
        if filtsize > 1:
            regmodel.filter_img(filtsize=filtsize, filtertype=filtertype,
                                Tests=False)
    if doBin:
        regmodel.bin_img(boxsize=binsize, stat=filtertype)

    regmodel.get_model_poly2D(
        sampling=filtsize, pdegree=pdegree, useBin=doBin)

    return regmodel

def fclipsig(img, clipsigma): 
    return stats.sigmaclip(img, clipsigma, clipsigma).clipped

def f_get_corrmap(sq1, sq2, N, submodel=False, estimator='median', clipsigma=4.,
        debug=False):
    """ """
    #from time import time

    # BY-PASS on TESTS
    # corrmap = np.zeros((N,N),dtype='float32')+0.01 # TESTS
    #corrmap[0,0] = 1.-0.01
    #mu = 1.
    #var = 1.
    # return corrmap, mu, var

    difimg = sq1 - sq2

    difimg -= np.median(difimg)  # ?

    if submodel:
        #t1 = time()
        model2d = get_model2d(difimg, pdegree=5, doBin=True, binsize=300, doFilter=False)
        #t2 = time()
        #print '%.1f seconds in computing 2D model...' % (t2-t1,)
        difimg -= model2d.imgmodel

    NAXIS1 = difimg.shape[0]
    NAXIS2 = difimg.shape[1]


    var = np.nanvar(fclipsig(difimg, clipsigma))
    mu = np.nanmedian(sq1)

    corrmap = np.zeros((N, N), dtype='float32')

    x = np.arange(0, NAXIS1 - N)
    y = np.arange(0, NAXIS2 - N)

    x2d, y2d = np.meshgrid(x, y, indexing='ij')


    if estimator == 'median':
        festimator = np.nanmedian
    elif estimator == 'mean':
        festimator = np.nanmean
    

    for i in range(N):
        for j in range(N):

            EXY = festimator(fclipsig(difimg[x2d, y2d] * difimg[x2d + i, y2d + j], 
                            clipsigma))
            EX = festimator(fclipsig(difimg[x2d, y2d], clipsigma))
            EY = festimator(fclipsig(difimg[x2d + i, y2d + j], clipsigma))

            corrmap[i, j] = (EXY - EX * EY) / var
            # corrmap[i,j] = j # test

    # better keep corrmap[0,0] for debugging purposes
    # corrmap[0,0] = 0. # not terribly interesting to see correlation of a pixel with itself

    if debug:
        stop()

    return corrmap, mu, var

def f_get_corrmap_v2(sq1, sq2, N, submodel=False, estimator='median', clipsigma=4.,
        debug=False):
    """ """
    #from time import time

    # BY-PASS on TESTS
    # corrmap = np.zeros((N,N),dtype='float32')+0.01 # TESTS
    #corrmap[0,0] = 1.-0.01
    #mu = 1.
    #var = 1.
    # return corrmap, mu, var

    difimg = sq1 - sq2

    difimg -= np.median(difimg)  # ?

    if submodel:
        #t1 = time()
        model2d = get_model2d(difimg, pdegree=5, doBin=True, binsize=300, doFilter=False)
        #t2 = time()
        #print '%.1f seconds in computing 2D model...' % (t2-t1,)
        difimg -= model2d.imgmodel

    NAXIS1 = difimg.shape[0]
    NAXIS2 = difimg.shape[1]

    var = np.nanvar(fclipsig(difimg, clipsigma))
    mu = np.nanmedian(sq1)

    corrmap = np.zeros((N, N), dtype='float32')

    x = np.arange(0, NAXIS1 - N)
    y = np.arange(0, NAXIS2 - N)

    x2d, y2d = np.meshgrid(x, y, indexing='ij')


    if estimator == 'median':
        festimator = np.nanmedian
    elif estimator == 'mean':
        festimator = np.nanmean
    
    difimg0 = difimg[x2d, y2d].copy()
    difimg0 -= np.nanmean(difimg[x2d, y2d])

    for i in range(N):
        for j in range(N):

            difimgij = difimg[x2d + i, y2d + j]
            difimgij -= np.nanmean(difimgij)

            corrmap[i,j] = festimator(fclipsig(difimg0 * difimgij, 
                    clipsigma)) / var


    if debug:
        stop()

    return corrmap, mu, var

def get_sigmaclipcorr(var, clipsigma, estimator, dims=None):
    """ """

    if dims is None:
        dims = (ccdmod.NcolsCCD,ccdmod.NrowsCCD)

    QDiffnoise = np.random.normal(
            loc=0.0, scale=np.sqrt(var*2.), 
            size=dims)

    sigvar = np.nanvar(QDiffnoise)

    return sigvar / var


def get_cov_maps(ccdobjList, Npix=4, vstart=0, vend=2066, 
        clipsigma=4., doTest=False, debug=False):
    """ """

    maxbadfrac = 0.2
    estimator = 'median'
    corefunc = f_get_corrmap
    doCorrClip = False
    Quads = ccdobjList[0].Quads

    Nframes = len(ccdobjList)
    Npairs = Nframes / 2

    tcorrmapv = np.zeros((Npix, Npix, Npairs), dtype='float32') + np.nan
    tmuv = np.zeros(Npairs, dtype='float32') + np.nan
    tvarv = np.zeros(Npairs, dtype='float32') + np.nan

    corrmapv = dict()
    for Q in Quads:
        corrmapv[Q] = tcorrmapv.copy()
    muv = dict()
    for Q in Quads:
        muv[Q] = tmuv.copy()
    varv = dict()
    for Q in Quads:
        varv[Q] = tvarv.copy()

    if not doTest:

        for iP in range(Npairs):

            print('Processing pair %i/%i' % (iP + 1, Npairs))

            ccd1 = ccdobjList[iP * 2]
            ccd2 = ccdobjList[iP * 2 + 1]

            for Q in Quads:

                sq1 = ccd1.extract_region(
                    Q, area='img', canonical=True, vstart=vstart,
                    vend=vend, extension=-1)
                sq2 = ccd2.extract_region(
                    Q, area='img', canonical=True, vstart=vstart,
                    vend=vend, extension=-1)

                badfrac1 = len(np.where(sq1.mask)[0]) / float(sq1.size)
                badfrac2 = len(np.where(sq2.mask)[0]) / float(sq2.size)

                if (badfrac1 > maxbadfrac) or (badfrac2 > maxbadfrac):
                    continue
                else:
                    try:
                        # if Q=='F': debug=True #TEST
                        corrmap, mu, var = corefunc(sq1, sq2, Npix, submodel=True,
                                                        estimator=estimator,
                                                        clipsigma=clipsigma,
                                                       debug=debug)

                        if doCorrClip:
                            corrclip = get_sigmaclipcorr(var, clipsigma, estimator,
                                dims = sq1.shape)
                            corrmap /= corrclip

                        corrmapv[Q][:, :, iP] = corrmap.copy()
                        muv[Q][iP] = mu
                        varv[Q][iP] = var

                        



                    except BaseException:
                        pass

    av_var = dict()
    av_mu = dict()
    av_corrmap = dict()

    for Q in Quads:
        av_corrmap[Q] = np.nanmean(corrmapv[Q], axis=2)
        av_var[Q] = np.nanmean(varv[Q])
        av_mu[Q] = np.nanmean(muv[Q])

    cov_dict = dict(varv=varv, corrmapv=corrmapv, muv=muv,
                    av_mu=av_mu, av_var=av_var,
                    av_corrmap=av_corrmap)

    cov_dict['Npairs'] = Npairs
    #cov_dict['Files'] = Files

    return cov_dict
