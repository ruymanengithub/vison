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

def fclipsig(img, clipsigma=None):
    
    if clipsigma is not None:
        return stats.sigmaclip(img, clipsigma, clipsigma).clipped
    else:
        return img.flatten().copy()

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

    difimg -= np.nanmedian(difimg)  # ?

    if submodel:
        #t1 = time()
        model2d = get_model2d(difimg, pdegree=5, doBin=True, binsize=300, doFilter=False)
        #t2 = time()
        #print '%.1f seconds in computing 2D model...' % (t2-t1,)
        difimg -= model2d.imgmodel

    NAXIS1 = difimg.shape[0]
    NAXIS2 = difimg.shape[1]


    vardif = np.nanvar(fclipsig(difimg, clipsigma))
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

            corrmap[i, j] = (EXY - EX * EY) / vardif
            # corrmap[i,j] = j # test

    # better keep corrmap[0,0] for debugging purposes
    # corrmap[0,0] = 0. # not terribly interesting to see correlation of a pixel with itself

    if debug:
        stop()

    return corrmap, mu, vardif

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

    difimg -= np.nanmedian(difimg)  # ?

    if submodel:
        #t1 = time()
        model2d = get_model2d(difimg, pdegree=5, doBin=True, binsize=300, doFilter=False)
        #t2 = time()
        #print '%.1f seconds in computing 2D model...' % (t2-t1,)
        difimg -= model2d.imgmodel

    NAXIS1 = difimg.shape[0]
    NAXIS2 = difimg.shape[1]
    
    ixnan = np.where(difimg.mask==True)
    msksq1 = sq1.data.copy()
    msksq1[ixnan] = np.nan
    msksq2 = sq2.data.copy()
    msksq2[ixnan] = np.nan
    mskdifimg = difimg.data.copy()
    mskdifimg[ixnan] = np.nan

    _mskdifimg = mskdifimg.copy()
    _mskdifimg[ixnan] = clipsigma * 10 * np.nanstd(mskdifimg)

    vardif = np.nanvar(fclipsig(_mskdifimg, clipsigma))
    mu = np.nanmean([np.nanmedian(msksq1),np.nanmedian(msksq2)])

    corrmap = np.zeros((N, N), dtype='float32')

    x = np.arange(0, NAXIS1 - N)
    y = np.arange(0, NAXIS2 - N)

    x2d, y2d = np.meshgrid(x, y, indexing='ij')


    if estimator == 'median':
        festimator = np.nanmedian
    elif estimator == 'mean':
        festimator = np.nanmean
    
    difimg0 = mskdifimg[x2d, y2d].copy()
    difimg0 -= np.nanmean(difimg0)

    for i in range(N):
        for j in range(N):

            difimgij = mskdifimg[x2d + i, y2d + j]
            difimgij -= np.nanmean(difimgij)

            product = difimg0 * difimgij
            product[np.where(np.isnan(product))] = clipsigma*10.*vardif

            corrmap[i,j] = festimator(fclipsig(product, 
                    clipsigma)) / vardif

    if debug:
        stop()
    
    return corrmap, mu, vardif


def f_get_corrmap_tests(sq1, sq2, N, submodel=False, estimator='median', clipsigma=4.,
        debug=False):
    """ """
    from time import time
    #from astropy.io import fits as fts
    # BY-PASS on TESTS
    # corrmap = np.zeros((N,N),dtype='float32')+0.01 # TESTS
    #corrmap[0,0] = 1.-0.01
    #mu = 1.
    #var = 1.
    # return corrmap, mu, var

    #difimg = sq1 - sq2

    #difimg -= np.nanmedian(difimg)  # ?

    def f_submodel(img):
        model2d = get_model2d(img, pdegree=5, doBin=True, binsize=300, doFilter=False)
        img -= model2d.imgmodel
        fluence = np.nanmean(model2d.imgmodel)
        return img, fluence

    if submodel:

        t1 = time()
        ssq1, flu1 = f_submodel(sq1)
        ssq2, flu2 = f_submodel(sq2)
        t2 = time()
        print(('%.1f seconds in computing 2D models...' % (t2-t1,)))

    mu = np.mean([flu1,flu2])

    difimg = ssq1 * mu/flu1 - ssq2 * (mu/flu2)
    
    NAXIS1 = difimg.shape[0]
    NAXIS2 = difimg.shape[1]

    mskdifimg = difimg.data.copy()
    ixnan = np.where(difimg.mask==True)
    mskdifimg[ixnan] = np.nan


    _mskdifimg = mskdifimg.copy()
    _mskdifimg[ixnan] = clipsigma * 10 * np.nanstd(mskdifimg) 
    # setting masked values to an arbitrarily large value so they get clipped
    # clipping can't work with nans..

    vardif = np.nanvar(fclipsig(_mskdifimg, clipsigma))

    corrmap = np.zeros((N, N), dtype='float32')

    x = np.arange(0, NAXIS1 - N)
    y = np.arange(0, NAXIS2 - N)

    x2d, y2d = np.meshgrid(x, y, indexing='ij')

    if estimator == 'median':
        festimator = np.nanmedian
    elif estimator == 'mean':
        festimator = np.nanmean
    
    difimg0 = mskdifimg[x2d, y2d].copy()
    difimg0 -= np.nanmean(difimg0)

    for i in range(N):
        for j in range(N):

            difimgij = mskdifimg[x2d + i, y2d + j]
            difimgij -= np.nanmean(difimgij)

            product = difimg0 * difimgij
            product[np.where(np.isnan(product))] = clipsigma*10.*vardif

            corrmap[i,j] = festimator(fclipsig(product, 
                    clipsigma)) / vardif

    if debug:
        stop()

    return corrmap, mu, vardif

def get_sigmaclipcorr(vardif, clipsigma, estimator, dims=None):
    """ """

    if dims is None:
        dims = (ccdmod.NcolsCCD,ccdmod.NrowsCCD)

    if estimator == 'median':
        festimator = np.nanmedian
    elif estimator == 'mean':
        festimator = np.nanmean

    QDiffnoise = np.random.normal(
            loc=0.0, scale=vardif**0.5, 
            size=dims)

    
    QDiffnoise -= np.nanmean(QDiffnoise)

    clipvar = festimator(fclipsig(QDiffnoise * QDiffnoise, 
                    clipsigma))
    clipcorr = clipvar / vardif
    #print('var= %.3f, clipcorr = %.3f' % (vardif, clipcorr))

    return clipcorr


def get_cov_maps(ccdobjList, Npix=4, vstart=0, vend=2066, 
        clipsigma=4., covfunc='ver1', doBiasCorr=False, 
        central='median', doTest=False, debug=False):
    """ """

    maxbadfrac = 0.2
    estimator = central
    corefuncs = dict(ver1=f_get_corrmap,
        ver2=f_get_corrmap_v2,
        tests=f_get_corrmap_tests)
    corefunc = corefuncs[covfunc]
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


    nomshape = (ccdobjList[0].NcolsCCD,ccdobjList[0].NrowsCCD)
    if doBiasCorr:
        corrclip = get_sigmaclipcorr(1., clipsigma, estimator,
            dims = nomshape)
    else:
        corrclip = 1.

    if not doTest:

        for iP in range(Npairs):

            print(('Processing pair %i/%i' % (iP + 1, Npairs)))

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
                        #mskdsq1 = sq1.data.copy()
                        #mskdsq1[np.where(sq1.mask==True)] = np.nan

                        #mskdsq2 = sq2.data.copy()
                        #mskdsq2[np.where(sq2.mask==True)] = np.nan                        

                        corrmap, mu, vardif = corefunc(sq1, sq2, Npix, 
                                                submodel=True,
                                                estimator=estimator,
                                                clipsigma=clipsigma,
                                                debug=debug)

                        if doBiasCorr:
                            corrmap /= corrclip
                        
                        corrmapv[Q][:, :, iP] = corrmap.copy()
                        muv[Q][iP] = mu
                        varv[Q][iP] = vardif / 2.

                        #stop()

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
