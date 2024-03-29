#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to PSF0X.

Created on Tue Jan 30 16:31:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict
import string as st
import copy
from matplotlib import cm

from vison.datamodel import cdp
from vison.plot import figclasses
from vison.plot import trends
from vison.point import Paux
# END IMPORT


def get_check_offsets_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(stats=['offset_pre', 'offset_ove'],
                trendaxis='time',
                figname='%s_offset_vs_time.png' % (test,),
                caption='%s: offset vs. time.' % (ntest,),
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: offsets' % ntest,
                          ylim=trends.offset_lims))


def get_check_deltaoff_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(
        stats=[
            'deltaoff_pre', 'deltaoff_ove'], trendaxis='time', figname='%s_deltaoff_vs_time.png' %
        (test,), caption='%s: $\delta$offset vs. time. Offset value in each frame minus the average value.' %
        (ntest,), meta=dict(
            doLegend=True, doNiceXDate=True, suptitle='%s-checks: delta-offsets' %
            ntest, ylim=[
                -10., 10.]))


def get_check_std_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(stats=['std_pre', 'std_ove'],
                trendaxis='time',
                figname='%s_std_vs_time.png' % test,
                caption='%s: std vs. time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=True,
                          suptitle='%s-checks: std' % ntest,
                          ylim=trends.RON_lims))


def get_check_bgd_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(stats=['bgd_img'],
                trendaxis='time',
                figname='%s_bgd_vs_time.png' % test,
                caption='%s: Background vs. time.' % ntest,
                meta=dict(doLegend=False,
                          doNiceXDate=True,
                          suptitle='%s-checks: BGD' % ntest))


def get_check_flu_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(stats=['chk_fluence'],
                trendaxis='exptime',
                figname='%s_flu_vs_exptime.png' % test,
                caption='%s: Fluence vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: Fluence' % ntest,
                          xlabel='seconds',
                          ylabel='Flu.[ADU]'))


def get_check_fwhmx_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(stats=['chk_fwhmx'],
                trendaxis='exptime',
                figname='%s_fwhmx_vs_exptime.png' % test,
                caption='%s: FWHM(x) vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: FWHM(x)' % ntest,
                          xlabel='seconds',
                          ylabel='FWHMx [pix]'))


def get_check_fwhmy_dict(test):
    
    ntest = test.replace('_', '\_')
    return dict(stats=['chk_fwhmy'],
                trendaxis='exptime',
                figname='%s_fwhmy_vs_exptime.png' % test,
                caption='%s: FWHM(y) vs. Exposure Time.' % ntest,
                meta=dict(doLegend=True,
                          doNiceXDate=False,
                          suptitle='%s-checks: FWHM(y)' % ntest,
                          xlabel='seconds',
                          ylabel='FWHMy [pix]'))



def get_crosstalk_dict(test, figtype):
    tcaption = '%s: Cross-Talk [%s]. Green means positive cross-talk, red means negative cross-talk' +\
        ' (does not mean compliance/non-compliance). Pale colours mean less accurate results.'
    
    ntest = test.replace('_', '\_')
    crosstalk_dict = dict(
        figname='%s_crosstalk_%s.png' % (test, figtype),
        caption=tcaption % (ntest, figtype),
        meta=dict(),
        data=None
    )
    return crosstalk_dict

def get_spotsposter_dict(test, BFE=True):
    """ """
    if BFE:
        figtype ='noBFE'
        tcaption = '%s: Spots Poster, BFE not corrected. Log scale.'
    else:
        figtype = 'withBFE'
        tcaption = '%s: Spots Poster, BFE corrected using G+15. Log scale.'

    
    ntest = test.replace('_', '\_')
    sp_dict = dict(
        figname='%s_spotsposter_%s.png' % (test, figtype),
        caption=tcaption % (ntest,),
        meta=dict(doColorbar=True,
            corekwargs=dict(cmap=cm.gray, 
                            aspect='auto',  
                            # norm=None,
                            origin='lower left')),
        data=None
        )
    return sp_dict

def get_FWHM_v_flu_dict(test, fwhmkey):
    """ """
    ntest = test.replace('_', '\_')

    fdict = dict(
    figname='%s_%s_v_flu.png' % (test, fwhmkey),
    caption='%s: Gaussian-fit %s vs. Peak Fluence.' % 
        (ntest,fwhmkey.upper()),
    meta=dict(doLegend=True,
              ylabel='%s, [pix]' % fwhmkey,
              xlabel=r'$I_{0}\ [10\ kADU]$',
              ylim = [0.75,2.5],
              xlim = [0.,6.5],
              corekwargs=dict(
                  noBFE=dict(marker='', linestyle='--', color='b'),
                  BFE=dict(marker='', linestyle='-', color='r'),
                  ideal=dict(marker='',linestyle=':',color='k')),
              suptitle='%s: gaussian-fit %s in pixels vs. Fluence' %\
                    (ntest, fwhmkey))
        )
    return fdict


def get_skew_dict(test, vswhat, direction):
    """vswhat: 'position', 'fluence' 
    direction: 'x', 'y' 
    """
    ntest = test.replace('_', '\_')

    if vswhat == 'position':
        xlabel = 'CCD-%s Pos. [pix]' % direction
        suptitle = '%s: gaussian-fit res. %s-skew vs. %s-%s' %\
                    (ntest, direction, direction, vswhat)
        caption = '%s: Gaussian-fit residuals %s-skewness vs. %s-%s.' % \
        (ntest, direction, direction, vswhat)
    elif vswhat == 'fluence':
        xlabel = 'Fluence [kADU]'
        suptitle = '%s: gaussian-fit res. %s-skew vs. %s' %\
                    (ntest, direction, vswhat)
        caption = '%s: Gaussian-fit residuals %s-skewness vs. %s.' % \
        (ntest, direction, vswhat)

    fdict = dict(
    figname='%s_skew_%s_vs_%s.png' % (test, direction, vswhat),
    caption=caption,
    meta=dict(doLegend=True,
              doYErrbars=True,
              ylabel='%s-skewness [adim.]' % direction,
              xlabel=xlabel,
              #ylim = [0.75,2.5],
              #xlim = [0.,6.5],
              corekwargs=dict(
                  noBFE=dict(marker='.', linestyle='', color='b'),
                  BFE=dict(marker='.', linestyle='', color='r')),
              suptitle=suptitle)
        )
    return fdict


def get_PSF0Xfigs(test):
    PSF0Xfigs = dict()
    PSF0Xfigs['PSF0Xchecks_offsets'] = [
        trends.Fig_Basic_Checkstat, get_check_offsets_dict(test)]
    PSF0Xfigs['PSF0Xchecks_deltaoff'] = [
        trends.Fig_Basic_Checkstat, get_check_deltaoff_dict(test)]
    PSF0Xfigs['PSF0Xchecks_stds'] = [
        trends.Fig_Basic_Checkstat, get_check_std_dict(test)]
    PSF0Xfigs['PSF0Xchecks_bgd'] = [
        trends.Fig_Basic_Checkstat, get_check_bgd_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fluence'] = [
        trends.Fig_Basic_Checkstat, get_check_flu_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhmx'] = [
        trends.Fig_Basic_Checkstat, get_check_fwhmx_dict(test)]
    PSF0Xfigs['PSF0Xchecks_fwhmy'] = [
        trends.Fig_Basic_Checkstat, get_check_fwhmy_dict(test)]
    PSF0Xfigs['PSF0X_crosstalk_ADU'] = [
        figclasses.Fig_Husk, get_crosstalk_dict(test, 'ADU')]
    PSF0Xfigs['PSF0X_crosstalk_RATIO'] = [
        figclasses.Fig_Husk, get_crosstalk_dict(test, 'RATIO')]
    PSF0Xfigs['SpotsPoster'] = [
        figclasses.Fig_ImgShow, get_spotsposter_dict(test, BFE=False)]
    PSF0Xfigs['SpotsPosterNOBFE'] = [
        figclasses.Fig_ImgShow, get_spotsposter_dict(test, BFE=True)]
    PSF0Xfigs['PSF0X_fwhmx_v_flu'] = [
        figclasses.Fig_Beam2DPlot, get_FWHM_v_flu_dict(test, fwhmkey='fwhmx')]
    PSF0Xfigs['PSF0X_fwhmy_v_flu'] = [
        figclasses.Fig_Beam2DPlot, get_FWHM_v_flu_dict(test, fwhmkey='fwhmy')]


    PSF0Xfigs['PSF0X_skew_dirx_vs_fluence'] = [
        figclasses.Fig_Beam2DPlot, get_skew_dict(test, 'fluence', 'x')]
    PSF0Xfigs['PSF0X_skew_diry_vs_fluence'] = [
        figclasses.Fig_Beam2DPlot, get_skew_dict(test, 'fluence', 'y')]

    PSF0Xfigs['PSF0X_skew_dirx_vs_pos'] = [
        figclasses.Fig_Beam2DPlot, get_skew_dict(test, 'position', 'x')]
    PSF0Xfigs['PSF0X_skew_diry_vs_pos'] = [
        figclasses.Fig_Beam2DPlot, get_skew_dict(test, 'position', 'y')]    


    PSF0Xfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]
    return PSF0Xfigs


def get_PSF01_PANCHRO_figs():
    PSF01_PANCHRO_figs = dict()
    return PSF01_PANCHRO_figs


def get_CDP_lib(test):
    """ """
    CDP_lib = OrderedDict()

    CDP_lib['SPOTS'] = cdp.CDP()
    CDP_lib['SPOTS'].rootname = 'Spots'

    CDP_lib['SPOTS_NOBFE'] = cdp.CDP()
    CDP_lib['SPOTS_NOBFE'].rootname = 'Spots_nobfe'

    CDP_lib['RAW_CTALK'] = cdp.CDP()
    CDP_lib['RAW_CTALK'].rootname = 'Raw_crosstalk'

    CDP_lib['CTALK'] = cdp.CDP()
    CDP_lib['CTALK'].rootname = 'crosstalk'

    CDP_lib.update(Paux.get_CDP_lib())

    return CDP_lib



def _f_xy_bin(x,y,Nbins=3):
    """ """
    from sklearn.cluster import KMeans
    
    xresh = x.reshape(len(x),-1)
    kmeansRes = KMeans(n_clusters=Nbins, verbose=0).fit(xresh)

    xbin = np.zeros(Nbins,dtype='float32')
    ybin = np.zeros(Nbins,dtype='float32')
    ybinsig = np.zeros(Nbins,dtype='float32')
    vNbin = np.zeros(Nbins,dtype='int32')

    for i in range(Nbins):
        ixsel = np.where(kmeansRes.labels_ == i)
        xbin[i] = np.mean(x[ixsel])
        ybin[i] = np.mean(y[ixsel])
        ybinsig[i] = np.std(x[ixsel])
        vNbin[i] = len(ixsel[0])

    ixorder = np.argsort(xbin)

    xbin = xbin[ixorder]
    ybin = ybin[ixorder]
    ybinsig = ybinsig[ixorder]
    vNbin = vNbin[ixorder]

    return (xbin, ybin, ybinsig, vNbin)
    