#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

vison pipeline: Classes to do plots.

Created on Mon Nov 13 17:54:08 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""
# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from collections import OrderedDict
from vison.datamodel.HKtools import format_date
import copy
import itertools
import os
import string as st

from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 12
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['image.interpolation'] = 'none'
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
# END IMPORT


def _squiggle_xy(a, b, c, d, i=np.arange(0.0, 2*np.pi, 0.05)):
    """Dummy function used for Tests only."""
    return np.sin(i*a)*np.cos(i*b), np.sin(i*c)*np.cos(i*d)


class BasicPlot(object):

    def __init__(self, **kwargs):

        self.fig = None
        self.figsize = (9, 9)
        if 'fig' in kwargs:
            self.fig = kwargs['fig']
        self.axs = []

    def populate_axes(self):
        raise NotImplementedError("Subclass must implement abstract method")

    def plt_trimmer(self):
        raise NotImplementedError("Subclass must implement abstract method")

    def init_fig(self):
        plt.close('all')
        self.fig = plt.figure(figsize=self.figsize)

    def render(self, figname=''):

        self.init_fig()

        self.populate_axes()

        self.plt_trimmer()

        if figname == '':
            plt.show()
        else:
            plt.savefig(figname)

        plt.close()


class CCD2DPlot(BasicPlot):

    def __init__(self, data, **kwargs):

        super(CCD2DPlot, self).__init__(**kwargs)

        defaults = dict(suptitle='', qtitles=dict(E='E', F='F', G='G', H='H'),
                        doLegend=False,
                        doNiceXDate=False)

        if 'meta' in kwargs:
            meta = kwargs['meta']
        else:
            meta = dict()

        self.figsize = (8, 8)
        self.Quads = ['E', 'F', 'H', 'G']
        self.data = copy.deepcopy(data)
        self.meta = dict()
        self.meta.update(defaults)
        self.meta.update(meta)

        self.handles = []
        self.labels = []

    def populate_axes(self):

        qtitles = self.meta['qtitles']

        self.axs = []
        for iQ, Q in enumerate(self.Quads):

            self.axs.append(self.fig.add_subplot(2, 2, iQ+1))

            try:
                xkeys = self.data[Q]['x'].keys()
            except AttributeError:
                xkeys = None

            if xkeys is not None:
                ykeys = self.data[Q]['y'].keys()
                isconsistent = np.all([xkeys[i] == ykeys[i]
                                       for i in range(len(xkeys))])
                assert (len(xkeys) == len(ykeys)) and isconsistent

                for key in xkeys:
                    xarr = self.data[Q]['x'][key]
                    yarr = np.ones_like(self.data[Q]['y'][key])

                    handle = self.axs[-1].plot(xarr, yarr, label=key)
                    if iQ == 0:
                        self.handles += handle
                        self.labels.append(key)
            else:
                xarr = self.data[Q]['x']
                yarr = self.data[Q]['y']
                self.axs[-1].plot(xarr, yarr)

            self.axs[-1].set_title(qtitles[Q])

            if self.meta['doNiceXDate']:
                _xticks = self.axs[-1].get_xticks()
                if len(_xticks) > 6:
                    self.axs[-1].set_xticks(_xticks[::2])

    def plt_trimmer(self):

        if self.meta['doLegend']:
            plt.figlegend(self.handles, self.labels, loc='center right')

        if self.meta['doNiceXDate']:

            plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
            self.fig.autofmt_xdate()

        plt.subplots_adjust(top=0.80)
        plt.suptitle(self.meta['suptitle'])

        plt.tight_layout()

        if self.meta['doLegend']:
            plt.subplots_adjust(right=0.85)


class BeamPlot(BasicPlot):
    """ """

    def __init__(self, data, **kwargs):

        super(BeamPlot, self).__init__(**kwargs)

        meta = dict(suptitle='', ccdtitles=dict(CCD1='CCD1', CCD2='CCD2', CCD3='CCD3'),
                    doLegend=False,
                    doNiceXDate=False)
        meta.update(kwargs)

        self.figsize = (15, 6)
        self.Quads = ['E', 'F', 'H', 'G']
        self.CCDs = [1, 2, 3]
        self.data = copy.deepcopy(data)
        self.meta = meta.copy()
        self.handles = []
        self.labels = []
        self.fig = None
        self.axs = dict()

    def init_fig(self):
        self._init_fig_and_axes()

    def _init_fig_and_axes(self):
        """ """
        plt.close('all')
        fig, axsarr = plt.subplots(
            2, 6, sharex=True, sharey=True, figsize=self.figsize)
        self.fig = fig

        # initialisation of self.axs

        for CCD in self.CCDs:
            CCDkey = 'CCD%i' % CCD
            self.axs[CCDkey] = dict()
            for Q in self.Quads:
                self.axs[CCDkey][Q] = None

        self.axs['CCD1']['E'] = axsarr[0, 0]

        plotlist = [item for item in itertools.product(self.CCDs, ['E', 'F'])] +\
                   [item for item in itertools.product(self.CCDs, ['H', 'G'])]

        for k in range(1, len(plotlist)+1):
            CCDkey = 'CCD%i' % plotlist[k-1][0]
            Q = plotlist[k-1][1]
            self.axs[CCDkey][Q] = axsarr.flatten()[k-1]

    def populate_axes(self):
        """ """

        for CCD in self.CCDs:
            CCDkey = 'CCD%i' % CCD
            for Q in self.Quads:

                ax = self.axs[CCDkey][Q]
                CQdict = self.data[CCDkey][Q]

                try:
                    xkeys = CQdict['x'].keys()
                except AttributeError:
                    xkeys = None

                if xkeys is not None:
                    ykeys = CQdict['y'].keys()
                    isconsistent = np.all([xkeys[i] == ykeys[i]
                                           for i in range(len(xkeys))])
                    assert (len(xkeys) == len(ykeys)) and isconsistent

                    for key in xkeys:
                        handle, label = self._ax_core_funct(ax, CQdict, key)
                        if Q == 'E' and CCD == 1:
                            self.handles += handle
                            self.labels.append(label)
                else:
                    _, _ = self._ax_core_funct(ax, CQdict)

                if Q in ['E', 'H']:
                    ax.text(0.05, 0.9, Q, horizontalalignment='left',
                            transform=self.axs[CCDkey][Q].transAxes)
                elif Q in ['F', 'G']:
                    ax.text(0.9, 0.9, Q, horizontalalignment='right',
                            transform=self.axs[CCDkey][Q].transAxes)

                if Q == 'E':
                    ax.set_title(CCDkey, x=1)

                if self.meta['doNiceXDate']:
                    _xticks = ax.get_xticks()
                    if len(_xticks) > 6:
                        ax.set_xticks(_xticks[::2])

                if 'xlabel' in self.meta and Q in ['H', 'G']:
                    ax.set_xlabel(self.meta['xlabel'])
                if 'ylabel' in self.meta and Q in ['E', 'H'] and CCD == 1:
                    ax.set_ylabel(self.meta['ylabel'])
                
                if 'ylim' in self.meta:
                    ax.set_ylim(self.meta['ylim'])

            # self.axs[CCDkey][Q].locator_params(nticks=4,axis='x')

    def plt_trimmer(self):

        for CCD in self.CCDs:
            for Q in ['E', 'F']:
                plt.setp(self.axs['CCD%i' % CCD]
                         [Q].get_xticklabels(), visible=False)
            if CCD != 1:
                for Q in self.Quads:
                    plt.setp(self.axs['CCD%i' % CCD]
                             [Q].get_yticklabels(), visible=False)
            else:
                for Q in ['F', 'G']:
                    plt.setp(self.axs['CCD%i' % CCD]
                             [Q].get_yticklabels(), visible=False)

        if self.meta['doLegend']:
            plt.figlegend(self.handles, self.labels, loc='center right')

        if self.meta['doNiceXDate']:
            plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
            self.fig.autofmt_xdate()
            # plt.locator_params(nticks=4,axis='x',prune='both')

        plt.locator_params(axis='y', nticks=7, prune='both')
        if not self.meta['doNiceXDate']:
            plt.locator_params(axis='x', nticks=4, prune='both')

        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(wspace=0.0)

        plt.margins(0.05)

        plt.suptitle(self.meta['suptitle'])
        # plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        if self.meta['doLegend']:
            plt.subplots_adjust(right=0.85)


class BeamPlotYvX(BeamPlot):

    def _ax_core_funct(self, ax, CQdict, key=''):

        if key != '':
            xarr = CQdict['x'][key]
            yarr = CQdict['y'][key]
            label = st.replace(key, '_', '\_')
            handle = ax.plot(xarr, yarr, '-', label=label)
        else:
            xarr = CQdict['x']
            yarr = CQdict['y']
            ax.plot(xarr, yarr, '_')
            handle, label = None, None

        return handle, label

class BeamImgShow(BeamPlot):
    
    def _ax_core_funct(self, ax, CQdict):        
        ax.imshow(CQdict['img'])
        


class Beam1DHist(BeamPlot):
    """ """

    def _ax_core_funct(self, ax, CQdict, key=''):

        hist_kwargs = dict(range=None, density=False, weights=None,
                           cumulative=False, histtype='bar', align='mid',
                           orientation='vertical', log=False)
        for key in hist_kwargs.keys():
            if key in self.meta:
                hist_kwargs[key] = self.meta[key]

        if key != '':
            label = st.replace(key, '_', '\_')
            bins = CQdict['bins'][key]
            x = CQdict['x'][key]
            hist_kwargs['label'] = label
        else:
            bins = CQdict['bins']
            x = CQdict['x']
            label = None

        _, _, patch = ax.hist(x, bins, **hist_kwargs)

        return patch, label


class ImgShow(BasicPlot):

    def __init__(self, data, **kwargs):

        super(ImgShow, self).__init__(**kwargs)

        defaults = dict(title='')

        self.figsize = (7, 7)
        self.data = data
        self.meta = dict()
        self.meta.update(defaults)
        self.meta.update(kwargs)

    def populate_axes(self):
        """ """
        self.axs = [self.fig.add_subplot(111)]
        self.axs[0].imshow(self.data)
        self.axs[0].set_title(self.meta['title'])

    def plt_trimmer(self):
        """ """
        return


def testBeam2DPlot():

    ccddict = dict()
    for iQ, Q in enumerate(['E', 'F', 'G', 'H']):
        _x, _y = _squiggle_xy(iQ+1, iQ+1, iQ+2, iQ+2)
        xdict = OrderedDict(foo=_x, bar=_x*2.)
        ydict = OrderedDict(foo=_y, bar=_y*2.)
        ccddict[Q] = dict(x=copy.deepcopy(xdict), y=copy.deepcopy(ydict))

    data = dict(CCD1=copy.deepcopy(ccddict),
                CCD2=copy.deepcopy(ccddict),
                CCD3=copy.deepcopy(ccddict))

    meta = dict(suptitle='Test', doLegend=True)
    beam2dplot = BeamPlotYvX(data, meta=meta)

    beam2dplot.render(figname='')


def testBeam2ImgShow():
    """ """
    
    ccddict = dict()
    for iQ, Q in enumerate(['E', 'F', 'G', 'H']):
        #_x, _y = _squiggle_xy(iQ+1, iQ+1, iQ+2, iQ+2)
        #xdict = OrderedDict(foo=_x, bar=_x*2.)
        #ydict = OrderedDict(foo=_y, bar=_y*2.)
        x = np.linspace(0.,1.,10)
        y = np.linspace(0.,1.,10)
        xx, yy = np.meshgrid(x,y)
        img = np.sqrt(xx**2.+yy**2.)
        ccddict[Q] = dict(img=img)

    data = dict(CCD1=copy.deepcopy(ccddict),
                CCD2=copy.deepcopy(ccddict),
                CCD3=copy.deepcopy(ccddict))

    meta = dict(suptitle='Test')
    beamimgshow = BeamImgShow(data, meta=meta)

    beamimgshow.render(figname='')



if __name__ == '__main__':

    #testBeam2DPlot()
    testBeam2ImgShow()
