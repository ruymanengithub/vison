#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Classes to handle FPA-level plots

Created on Fri Jul  5 16:16:20 2019

@author: raf
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

from vison.fpa import fpa as fpamod

import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot as plt
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

from vison.plot.baseplotclasses import BasicPlot
# END IMPORT



class FpaPlot(BasicPlot):
    """ """

    def __init__(self, data, **kwargs):

        super(FpaPlot, self).__init__(**kwargs)
        
        self.figsize = (15,15)

        meta = dict(suptitle='', 
                    #ccdtitles=dict(CCD1='CCD1', CCD2='CCD2', CCD3='CCD3'),
                    doLegend=False,
                    doColorbar=False,
                    doNiceXDate=False,
                    labelCCDs=False,
                    doRotateXLabels=False)
        meta.update(kwargs)
        
        self.NSlices = 6
        self.Ncols = 6 
        self.Quads = ['E', 'F', 'H', 'G']
        self.CCDs = [1, 2, 3]
        
        self.data = copy.deepcopy(data)
        self.meta = dict()
        self.meta.update(meta)
        self.handles = []
        self.labels = []
        self.fig = None
        self.axs = dict()
        self.axsarr = []
        
        self.corekwargs = dict()
        if 'corekwargs' in kwargs:
            self.corekwargs.update(kwargs['corekwargs']) 

    def init_fig(self):
        self._init_fig_and_axes()

    def _init_fig_and_axes(self):
        """ """
        plt.close('all')
        self.fig, self.axsarr = plt.subplots(
            6, 6, sharex=True, sharey=True, figsize=self.figsize)
                
        NY = self.Ncols
        NX = self.NSlices
        
        # initialisation of self.axs
        
        for jY in range(NY):
            for iX in range(NX):    
                Ckey = 'C_%i%i' % (jY+1,iX+1)
                self.axs[Ckey] = self.axsarr[NY-1-jY,iX]

            
    def _ax_core_funct(self, ax, CQdict, key=''):
        """ """
        raise NotImplementedError("Subclass must implement abstract method")

    def populate_axes(self):
        """ """
        
        try:
            labelkeys = self.data['labelkeys']
        except KeyError:
            labelkeys = []
            
        NX = self.Ncols
        NY = self.NSlices
        
        for jY in range(NY):
            for iX in range(NX):   
                Ckey = 'C_%i%i' % (jY+1,iX+1)
                
                ax = self.axs[Ckey]
                CCdict = self.data[Ckey]
                

                if len(labelkeys) > 0:
                    for labelkey in labelkeys:
                        handle, label = self._ax_core_funct(ax, CCdict, labelkey)
                        if Ckey == 'C_11':
                            self.handles += handle
                            self.labels.append(label)
                else:
                    _, _ = self._ax_core_funct(ax, CCdict)
                
                
                if self.meta['labelCCDs']:
                    rCkey = st.replace(Ckey,'_','')
                    ax.text(0.05,0.85,rCkey,horizontalalignment='left',
                            transform=self.axs[Ckey].transAxes)


                if self.meta['doNiceXDate']:
                    _xticks = ax.get_xticks()
                    if len(_xticks) > 6:
                        ax.set_xticks(_xticks[::2])

                if 'xlabel' in self.meta and jY==0:
                    ax.set_xlabel(self.meta['xlabel'])
                if 'ylabel' in self.meta and iX==0:
                    ax.set_ylabel(self.meta['ylabel'])
                
                if 'ylim' in self.meta:
                    ax.set_ylim(self.meta['ylim'])
                if 'xlim' in self.meta:
                    ax.set_xlim(self.meta['xlim'])

        
        # self.axs[CCDkey][Q].locator_params(nticks=4,axis='x')

    def plt_trimmer(self):
        
        
        for jY in range(self.Ncols):
            for iX in range(self.NSlices):            
            
                
                Ckey = 'C_%i%i' % (jY+1,iX+1)
                
                if jY>0:
                    plt.setp(self.axs[Ckey].get_xticklabels(),visible=False)
                
                
                if iX>0:
                    plt.setp(self.axs[Ckey].get_yticklabels(),visible=False)


                if self.meta['doRotateXLabels']:
                    for tick in self.axs[Ckey].get_xticklabels():
                        tick.set_rotation(45)
            
        if self.meta['doLegend']:
            plt.figlegend(self.handles, self.labels, loc='center right')
            
        if self.meta['doNiceXDate']:
            plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
            self.fig.autofmt_xdate()
            # plt.locator_params(nticks=4,axis='x',prune='both')

        plt.locator_params(axis='y', nbins=5, prune='both')
        
        #plt.locator_params(axis='y',prune='both')
        if not self.meta['doNiceXDate']:
            try: plt.locator_params(axis='x', nbins=4, prune='both')
            except: pass
            
        
        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(wspace=0.0)

        plt.margins(0.05)


        plt.suptitle(self.meta['suptitle'])
        # plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        if self.meta['doLegend']:
            plt.subplots_adjust(right=0.85)

        if self.meta['doColorbar']:
            #cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
            #plt.colorbar(cax=cbar_ax, mappable=self.mappables[0],orientation='vertical')
            self.fig.colorbar(self.mappables[0], ax=self.axsarr.flatten().tolist(), 
                              orientation='vertical', fraction=.1)

class FpaPlotYvsX(FpaPlot):
    
    def _ax_core_funct(self, ax, Cdict, key=''):
        
        ckwargs = self.corekwargs.copy()

        if key != '':
            
            xarr = Cdict['x'][key]
            yarr = Cdict['y'][key]
            
            label = st.replace(key, '_', '\_')
            kwargs=dict(label=label,marker='.',linestyle='')
            if key in ckwargs:
                kwargs.update(ckwargs[key])
            else:
                kwargs.update(ckwargs)
            handle = ax.plot(xarr, yarr, **kwargs)
        else:
            xarr = Cdict['x']
            yarr = Cdict['y']
            kwargs=dict(marker='.',linestyle='')
            kwargs.update(ckwargs)
            ax.plot(xarr, yarr, **kwargs)
            handle, label = None, None

        return handle, label


class FpaImgShow(FpaPlot):
    
    mappables = []
    
    def _ax_core_funct(self, ax, BCdict):
        
        internals = dict(origin='lower left')
        ckwargs = self.corekwargs.copy()
        internals.update(ckwargs)
        self.mappables.append(ax.imshow(BCdict['img'], **internals))
        handle, label = None, None
        return handle, label



class FpaHeatMap(BasicPlot):
    """ """
    
    def __init__(self, data, **kwargs):

        super(FpaHeatMap, self).__init__(**kwargs)
        
        self.figsize = (8,7)

        meta = dict(suptitle='',                     
                    doColorbar=False)
        meta.update(kwargs)
        
        self.NSlices = 6
        self.Ncols = 6 
        self.Quads = ['E', 'F', 'H', 'G']        
        
        self.data = copy.deepcopy(data)
        self.meta = dict()
        self.meta.update(meta)
        self.fig = None
        self.ax = None
        self.mappable = None
        
        self.corekwargs = dict()
        if 'corekwargs' in kwargs:
            self.corekwargs.update(kwargs['corekwargs']) 

    def init_fig(self):
        self._init_fig_and_axes()

    def _init_fig_and_axes(self):
        """ """
        plt.close('all')
        self.fig = plt.figure(figsize=self.figsize)        
        self.ax = self.fig.add_subplot(111)

    def _ax_core_funct(self,HM):
        """ """
        
        internals = dict(origin='lower left')
        ckwargs = self.corekwargs.copy()
        internals.update(ckwargs)
        
        self.mappable = self.ax.imshow(HM, **internals)
        
        return None
    
    def populate_axes(self):
        """ """
        
        NX = self.Ncols
        NY = self.NSlices
        NQ = len(self.Quads)
        
        HM = np.zeros((NX*NQ/2, NY*NQ/2))
        
        for jY in range(NY):
            for iX in range(NX):   
                Ckey = 'C_%i%i' % (jY+1,iX+1)
                
                CCdict = self.data[Ckey]
                
                ix = iX * NQ/2
                jy = jY * NQ/2
                
                if iX <=2:                
                    HM[jy+0,ix+0] = CCdict['E']
                    HM[jy+0,ix+1] = CCdict['F']
                    HM[jy+1,ix+1] = CCdict['G']
                    HM[jy+1,ix+0] = CCdict['H']
                else:
                    HM[jy+1,ix+1] = CCdict['E']
                    HM[jy+1,ix+0] = CCdict['F']
                    HM[jy+0,ix+0] = CCdict['G']
                    HM[jy+0,ix+1] = CCdict['H']
        
        
        self._ax_core_funct(HM)

        

    def plt_trimmer(self):
        
        NX = self.Ncols
        NY = self.NSlices
        NQ = len(self.Quads)
        
        for iX in range(NX):
            self.ax.axvline(x=iX*NQ/2+1.5,color='w')
            for jY in range(NY):
                self.ax.axhline(y=jY*NQ/2+1.5,color='w')
        
        self.ax.set_xticks(np.arange(6)*2+0.5)
        self.ax.set_xticklabels(['C1','C2','C3','C4','C5','C6'])
        
        self.ax.set_yticks(np.arange(6)*2+0.5)
        self.ax.set_yticklabels(['B1','B2','B3','B4','B5','B6'])                



        plt.suptitle(self.meta['suptitle'])
        # plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        if self.meta['doColorbar']:
            #cbar_ax = self.fig.add_axes([0.85, 0.15, 0.05, 0.7])
            #plt.colorbar(cax=cbar_ax, mappable=self.mappables[0],orientation='vertical')
            self.fig.colorbar(self.mappable, ax=self.ax, 
                              orientation='vertical', fraction=.1)    
    

class FpaFindingChart(FpaHeatMap):
    """ """
    
    def __init__(self, **kwargs):

        super(FpaFindingChart, self).__init__(data=None,**kwargs)
        

        meta = dict()
        meta.update(kwargs)
                
        self.meta.update(meta)        
        

    def init_fig(self):
        self._init_fig_and_axes()

    def _init_fig_and_axes(self):
        """ """
        plt.close('all')
        self.fig = plt.figure(figsize=self.figsize)        
        self.ax = self.fig.add_subplot(111)


    
    def populate_axes(self):
        """ """
        
        NX = self.Ncols
        NY = self.NSlices
        NQ = len(self.Quads)
        
        
        # LABELLING QUADRANTS
        
        for jY in range(NY):
            for iX in range(NX):   
                
                ix = iX * NQ/2
                jy = jY * NQ/2
                
                if iX <=2:   
                    items = [[jy+0,ix+0,'E'],
                             [jy+0,ix+1,'F'],
                             [jy+1,ix+1,'G'],
                             [jy+1,ix+0,'H']]
                else:   
                    items = [[jy+1,ix+1,'E'],
                             [jy+1,ix+0,'F'],
                             [jy+0,ix+0,'G'],
                             [jy+0,ix+1,'H']]
                
                for item in items:
                        self.ax.text(item[1],item[0],'%i%i%s' % (jY+1,iX+1,item[2]),
                                     fontsize=7,
                                     horizontalalignment='center',
                                     verticalalignment='center')    


        # LABELLING BLOCKS
        
        for jY in range(NY):
            
            for iX in [1,4]:
                
                ix = iX * NQ/2 + 0.5
                jy = jY * NQ/2 + 0.5
                
                Ckey = 'C_%i%i' % (jY+1,iX+1)
                #print(Ckey)
                blockname = fpamod.FPA_MAP[Ckey][0]
                
                blockID = '%s/FM%i' % (blockname,fpamod.BLOCK_SNs[blockname])
                
                self.ax.text(ix,jy,blockID,
                                     fontsize=24,
                                     alpha=0.2,
                                     horizontalalignment='center',
                                     verticalalignment='center')    

    def plt_trimmer(self):
        
        NY = self.Ncols
        NX = self.NSlices
        NQ = len(self.Quads)
        
        for iX in range(NX):
            self.ax.axvline(x=iX*NQ/2+.5,linestyle='--',linewidth=0.5,color='k')
            self.ax.axvline(x=iX*NQ/2+1.5,color='k')
            for jY in range(NY):
                self.ax.axhline(y=jY*NQ/2+.5,linestyle='--',linewidth=0.5,color='k')
                self.ax.axhline(y=jY*NQ/2+1.5,color='k')
        
        self.ax.set_xlim([-0.5,11.5])
        self.ax.set_ylim([-0.5,11.5])
        
        self.ax.set_xticks(np.arange(6)*2+0.5)
        self.ax.set_xticklabels(['CCD1','CCD2','CCD3','CCD3','CCD2','CCD1'])
        
        self.ax.set_yticks(np.arange(6)*2+0.5)
        self.ax.set_yticklabels(['B1','B2','B3','B4','B5','B6'])                


        plt.suptitle(self.meta['suptitle'])
        # plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        
    