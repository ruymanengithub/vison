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


class Fig(object):
    
    def __init__(self,filename=''):
        
        self.filename = filename
        self.texfraction = 1.0
        self.caption = ''
        
    def configure(self,**kwargs):
        """ """
        raise NotImplementedError('Subclass implements abstract method')
        
    
    def plot(self,figname=''):
        """ """
        raise NotImplementedError('Subclass implements abstract method')
    

class BasicPlot(object):
    
    def __init__(self,**kwargs):
        
        self.fig = None
        self.figsize = (9,9)
        if 'fig' in kwargs: self.fig = kwargs['fig']
        self.axs = []
     
    def axmethod(self):
        raise NotImplementedError("Subclass must implement abstract method")
    def plt_trimmer(self):
        raise NotImplementedError("Subclass must implement abstract method")
    
    def render(self,figname=''):
        
        
        self.fig = plt.figure(figsize=self.figsize)
        
        self.axmethod()
        
        self.plt_trimmer()
        
        if figname == '':
            plt.show()
        else:
            plt.savefig(figname)
        
        plt.close()
            

class CCD2DPlot(BasicPlot):
    
    
    def __init__(self,data,**kwargs):
        
        super(CCD2DPlot,self).__init__(**kwargs)
        
        defaults = dict(suptitle='',qtitles=dict(E='E',F='F',G='G',H='H'),
                        doLegend=False,
                        doNiceXDate=False)
        
        if 'meta' in kwargs: meta = kwargs['meta']
        else: meta = dict()
        
        self.figsize=(8,8)
        self.Quads = ['E','F','H','G']
        self.data = data
        self.meta = dict()
        self.meta.update(defaults)
        self.meta.update(meta)
    
        self.handles = []
        self.labels = []
    
    
    def axmethod(self):
        
        qtitles = self.meta['qtitles']
        
        self.axs = []
        for iQ, Q in enumerate(self.Quads):
            
            self.axs.append(self.fig.add_subplot(2,2,iQ+1))
            
            try: 
                xkeys = self.data[Q]['x'].keys()
            except AttributeError:
                xkeys = None
                
            if xkeys is not None:
                ykeys = self.data[Q]['y'].keys()
                isconsistent = np.all([xkeys[i] == ykeys[i] for i in range(len(xkeys))])
                assert (len(xkeys) == len(ykeys)) and isconsistent
                
                for key in xkeys:
                    xarr = self.data[Q]['x'][key]
                    yarr = np.ones_like(self.data[Q]['y'][key])
                    
                    handle = self.axs[-1].plot(xarr,yarr,label=key)
                    if iQ==0:
                        self.handles += handle
                        self.labels.append(key)
            else:
                xarr = self.data[Q]['x']
                yarr = self.data[Q]['y']
                self.axs[-1].plot(xarr,yarr)
            
            self.axs[-1].set_title(qtitles[Q])
        
    def plt_trimmer(self):
        
        if self.meta['doLegend']:
            plt.figlegend(self.handles,self.labels,loc='center right')
        
        
        if self.meta['doNiceXDate']:
            
            plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))   
            self.fig.autofmt_xdate()
        
        
        plt.suptitle(self.meta['suptitle'])
        plt.tight_layout()
        plt.subplots_adjust(top=0.90)
        if self.meta['doLegend']:
            plt.subplots_adjust(right=0.85)

class Beam2DPlot(BasicPlot):
    
    
    def __init__(self,data,**kwargs):
        
        super(Beam2DPlot,self).__init__(**kwargs)
        
        defaults = dict(suptitle='',ccdtitles=dict(CCD1='CCD1',CCD2='CCD2',CCD3='CCD3'),
                        doLegend=False,
                        doNiceXDate=False)
        
        if 'meta' in kwargs: meta = kwargs['meta']
        else: meta = dict()
        
        self.figsize=(15,6)
        self.Quads = ['E','F','H','G']
        self.CCDs = [1,2,3]
        self.data = data
        self.meta = dict()
        self.meta.update(defaults)
        self.meta.update(meta)
    
        self.handles = []
        self.labels = []
    
    
    def axmethod(self):
        """ 
        
        TODO:
            3 CCDs with 4 Quadrants each
            share x and y axes within CCDs, but with gaps between CCDs
            allow for optional legend
            allow for optional "nice" date x-axis
        
        
        """
        #qtitles = self.meta['qtitles']
        
        self.axs = dict()
        
        gss = dict()
        
        gss['CCD1'] = plt.GridSpec(2,6,hspace=0,wspace=0,left=0.01,right=0.30,top=0.85,bottom=0.15)
        gss['CCD2']  = plt.GridSpec(2,6,hspace=0,wspace=0,left=0.31,right=0.6,top=0.85,bottom=0.15)
        gss['CCD3']  = plt.GridSpec(2,6,hspace=0,wspace=0,left=0.61,right=0.90,top=0.85,bottom=0.15)
        
        vQs = (['E','F'],['G','H'])
        hQs = (['E','H'],['F','G'])
        
        for iCCD,CCD in enumerate(self.CCDs):
            
            CCDkey = 'CCD%i' % CCD
            
            self.axs[CCDkey] = dict()
            
            
            for iQ, Q in enumerate(self.Quads):
                
                ip = [ix for ix in range(2) if Q in vQs[ix]][0]
                jp = [ix for ix in range(2) if Q in hQs[ix]][0]
                
                if (Q == 'E') and (CCD == 1):
                    shareax = None
                else:
                    shareax = self.axs['CCD1']['E']
                    
                
                self.axs[CCDkey][Q] = self.fig.add_subplot(gss[CCDkey][ip,jp],sharex=shareax,sharey=shareax)
                
                try: 
                    xkeys = self.data[CCDkey][Q]['x'].keys()
                except AttributeError:
                    xkeys = None
                    
                if xkeys is not None:
                    ykeys = self.data[CCDkey][Q]['y'].keys()
                    isconsistent = np.all([xkeys[i] == ykeys[i] for i in range(len(xkeys))])
                    assert (len(xkeys) == len(ykeys)) and isconsistent
                    
                    for key in xkeys:
                        xarr = self.data[CCDkey][Q]['x'][key]
                        yarr = np.ones_like(self.data[CCDkey][Q]['y'][key])
                        
                        handle = self.axs[CCDkey][Q].plot(xarr,yarr,label=key)
                        if Q=='E' and CCD==1:
                            self.handles += handle
                            self.labels.append(key)
                else:
                    xarr = self.data[CCDkey][Q]['x']
                    yarr = self.data[CCDkey][Q]['y']
                    self.axs[CCDkey][Q].plot(xarr,yarr)
                
                if Q in ['E','H']:
                    self.axs[CCDkey][Q].text(0.05,0.9,Q,horizontalalignment='left',
                            transform=self.axs[CCDkey][Q].transAxes)
                else:
                    self.axs[CCDkey][Q].text(0.90,0.9,Q,horizontalalignment='left',
                            transform=self.axs[CCDkey][Q].transAxes)
                
                if Q == 'E': self.axs[CCDkey][Q].set_title(CCDkey,x=1)
                
                
                #self.axs[-1].set_title(qtitles[Q])
        
    def plt_trimmer(self):
        
        
        for CCD in self.CCDs:
            for Q in ['E','F']:
                plt.setp(self.axs['CCD%i' % CCD][Q].get_xticklabels(),visible=False)
            for Q in ['F','G']:
                plt.setp(self.axs['CCD%i' % CCD][Q].get_yticklabels(),visible=False)
        
        if self.meta['doLegend']:
            plt.figlegend(self.handles,self.labels,loc='center right')
        
        
        if self.meta['doNiceXDate']:
            
            plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))   
            self.fig.autofmt_xdate()
        
        
        plt.locator_params(axis='y',nticks=6,prune='both')
        plt.locator_params(axis='x',nticks=4,prune='both')
        
        plt.suptitle(self.meta['suptitle'])
        #plt.tight_layout()
        #plt.subplots_adjust(top=0.90)
        #if self.meta['doLegend']:
        #    plt.subplots_adjust(right=0.85)


    
def testBeam2DPlot():
    
    x = np.arange(10,dtype='float32')
    yE = x*2.
    yF = x*3.
    yG = x*4.
    yH = x*5.
    
    xdict = OrderedDict(foo=x,bar=x)
    yEdict = OrderedDict(foo=yE,bar=yE+5.)
    yFdict = OrderedDict(foo=yF,bar=yF+5.)
    yGdict = OrderedDict(foo=yG,bar=yG+5.)
    yHdict = OrderedDict(foo=yH,bar=yH+5.)
    
    ccddict = dict(E=dict(x=xdict,y=yEdict),
                F=dict(x=xdict,y=yFdict),
                G=dict(x=xdict,y=yGdict),
                H=dict(x=xdict,y=yHdict))
    
    data = dict(CCD1=copy.deepcopy(ccddict),
                CCD2=copy.deepcopy(ccddict),
                CCD3=copy.deepcopy(ccddict))
    
    meta = dict(suptitle='Test',doLegend=True)
    beam2dplot = Beam2DPlot(data,meta=meta)
    
    beam2dplot.render(figname='')
    
    stop()
    

if __name__ == '__main__':
    #testCCD2DPlot()
    
    testBeam2DPlot()
    
