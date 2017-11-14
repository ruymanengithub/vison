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
    
    def plot(self,*args,**kwargs):
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
        
    
def test():
    
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
    
    data = dict(E=dict(x=xdict,y=yEdict),
                F=dict(x=xdict,y=yFdict),
                G=dict(x=xdict,y=yGdict),
                H=dict(x=xdict,y=yHdict))
    
    meta = dict(suptitle='Test',doLegend=True)
    ccd2dplot = CCD2DPlot(data,meta=meta)
    
    ccd2dplot.render(figname='')
    
    stop()
    

if __name__ == '__main__':
    test()
