#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:32:19 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from scipy import optimize as opt
import copy

from pylab import plot,show
import matplotlib
#matplotlib.use('pdf')
#matplotlib.rc('text', usetex=True)
matplotlib.rcParams['font.size'] = 17
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('axes', linewidth=1.1)
matplotlib.rcParams['legend.fontsize'] = 11
matplotlib.rcParams['legend.handlelength'] = 3
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['image.interpolation'] = 'none'
import matplotlib.pyplot as plt

# END IMPORT



def myplot_NL(datasets,xlabel='',ylabel='',ylim=[],title='',figname=''):
    """ """
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for dkey in datasets.keys():
        dataset = datasets[dkey]
        x = dataset['x']
        y = dataset['y']
        try: ls = dataset['ls']
        except: ls = ''
        try: color = dataset['color']
        except: color = None
        try: marker = dataset['marker']
        except: marker = None
        
        ax.plot(x,y,color=color,ls=ls,marker=marker,label=dkey)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if len(ylim)>0:
        ax.set_ylim(ylim)
    
    handles,labels = ax.get_legend_handles_labels()
    
    ax.legend(handles,labels,loc='best')
    
    ax.set_title(title)
    
    plt.tight_layout()
    
    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    plt.close()
    


def f_pol_byorigin(x,*theta):
    """ """
    res = np.zeros_like(x)
    for i in range(len(theta)):
        try: res += theta[i]*x**(i+1)
        except: stop()
    
    return res

def fit_pol_byorigin(x,y,deg=2,sigma=None):
    """ """
    
    if sigma is not None:
        assert len(y) == len(sigma)
    
    p0 = np.zeros(deg,dtype='float32')
    ixnoinf = np.where(~np.isinf(x) & ~np.isnan(x) & ~np.isinf(y) & ~np.isnan(y))
    

    pfit,_ = opt.curve_fit(f_pol_byorigin, x[ixnoinf],y[ixnoinf],
                        p0=p0,method='lm',sigma=sigma[ixnoinf],absolute_sigma=False)
    
    return pfit
    
def find_NL_pol(x,y,deg=2,sigma=None,Full=False, debug=False):
    """ """
    
    pfitLin = np.polyfit(x, y, deg=1, w=sigma)
    yLin = np.polyval(pfitLin,x)
    pfitLin_nozero = copy.deepcopy(pfitLin)
    pfitLin_nozero[-1] = 0.
    yLin_nozero = np.polyval(pfitLin_nozero,x)

    NLin = (y-yLin)/yLin_nozero
    
    pfitNLin = np.polyfit(x,NLin,deg=deg,w=sigma)
    
    if debug:
        stop()
    
    if Full:
        return pfitNLin,NLin
    else:
        return pfitNLin
    
