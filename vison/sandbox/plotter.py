#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

vison pipeline: Classes to do plots.

Created on Mon Nov 13 17:54:08 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""
# IMPORT STUFF
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


class PData(object):
    """ """
    
    

class BasicPlot(object):
    
    def __init__(self,*args,**kwargs):
        
        self.fig = None
        self.axs = []
        self.axmethod = None
        
    
    def render(self,figname=''):
        
        
        self.fig = plt.figure(figsize=self.figsize)
        
        self.axs = self.axmethod()
        
        self.plt_trimmer(plt)
        
        if figname == '':
            plt.show()
        else:
            plt.savefig(figname)
        
        plt.close()
            

    
def test():
    pass
    
    

if __name__ == '__main__':
    test()
