#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in PTC analysis.

Created on Thu Sep 14 16:29:36 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

#IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
#END IMPORT

def fitPTC(means,var): #,doplot=False,savefig=''):
    """ """
    
    
    order = np.argsort(means)
    means = np.array(means)[order]
    var = np.array(var)[order]   
    
    ixmaxvar = np.argmax(var)
    maxvar = var[ixmaxvar]
    
    ixsel = np.where((means < means[ixmaxvar]) & (var < maxvar*0.95))
    fvar = var[ixsel]
    fmeans = means[ixsel]
    
    
    res = np.polyfit(fmeans,fvar,2,full=False,cov=True)
    
    p = res[0]
    V = res[1]
    ep = np.sqrt(np.diag(V))
    
    g = 1./p[1]
    cuadterm = p[0]
    rn = p[2]
    
    badresult=False
    dynrange = fmeans.max()/fmeans.min()
    if (dynrange < 3) or fmeans.max()<2.E4:
        badresult=True
    
    fitresults = dict(fit=p,efit=ep,gain=g,cuadterm=cuadterm,rn=rn,badresult=badresult)
    
#    pmeans = np.linspace(0,fmeans.max(),100)
#    pvar  = np.polyval(p,pmeans)
#    
#    if doplot: 
#        fig = plt.figure(1)
#        ax = fig.add_subplot(111)
#        ax.plot(means,var,'ro')
#        ax.plot(fmeans,fvar,'bo')
#        ax.plot(pmeans,pvar,'b--')
#        ax.set_xlabel('mean [ADU]')
#        ax.set_ylabel('var [ADU^2]')
#        ax.set_title(savefig)
#        ax.set_xlim((0,50000.))
#        ax.set_ylim((0.,20000.))
#        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
#             ax.get_xticklabels() + ax.get_yticklabels()):
#            item.set_fontsize(16)
#        plt.tight_layout()
#        if savefig != '':
#            plt.savefig(savefig)
#        else:
#            plt.show()
#        plt.close('all')

    return fitresults


