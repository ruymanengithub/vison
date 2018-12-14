#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

ROE-TAB Calibration
Module to fit the Non-Linearity waveform using Levenberg-Marquardt

Created on Wed Mar 21 14:54:00 2018

:author: raf

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
#from multiprocessing import Pool
#import emcee
#import corner
import scipy.optimize as opt

from pylab import plot,show
import matplotlib
#matplotlib.use('Agg')
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

import Wave_bayes as Wbayes

# END IMPORT


def f_fit(X,*theta):
    """ """
    N = len(X)
    return Wbayes.waveform_generator(theta,N)



def fitWaveformNLSQ(Varr, Vrange, pixTx, Nlevels,wpix=None,doPlot=False):
    """
    'Uses non-linear least-squares to fit a function, f, to data'.
    
    """
    var = 0.05**2.
    
    # theta: phase, ground, wref, wpix, a_i {i=1,N}
    parnames = ['phi','gnd','wref','wpix']
    for i in range(1,Nlevels+1): parnames += ['a%i' % i]
    
    bounds_listoftuples = Wbayes.get_priors(Vrange,pixTx,Nlevels)
    
    lo_bounds = [item[0] for item in bounds_listoftuples]
    hi_bounds = [item[1] for item in bounds_listoftuples]
    bounds = (lo_bounds,hi_bounds)
    
    if wpix is not None:
        ixwpix = parnames.index('wpix')
        bounds[0][ixwpix] = wpix*0.999
        bounds[1][ixwpix] = wpix*1.001
    
    
    results = dict()
    
    N = len(Varr)
    
    sigmarr = np.ones(N,dtype='float32') * var**0.5
    
    ndim = len(lo_bounds)
    p0 = np.zeros(ndim)
    for ip in range(ndim):
        try: p0[ip] = np.random.uniform(lo_bounds[ip],hi_bounds[ip],1)
        except: stop()
        
    params_fit,pcov = opt.curve_fit(f_fit, np.arange(N,dtype='float32'),Varr,bounds=bounds,
                        p0=p0,method='dogbox',sigma=sigmarr,
                        max_nfev=1E8)
        
    model = f_fit(np.arange(N,dtype='float32'),*params_fit)
    
    if doPlot:
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(Varr)
        ax.plot(model,'k--')
        plt.show()
        plt.close()
    
    
    # packing results
    
    # peak, center_x, center_y, radius, focus, width_x, width_y = params_fit
    
    for ip,parname in enumerate(parnames):
        results[parname] = params_fit[ip]
    
    
    return results


def test0():
    """ """
    
    import ROE_LinCalib as RLC
    import os
    from astropy.io import ascii
    
    
    Nlevels = 17
    pixT = 14.245E-6 # s
    SampInter = 100.E-9
    datapath = 'CH1_Top_E_F_Pix_bounce/With_17_levels/'
    datafile = 'LinCalib_FQM_ROE_15Mar18.txt'
    wpix = None
    
    pixTx = int(np.round(pixT / SampInter))
    toi = pixT/3.
    #SampInter = 100.E-9 # s
    
    filt_kernel = int(np.round(toi/SampInter)/4.)
    
    
    indata = ascii.read(datafile)
    
    WFList = indata['WAVEFORM'].data.copy()
    
    WFf = os.path.join(datapath,WFList[0])
    
    timex,rV = RLC.load_WF(WFf,chkNsamp=1.E5,chkSampInter=SampInter)
                
    fV = RLC.filter_Voltage(rV,filt_kernel)
    
    Vrange = (fV.min(),fV.max())
    
    res = fitWaveformNLSQ(fV, Vrange, pixTx, Nlevels, wpix=wpix,doPlot=True)
    
    stop()
    
    
        
    
if __name__ == '__main__':
    
    test0()
