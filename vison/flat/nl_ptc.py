# -*- coding: utf-8 -*-
"""


Module with tools used in NL analysis, from PTC meta-analysis.

Created on Wed Oct 14 16:40:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from scipy import stats
import copy
import os
from scipy import interpolate as interp
from matplotlib import pyplot as plt

from vison.support.files import cPickleRead
from . import PTC0Xaux
from scipy import optimize as opt
# END IMPORT

def get_mse_var(samples):
    """ 
    https://en.wikipedia.org/wiki/Mean_squared_error#Variance

    'In statistics, the mean squared error (MSE)[1][2] or mean squared deviation (MSD) of 
    an estimator (of a procedure for estimating an unobserved quantity) measures the 
    average of the squares of the errors—that is, the average squared difference between 
    the estimated values and the actual value. 

    The MSE is a measure of the quality of an estimator —it is always non-negative, and 
    values closer to zero are better.'
    
    MSE(S2_n-1) = 1/n (gamma_2 + 2n/(n-1)) sigma_4
     
    """
    
    if isinstance(samples, np.ma.MaskedArray):
        flsamples = samples.data[np.where(~samples.mask)]
    else:
        flsamples = samples.flatten()

    
    n = flsamples.size
    MSEvar = np.nan
    if n>3:
        sigma = np.nanstd(flsamples,ddof=1)
        mu4 = stats.moment(flsamples,moment=4)
        
        gamma2 = mu4 / sigma**4.-3. # excess kurtosis
    
        MSEvar = 1./n * (gamma2 + 2*n/(n-1) ) * sigma**4.
    
    return MSEvar


def f_extract_PTC(self, ccdobjcol, medcol, varcol, evarcol, binfactor=1):
    """ """

    # HARDWIRED VALUES
    wpx = self.window['wpx']
    hpx = self.window['hpx']

    indices = copy.deepcopy(self.dd.indices)

    nObs, nCCD, nQuad = indices.shape[0:3]

    Quads = indices.get_vals('Quad')
    CCDs = indices.get_vals('CCD')

    tile_coos = dict()
    for Quad in Quads:
        tile_coos[Quad] = self.ccdcalc.get_tile_coos(Quad, wpx, hpx)
    Nsectors = tile_coos[Quads[0]]['Nsamps']
    sectornames = np.arange(Nsectors)

    Sindices = copy.deepcopy(self.dd.indices)
    if 'Sector' not in Sindices.names:
        Sindices.append(core.vIndex('Sector', vals=sectornames))

    # Initializing new columns

    valini = 0.
    self.dd.initColumn(medcol, Sindices, dtype='float32', valini=valini)
    self.dd.initColumn(varcol, Sindices, dtype='float32', valini=valini)
    self.dd.initColumn(evarcol, Sindices, dtype='float32', valini=valini)

    # labels should be the same accross CCDs. PATCH.
    label = self.dd.mx['label'][:, 0].copy()
    ulabels = np.unique(label)
    ObsIDs = self.dd.mx['ObsID'][:].copy()

    # Pairing ObsIDs

    self.dd.initColumn(
        'ObsID_pair', self.dd.mx['ObsID'].indices, dtype='int64', valini=0)


    for ulabel in ulabels:
        six = np.where(label == ulabel)
        nsix = len(six[0])
        if nsix % 2 ==0:
            ixeven = np.arange(0, nsix, 2)
            ixodd = np.arange(1, nsix, 2)

            self.dd.mx['ObsID_pair'][six[0][ixeven]] = ObsIDs[six[0][ixodd]]

    if not self.drill:

        # if self.proc_histo['Masked']:
        #    estimators = dict(median=np.ma.median,std=np.ma.std)
        # else:
        #    estimators = dict(median=np.median,std=np.std)

        dpath = self.inputs['subpaths']['ccdpickles']

        misspairs = []

        for iObs in range(nObs):

            _ObsID_pair = self.dd.mx['ObsID_pair'][iObs]
            if _ObsID_pair == 0:
                continue
            iObs_pair = np.where(ObsIDs == _ObsID_pair)[0][0]
            

            flu_i = self.dd.mx['flu_med_img'][iObs, ...].mean()
            flu_p = self.dd.mx['flu_med_img'][iObs_pair, ...].mean()

            flus = np.array([flu_i, flu_p])

            if flus.std() / flus.mean() > 0.1:

                self.dd.mx[medcol][iObs, ...] = np.nan
                self.dd.mx[varcol][iObs, ...] = np.nan

                misspairs.append((self.dd.mx['ObsID'][iObs], self.dd.mx['ObsID'][iObs_pair]))

                continue

            for jCCD, CCDk in enumerate(CCDs):

                ccdobj_odd_f = os.path.join(
                    dpath, '%s.pick' % self.dd.mx[ccdobjcol][iObs, jCCD])
                ccdobj_eve_f = os.path.join(
                    dpath, '%s.pick' % self.dd.mx[ccdobjcol][iObs_pair, jCCD])

                if (not os.path.exists(ccdobj_odd_f)) or \
                    (not os.path.exists(ccdobj_eve_f)):
                    continue

                ccdobj_odd = copy.deepcopy(
                    cPickleRead(ccdobj_odd_f))
                ccdobj_eve = copy.deepcopy(
                    cPickleRead(ccdobj_eve_f))

                evedata = ccdobj_eve.extensions[-1].data.copy()

                # easy way to subtract one image from the other
                ccdobj_sub = copy.deepcopy(ccdobj_odd)
                ccdobj_sub.sub_bias(evedata, extension=-1)

                for kQ in range(nQuad):

                    Quad = Quads[kQ]

                    _tile_coos = tile_coos[Quad]

                    if binfactor >1:

                        ccdobj_odd = PTC0Xaux.CCDclone(ccdobj_odd)

                        _meds = ccdobj_odd.get_tiles_stats(
                            Quad, _tile_coos, 'median', extension=-1, binfactor=binfactor)

                        # IT'S A SUBTRACTION, SO WE HAVE TO DIVIDE BY 2 THE VARIANCE!

                        ccdobj_sub = PTC0Xaux.CCDclone(ccdobj_sub)

                        _vars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, 'std', extension=-1, binfactor=binfactor)**2. / 2.

                        _evars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, statkey=None, 
                            estimator=get_mse_var,
                            extension=-1, binfactor=binfactor) / 2.**2.

                    else:
                        _meds = ccdobj_odd.get_tiles_stats(
                            Quad, _tile_coos, 'median', extension=-1)

                        # IT'S A SUBTRACTION, SO WE HAVE TO DIVIDE BY 2 THE VARIANCE!
                        _vars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, 'std', extension=-1)**2. / 2.

                        _evars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, statkey=None, 
                            estimator=f_extract_mse_var,
                            extension=-1) / 2.**2.

                    self.dd.mx[medcol][iObs, jCCD, kQ, :] = _meds.copy()
                    self.dd.mx[varcol][iObs, jCCD, kQ, :] = _vars.copy()
                    self.dd.mx[evarcol][iObs, jCCD, kQ, :] = _evars.copy()
        
        if len(misspairs) > 0:
            if self.report is not None:
                self.report.add_Text(
                    'Pairs with unequal fluence skipped: %s' %
                    misspairs.__repr__())
            if self.log is not None:
                self.log.info('Pairs with unequal fluence skipped: %s' % \
                    misspairs.__repr__())

fscale = 1./2.E5

def fnon_lin_zero(x,theta,scaled=False):
    """non-linearity function: polynomial with zero intercept"""
        
    if scaled: xp = copy.deepcopy(x)
    else: xp = x * fscale
    
    npar = len(theta)
    fval = np.zeros_like(xp)
    for i in range(npar):
        fval += theta[i] * xp**(i+1)
    
    if scaled: return fval
    else: return fval / fscale

def fnon_lin_pol(x,theta,scaled=False):
    """non-linearity function: just a polynomial"""
        
    if scaled: xp = copy.deepcopy(x)
    else: xp = x * fscale
    
    npar = len(theta)
    fval = np.zeros_like(xp)
    for i in range(npar):
        fval += theta[i] * xp**i
    
    if scaled: return fval
    else: return fval / fscale

def fcorr_lin_zero(x,theta,scaled=False):
    """Returns linear mu, for input array of non-lin mu, 
    and non-lin model parameters theta."""

    if not scaled: xp = x * fscale
    else: xp = copy.deepcopy(x)

    xx = np.linspace(0.,1.3,100)
    yy = fnon_lin_zero(xx,theta,scaled=True)

    #print yy.min(),yy.max(),min(xp),max(xp)
    xpl = interp.interp1d(yy,xx,kind='cubic')(xp).copy()

    if scaled: return xpl
    else: return xpl / fscale

def fcorr_lin_pol(x,theta,scaled=False):
    """Returns linear mu, for input array of non-lin mu, 
    and non-lin model parameters theta."""

    if not scaled: xp = x * fscale
    else: xp = copy.deepcopy(x)

    xx = np.linspace(0.,1.3,100)
    yy = fnon_lin_pol(xx,theta,scaled=True)

    #print yy.min(),yy.max(),min(xp),max(xp)
    xpl = interp.interp1d(yy,xx,kind='cubic')(xp).copy()

    if scaled: return xpl
    else: return xpl / fscale

def fder_non_lin_zero(x,theta,scaled=False):
    """derivative of the non-linearity function"""
        
    if scaled: xp = copy(x)
    else: xp = x * fscale
    
    npar = len(theta)
    fval = np.zeros_like(xp)
    for i in range(npar):
        fval += (i+1)*theta[i] * xp**(i)
    
    if scaled: return fval
    else: return fval / fscale
    

def fder_non_lin_pol(x,theta,scaled=False):
    """derivative of the non-linearity function"""
        
    if scaled: xp = copy(x)
    else: xp = x * fscale
    
    npar = len(theta)
    fval = np.zeros_like(xp)
    for i in range(1,npar):
        fval += i*theta[i] * xp**(i-1.)
    
    if scaled: return fval
    else: return fval / fscale

def make_fitPTC_func(ron, binfactor=1,model='zero'):
    """Retrieves function that generates n-l variances, for a fixed value of RON,
    and binning factor."""

    if model == 'zero':
        _fcorr_lin = fcorr_lin_zero
        _fder_non_lin = fder_non_lin_zero
    elif model == 'pol':
        _fcorr_lin = fcorr_lin_pol
        _fder_non_lin = fder_non_lin_pol

    
    def func(mu_nle,*theta):
        """Retrieves the non-lin variances for a set of input mu_nle,
        given non-linearity parameters theta."""
        
        try:
            mu_le = _fcorr_lin(mu_nle,theta,scaled=False).copy()
        except Exception as inst:
            #print type(inst)
            #print(inst)
            #print(theta)
            return np.zeros_like(mu_nle)-np.inf
        
        corr_factor = _fder_non_lin(mu_le,theta,scaled=False)     
        #corr_factor = fnon_lin(mu_le,theta,scaled=False)/mu_le
        var_nle = corr_factor * mu_le  + ron**2 # Poisson + gaussian_ro
        #print '\n',var_nle,mu_le

        var_nle /= binfactor**2.
        
        return var_nle

    return func

def plot_NLcurve(params_fit, model='zero'):
    """ """

    x = np.linspace(1.,2.E5,500)

    params_lin = copy.deepcopy(params_fit)
    if model == 'zero':
        params_lin[1:]=0. # null non-linear terms
        _fnon_lin = fnon_lin_zero
    elif model == 'pol':
        params_lin[2:]=0. # null non-linear terms
        _fnon_lin = fnon_lin_pol

    yout = _fnon_lin(x,params_fit,scaled=False)
    youtlin = _fnon_lin(x,params_lin,scaled=False)

    fig = plt.figure(figsize=(12,7))
    ax1 = fig.add_subplot(121)
    ax1.plot(x,yout-youtlin,'b-',label='NL-curve')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles,labels,loc='best')
    ax1.set_xlabel('Input [e]')
    ax1.set_ylabel('Delta-Response [e]')
    ax1.set_title('abs. non-lin')
    ax1.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
    ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    ax2 = fig.add_subplot(122)
    ax2.plot(x,(yout/youtlin-1.)*100.,'b-',label='NL-curve')
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles,labels,loc='best')
    ax2.set_xlabel('Input [e]')
    ax2.set_ylabel('Response [\%]')
    ax2.set_title('rel. non-lin')
    ax2.ticklabel_format(style='sci',axis='x',scilimits=(0,0))

    plt.tight_layout()
    plt.show()
    plt.close('all')

def forward_PTC_LM(indata, npol=6):
    """ """
    from matplotlib import pyplot as plt

    gain = indata['gain']
    ron = indata['ron']
    binfactor = indata['binfactor']
    mu_nle = indata['mu_nle']
    var_nle = indata['var_nle']
    evar_nle = indata['evar_nle']

    mu_nle *= gain
    var_nle *= gain**2.
    evar_nle *= gain**2.

    ixsort = mu_nle.argsort() # needed? should not harm, either

    mu_nle = mu_nle[ixsort]
    var_nle = var_nle[ixsort]
    evar_nle = evar_nle[ixsort]

    model = 'pol'


    doPlot1 = True
    if doPlot1:

        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        #ax.errorbar(mu_nle,var_nle*binfactor**2./mu_nle,
        #    yerr=evar_nle*binfactor**2./mu_nle,fmt='-.')
        ax.plot(mu_nle,var_nle*binfactor**2./mu_nle,'b.')
        ax.axhline(1.,ls='--',c='k')
        ax.set_xlabel(r'$\mu_{NL,e}$')
        ax.set_ylabel(r'$var_{NL,e}/\mu_{NL,e}$')
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-2,2))
        plt.tight_layout()
        plt.show()
        plt.close()

    stop()
    p0 = np.zeros(npol,dtype='float32')
    if model == 'zero':
        p0[0] = 1. # initial conditions: perfectly linear response
    elif model == 'pol':
        p0[1] = 1.

    fitfunc = make_fitPTC_func(ron, binfactor, model)

    popt,pcov = opt.curve_fit(fitfunc,mu_nle,var_nle,
        p0=p0,
        sigma=evar_nle,
        absolute_sigma=True,
        method='lm')
    
    #Get the best parameters and their respective errors and print best fits
    params_fit = popt
    errors_fit = np.sqrt(np.diagonal(pcov))
    print('Best Pars: ', params_fit)
    print('Error Pars: ', errors_fit)
    
    varle = fitfunc(mu_nle, *p0) # linear response
    varnle_best = fitfunc(mu_nle,*popt) # nl variances from mu_nle, given best fit model pars

    doPlot2 = True

    if doPlot2:

        fig2 = plt.figure()
        ax1 = fig2.add_subplot(111)
        ax1.plot(mu_nle, var_nle, 'k.', label='data')
        ax1.plot(mu_nle, varnle_best, 'r--', label='best-fit')
        ax1.plot(mu_nle, varle, 'g:', label='linear model')
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles, labels, loc='best')
        ax1.set_xlabel('mu\\_nle, [e]')
        ax1.set_ylabel('var [e2]')
        ax1.set_title('Fitting the var-mu relation')
        plt.show()
    plt.close()


    doPlotNL = True

    if doPlotNL:

        plot_NLcurve(params_fit, model)


    stop()
