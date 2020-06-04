"""

Bayesian CCD Spot Measurements
==============================


NEEDSREVISION AND UPGRADES

TODO:
    - single fluence multi-spot fitting (possible?? all spots have different fluences... even for same exp-time)
    - single spot multi-fluence fitting
    - multi-spot, multi-fluence fitting (too many parameters?)


Azzollini: (18th September 2017)

WARNING: This is a modified version of a script (BayesPSF_clinic.py, tagged on
8/Jun/2017), used in an investigation of the fiduciality of the Bayesian fits
of the opto-mechanic-ccd PSF, during RUN 2 of the EM1a CCD Characterisation at
MSSL.


Sami-Matias Niemi:

"Analyse laboratory CCD PSF measurements by forward modelling to the data.

The methods used here seem to work reasonably well when the spot has been well centred. If however, the
spot is e.g. 0.3 pixels off then estimating the amplitude of the Airy disc becomes rather difficult.
Unfortunately this affects the following CCD PSF estimates as well and can lead to models that are
rather far from the truth. Also, if when the width of the CCD PSF kernel becomes narrow, say 0.2, pixels
it is very difficult to recover. This most likely results from an inadequate sampling. In this case it might
be more appropriate to use "cross"-type kernel.

Because the amplitude can be very tricky to estimate, the version 1.4 (and later) implement a meta-parameter
called peak, which is the peak counts in the image. This is then converted to the amplitude by using the centroid
estimate. Because the centroids are fitted for each image, the amplitude can also vary in the joint fits. This
seems to improve the joint fitting constrains. Note however that this does couple the radius of the Airy disc
as well, because the amplitude estimate uses the x, y, and radius information as well.

One question to address is how the smoothing of the Airy disc is done. So far I have assumed that the Gaussian that
represents defocus should be centred at the same location as the Airy disc. However, if the displacement if the
Airy disc is large, then the defocus term will move the Airy disc to the direction of the displacement and make
it more asymmetric. Another option is to use scipy.ndimage.filters.gaussian_filter which simply applies Gaussian
smoothing to the input image. Based on the testing carried out this does not seem to make a huge difference. The
latter (smoothing) will lead to more or less the same CCD PSF estimates, albeit with slightly higher residuals.
We therefore adopt a Gaussian kernel that is centred with the Airy disc."


:requires: NumPy
:requires: SciPy
:requires: astropy
:requires: matplotlib
:requires: VISsim-Python
:requires: emcee
:requires: sklearn
:requires: multiprocessing


:author: R. Azzollini, Sami-Matias Niemi

"""
import matplotlib
# matplotlib.use('Agg')
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
import numpy as np
import emcee

import scipy.ndimage.measurements as m
from scipy import signal
from scipy.special import j1, jn_zeros
from astropy.modeling import models  # , fitting
import corner

from multiprocessing import Pool
from astropy.io import fits as fts
from pdb import set_trace as stop

import datetime

__author__ = 'R. Azzollini, Sami-Matias Niemi'


model_pars = dict(singauss=['peak', 'center_x', 'center_y', 'sigmax', 'sigmay'],
                  airyngauss=['peak', 'center_x', 'center_y', 'radius', 'focus'],
                  ccdgauss=['peak', 'center_x', 'center_y', 'radius', 'focus', 'sigmax', 'sigmay'],
                  ccdexp=['peak', 'center_x', 'center_y', 'radius', 'focus', 'sigmax', 'sigmay'],
                  ccdtriangle=['peak', 'center_x', 'center_y', 'radius', 'focus', 'weightx', 'weighty'])

model_priors = dict(singauss=[None, [7., 14.], [7., 14.], [0.1, 2.], [0.1, 2.]],
                    airyngauss=[None, [7., 14.], [7., 14.], [0.2, 2.], [0.2, 2.]],
                    ccdgauss=[None, [7., 14.], [7., 14.], [0.2, 2.], [0.2, 2.], [1.E-6, 0.6], [1.E-6, 0.6]],
                    ccdexp=[None, [7., 14.], [7., 14.], [0.2, 2.], [0.2, 2.], [1.E-6, 0.6], [1.E-6, 0.6]],
                    ccdtriangle=[None, [7., 14.], [7., 14.], [0.05, 2.], [0.05, 2.], [1.E-8, 1.E-2], [1.E-8, 1.E-2]])


def save_stamps_fits(outfits, stampdict, gof, modeltype, gain, burn, run):
    """ """

    hdu = fts.PrimaryHDU()
    hdu.header.add_history('PSF02_bayes.py at %s' %
                           (datetime.datetime.isoformat(datetime.datetime.now())))

    hdu.header.add_history('Extension 1: STAMP')
    hdu.header.add_history('Extension 2: MODEL')
    hdu.header.add_history('Extension 3: RESIDUAL')
    hdu.header.add_history('Extension 4: RESIDUALSQ')

    hdulist = fts.HDUList([hdu])

    hdustamp = fts.ImageHDU(stampdict['stamp'].transpose())
    hdustamp.header['EXTNAME'] = 'STAMP'
    hdulist.append(hdustamp)

    #hdukernel = fts.ImageHDU(stampdict['kernel'].transpose())
    #hdukernel.header['EXTNAME'] = 'KERNEL'
    # hdulist.append(hdukernel)

    hdumodel = fts.ImageHDU(stampdict['model'].transpose())
    hdumodel.header['EXTNAME'] = 'MODEL'
    hdulist.append(hdumodel)

    hdures = fts.ImageHDU(stampdict['residual'].transpose())
    hdures.header['EXTNAME'] = 'RESIDUAL'
    hdulist.append(hdures)

    hduressq = fts.ImageHDU(stampdict['residualSQ'].transpose())
    hduressq.header['EXTNAME'] = 'RESIDUALSQ'
    hdulist.append(hduressq)

    hdulist.writeto(outfits, clobber=True)


def log_posterior(theta, x, y, z, var, peakrange, spshape, modeltype):
    """
    Posterior probability: combines the prior and likelihood.
    """

    priors = model_priors[modeltype]
    lp = log_prior(theta, peakrange, priors)

    if not np.isfinite(lp):
        return -np.inf

    return lp + log_likelihood(theta, x, y, z, var, spshape, modeltype)


def log_prior(theta, peakrange, priors):
    """
    Generic model

    Priors, limit the values to a range but otherwise flat.
    """

    #priors = model_priors[model]

    areinrange = [(priors[i][0] < theta[i] < priors[i][1])
                  for i in range(1, len(theta))]
    disc = (peakrange[0] < theta[0] < peakrange[1]) * np.prod(areinrange)

    if disc:
        return 0.
    return -np.inf


def log_likelihood(theta, x, y, data, invar, spshape, modeltype):
    """
    Logarithm of the likelihood function.
    """

    model = generate_spot_model(theta, x, y, spshape, modeltype).flatten()
    # true for Gaussian errors
    #lnL = - 0.5 * np.sum((data - model)**2 / var)
    # Gary B. said that this should be from the model not data so recompute var (now contains rn**2)

    var = invar + model.copy()
    lnL = - 0.5 * np.sum((data - model)**2 / var)
    #lnL = - (np.size(var)*np.sum(np.log(var))) - (0.5 * np.sum((data - model)**2 / var))

    return lnL


def generate_singauss(theta, x, y, spshape, Full=False):
    """Single Gaussian Model"""

    # unpack the parameters
    peak, center_x, center_y, sigmax, sigmay = theta

    f = models.Gaussian2D(peak, center_x, center_y, sigmax, sigmay, 0.)

    model = f.evaluate(x, y, peak, center_x, center_y,
                       sigmax, sigmay, 0.).reshape(spshape)

    return model


def generate_airyngauss(theta, x, y, spshape, Full=False):
    """Airy Disk convolved with Gaussian de-focus kernel"""

    # unpack the parameters
    peak, center_x, center_y, radius, focus = theta

    # 1)Generate a model Airy disc
    #amplitude = _amplitudeFromPeak(peak, center_x, center_y, radius, x_0=int(spshape[0]/2.-0.5), y_0=int(spshape[1]/2.-0.5))
    amplitude = peak  # TESTS
    airy = models.AiryDisk2D(amplitude, center_x, center_y, radius)
    adata = airy.evaluate(x, y, amplitude, center_x,
                          center_y, radius).reshape(spshape)

    # 2)Apply Focus
    f = models.Gaussian2D(1., center_x, center_y, focus, focus, 0.)

    focusdata = f.evaluate(x, y, 1., center_x, center_y,
                           focus, focus, 0.).reshape(spshape)
    if focusdata.max() > 0.:
        focusdata /= focusdata.sum()  # normalization
        model = signal.convolve2d(adata, focusdata, mode='same')
    else:
        model = adata.copy()

    if not Full:
        return model
    else:
        return model, focusdata, adata


def generate_ccdgauss(theta, x, y, spshape, Full=False):
    """Airy Disk * gaussian de-focus kernel * gaussian ccd kernel"""

    # unpack the parameters
    peak, center_x, center_y, radius, focus, sigmax, sigmay = theta

    suboutputs = generate_airyngauss(theta[0:5], x, y, spshape, Full)

    if Full:
        model = suboutputs
    else:
        model, focusdata, adata = suboutputs

    # 3)Apply CCD diffusion, approximated with a Gaussian

    CCD = models.Gaussian2D(
        1., spshape[0] / 2. - 0.5, spshape[1] / 2. - 0.5, sigmax, sigmay, 0.)
    CCDdata = CCD.evaluate(
        x, y, 1., spshape[0] / 2. - 0.5, spshape[1] / 2. - 0.5, sigmax, sigmay, 0.).reshape(spshape)

    if CCDdata.max() > 0.:

        norm = 2. * np.pi * (sigmax**2. + sigmay**2.)
        CCDdata /= norm
        model = signal.convolve2d(model, CCDdata, mode='same')
    else:
        pass

    if not Full:
        return model
    else:
        return model, CCDdata, focusdata, adata


def generate_ccdexp(theta, x, y, spshape, Full=False):
    """Airy Disk * gaussian de-focus * exponential ccd kernel"""

    # unpack the parameters
    peak, center_x, center_y, radius, focus, sigmax, sigmay = theta

    suboutputs = generate_airyngauss(theta[0:5], x, y, spshape, Full)

    if Full:
        model = suboutputs
    else:
        model, focusdata, adata = suboutputs

    # 3)Apply CCD diffusion

    CCDdata = np.exp(-(np.abs((x - center_x)) / sigmax +
                       np.abs(y - center_y) / sigmay)).reshape(spshape)

    if CCDdata.max() > 0.:
        CCDdata /= CCDdata.sum()
        model = signal.convolve2d(model, CCDdata, mode='same')
    else:
        pass

    if not Full:
        return model
    else:
        return model, CCDdata, focusdata, adata


def generate_ccdtriangle(theta, x, y, spshape, Full=False):
    """Airy Disk * gaussian de-focus * 'pyramid' ccd kernel"""

    # unpack the parameters
    peak, center_x, center_y, radius, focus, weightx, weighty = theta

    suboutputs = generate_airyngauss(theta[0:5], x, y, spshape, Full)

    if Full:
        model = suboutputs
    else:
        model, focusdata, adata = suboutputs

    CCDdata = np.array([[0.0, weighty, 0.0],
                        [weightx, (1. - weighty - weighty -
                                   weightx - weightx), weightx],
                        [0.0, weighty, 0.0]])

    if CCDdata.max() > 0.:
        CCDdata /= CCDdata.sum()
        model = signal.convolve2d(model, CCDdata, mode='same')
    else:
        pass

    if not Full:
        return model
    else:
        return model, CCDdata, focusdata, adata


model_generators = dict(singauss=generate_singauss,
                        airyngauss=generate_airyngauss,
                        ccdgauss=generate_ccdgauss,
                        ccdexp=generate_ccdexp,
                        ccdtriangle=generate_ccdtriangle)


def generate_spot_model(theta, x, y, spshape, modeltype, Full=False):
    """ """
    generator = model_generators[modeltype]
    return generator(theta, x, y, spshape, Full)


def prepare_stamp(img, gain=3.5, size=10, spotx=100., spoty=100.):
    """ """

    oimg = img * gain

    # roughly the correct location - to avoid identifying e.g. cosmic rays
    data = oimg[spoty - (size * 3):spoty + (size * 3) + 1, spotx -
                (size * 3):spotx + (size * 3) + 1].copy()

    # maximum position within the cutout

    y, x = m.maximum_position(data)

    # spot and the peak pixel within the spot, this is also the CCD kernel position
    spot = data[y - size:y + size + 1, x - size:x + size + 1].copy()
    CCDy, CCDx = m.maximum_position(spot)
    print('CCD Kernel Position (within the postage stamp):', CCDx, CCDy)

    # bias estimate

    # bias = np.median(eimg[spoty-size: spoty+size, spotx-220:spotx-20]) #works for read o
    #rn = np.std(eimg[spoty-size: spoty+size, spotx-220:spotx-20])
    #bias = np.median(stats.sigma_clip(eimg,6))
    #rn = np.std(stats.sigma_clip(eimg,6))
    #rn = 4.5
    mask = np.zeros_like(oimg, dtype='bool')
    mask[y - 2 * size:y + 2 * size, x - 2 * size:x + 2 * size] = True
    bias = np.nanmedian(oimg[~mask])
    rn = np.std(oimg[~mask])

    #print 'Readnoise (e):', rn
    # if rn < 2. or rn > 6.:
    #    print 'NOTE: suspicious readout noise estimate...'
    print('ADC offset (e):', bias)

    # remove bias
    spot -= bias

    return spot, bias, rn


def forwardModel(
        spot,
        rn,
        stampkey='Unknown',
        modeltype='gauss',
        burn=500,
        run=700,
        cores=8,
        drill=False):
    """
    Forward models the spot data found from the input file. Can be used with simulated and real data.

    Notes:
    - emcee is run three times as it is important to have a good starting point for the final run.

    """
    print('\n\n\n')
    print('_' * 120)
    print('Processing: %s' % stampkey)
    # get data and convert to electrons

    parnames = model_pars[modeltype]
    priors = model_priors[modeltype]

    outputs = []
    results = dict()
    stampdict = dict()

    # save to file
    #fileIO.writeFITS(spot, stampkey+'_small.fits', int=False)
    # outputs.append(stampkey+'_small.fits')
    stampdict['stamp'] = spot.copy()

    # make a copy to generate error array
    fdata = spot.copy().flatten()
    # assume that uncertanties scale as sqrt of the values + readnoise
    #var = tmp.copy() + rn**2
    #var = np.zeros_like(fdata) + rn**2.
    # Gary B. said that actually this should be from the model or is biased,
    # so I only pass the readout noise part now

    # fit a simple model
    #print 'Least Squares Fitting...'
    #gauss = models.Gaussian2D(spot.max(), size, size, x_stddev=0.5, y_stddev=0.5)
    # gauss.theta.fixed = True  #fix angle
    #p_initG = gauss
    #fit_p = fitting.LevMarLSQFitter()
    #stopy, stopx = spot.shape
    #XG, YG = np.meshgrid(np.arange(0, stopx, 1), np.arange(0, stopy, 1))
    #pG = fit_p(p_initG, XG, YG, spot)
    #print 'Gaussian Model (variance): ',pG
    #modelG = pG(XG, YG)
    #var = modelG.flatten() + rn**2.

    #fileIO.writeFITS(model, stampkey+'BasicModel.fits', int=False)
    #fileIO.writeFITS(model - spot, stampkey+'BasicModelResidual.fits', int=False)

    # goodness of fit
    #gof = (1./(np.size(data) - 5.)) * np.sum((model.flatten() - data)**2 / var)
    #print 'GoF:', gof
    #print 'Done\n\n'

    # maximum value
    maxs = np.max(spot)
    peakrange = (0.1 * maxs, 2. * maxs)
    sums = np.sum(spot)

    print('Maximum Value:', maxs)
    print('Sum of the values:', sums)
    print('Peak Range:', peakrange)

    # MCMC based fitting
    print('Bayesian Model Fitting...')
    nwalkers = 1000

    # Initialize the sampler with the chosen specs.
    # Create the coordinates x and y
    x = np.arange(0, spot.shape[1])
    y = np.arange(0, spot.shape[0])
    # Put the coordinates in a mesh
    xx, yy = np.meshgrid(x, y)

    # Flatten the arrays
    xx = xx.flatten()
    yy = yy.flatten()

    print("Let's go Bayes...")

    ndim = len(parnames)

    # Choose an initial set of positions for the walkers - fairly large area not to bias the results
    p0 = np.zeros((nwalkers, ndim))
    #peak, center_x, center_y, radius, focus, width_x, width_y = theta
    p0[:, 0] = np.random.normal(
        maxs, maxs / 100., size=nwalkers)                 # peak value

    for ip in range(1, ndim):
        p0[:, ip] = np.random.uniform(priors[ip][0], priors[ip][1], nwalkers)

    if not drill:
        # stop()
        # initiate sampler

        if cores > 1:
            # A hack Dan gave me to not have ghost processes running as with threads keyword
            pool = Pool(cores)
            # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[xx,
            # yy, data, var, peakrange, spot.shape],
            sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior,
                                            args=[xx, yy, fdata, rn**2.,
                                                  peakrange, spot.shape, modeltype],
                                            pool=pool)
        else:
            sampler = emcee.EnsembleSampler(
                nwalkers, ndim, log_posterior, args=[
                    xx, yy, fdata, rn**2., peakrange, spot.shape, modeltype])

        # Run a burn-in and set new starting position
        print("Burning-in...")
        pos, prob, state = sampler.run_mcmc(p0, burn)
        maxprob_index = np.argmax(prob)
        params_fit = pos[maxprob_index]
        print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))
        print('Estimate:', params_fit)
        sampler.reset()

        print("Running MCMC...")
        pos, prob, state = sampler.run_mcmc(pos, run, rstate0=state)
        print("Mean acceptance fraction:", np.mean(sampler.acceptance_fraction))

        # Get the index with the highest probability
        maxprob_index = np.argmax(prob)

        # Get the best parameters and their respective errors and print best fits
        params_fit = pos[maxprob_index]
        errors_fit = [sampler.flatchain[:, i].std() for i in xrange(ndim)]
    else:
        params_fit = [maxs, -1., -1., 1.5, 0.6, 0.02, 0.03]
        errors_fit = np.zeros_like(params_fit)

    _printResults(params_fit, errors_fit)

    # Best fit model
    #peak, center_x, center_y, radius, focus, width_x, width_y = params_fit

    model = generate_spot_model(params_fit, xx, yy, spot.shape, modeltype)

    #amplitude = _amplitudeFromPeak(peak, center_x, center_y, radius, x_0=int(spot.shape[0]/2.-0.5), y_0=int(spot.shape[1]/2.-0.5))
    amplitude = params_fit[0]  # TESTS

    #stampdict['kernel'] = model*0.

    stampdict['model'] = model.copy()

    tmp = spot.copy()
    tmp[tmp + rn**2 < 0.] = 0.  # set highly negative values to zero
    var = tmp.copy() + rn**2

    # residuals
    stampdict['residual'] = (model - spot).copy()
    stampdict['residualSQ'] = ((model - spot)**2. / var).copy()

    # a simple goodness of fit
    gof = (1. / (np.size(spot) - ndim)) * \
        np.sum((model.flatten() - fdata)**2 / var.flatten())
    maxreldiff = np.max(np.abs(model.flatten() - fdata) / np.sqrt(var.flatten()))
    print('GoF:', gof, ' Maximum (abs(difference)/sigma) :', maxreldiff)
    if maxreldiff > 10 or gof > 4.:
        print('\nFIT UNLIKELY TO BE GOOD...\n')
    print('Amplitude estimate:', amplitude)

    # packing results

    # peak, center_x, center_y, radius, focus, width_x, width_y = params_fit

    for ip, parname in enumerate(parnames):
        results[parname] = params_fit[ip]
        results['e%s' % parname] = errors_fit[ip]

    results['GoF'] = gof

    if not drill:

        # plot
        samples = sampler.chain.reshape((-1, ndim))
        fig = corner.corner(samples, labels=parnames)
        fig.suptitle('%s, GoF=%.3f' % (modeltype, gof))
        fig.savefig(stampkey + '_Triangle.png')
        plt.close()

        outputs.append(stampkey + '_Triangle.png')

    outfits = stampkey + '_bayes.fits'
    save_stamps_fits(outfits, stampdict, gof, modeltype, 1., burn, run)

    outputs.append(outfits)

    if cores > 1:
        pool.close()

    return results, outputs


def save_stamps_fits(outfits, stampdict, gof, modeltype, gain, burn, run):
    """ """

    hdu = fts.PrimaryHDU()
    hdu.header.add_history('PSF02_bayes.py at %s' %
                           (datetime.datetime.isoformat(datetime.datetime.now())))

    hdu.header.add_history('Extension 1: STAMP')
    hdu.header.add_history('Extension 2: MODEL')
    hdu.header.add_history('Extension 3: RESIDUAL')
    hdu.header.add_history('Extension 4: RESIDUALSQ')

    hdulist = fts.HDUList([hdu])

    hdustamp = fts.ImageHDU(stampdict['stamp'].transpose())
    hdustamp.header['EXTNAME'] = 'STAMP'
    hdulist.append(hdustamp)

    #hdukernel = fts.ImageHDU(stampdict['kernel'].transpose())
    #hdukernel.header['EXTNAME'] = 'KERNEL'
    # hdulist.append(hdukernel)

    hdumodel = fts.ImageHDU(stampdict['model'].transpose())
    hdumodel.header['EXTNAME'] = 'MODEL'
    hdulist.append(hdumodel)

    hdures = fts.ImageHDU(stampdict['residual'].transpose())
    hdures.header['EXTNAME'] = 'RESIDUAL'
    hdulist.append(hdures)

    hduressq = fts.ImageHDU(stampdict['residualSQ'].transpose())
    hduressq.header['EXTNAME'] = 'RESIDUALSQ'
    hdulist.append(hduressq)

    hdulist.writeto(outfits, clobber=True)


def _printResults(best_params, errors):
    """
    Print basic results.
    """
    print(("=" * 60))
    print('Fitting with MCMC:')
    pars = ['peak', 'center_x', 'center_y',
            'radius', 'focus', 'width_x', 'width_y']
    print(('*' * 20 + ' Fitted parameters ' + '*' * 20))
    for name, value, sig in zip(pars, best_params, errors):
        print(("{:s} = {:e} +- {:e}" .format(name, value, sig)))
    print(("=" * 60))


def _printFWHM(sigma_x, sigma_y, sigma_xerr, sigma_yerr, req=10.8):
    """
    Print results and compare to the requirement at 800nm.
    """
    print(("=" * 60))
    print('FWHM (requirement %.1f microns):' % req)
    print(round(np.sqrt(_FWHMGauss(sigma_x) * _FWHMGauss(sigma_y)), 2), ' +/- ', \
        round(np.sqrt(_FWHMGauss(sigma_xerr) * _FWHMGauss(sigma_yerr)), 3), ' microns')
    print('x:', round(_FWHMGauss(sigma_x),
                      2), ' +/- ', round(_FWHMGauss(sigma_xerr), 3), ' microns')
    print('y:', round(_FWHMGauss(sigma_y),
                      2), ' +/- ', round(_FWHMGauss(sigma_yerr), 3), ' microns')
    print(("=" * 60))


def _FWHMGauss(sigma, pixel=12):
    """
    Returns the FWHM of a Gaussian with a given sigma.
    The returned values is in microns (pixel = 12microns).
    """
    return sigma * 2 * np.sqrt(2 * np.log(2)) * pixel


def _ellipticityFromGaussian(sigmax, sigmay):
    """
    Ellipticity
    """
    return np.abs((sigmax**2 - sigmay**2) / (sigmax**2 + sigmay**2))


def _ellipticityerr(sigmax, sigmay, sigmaxerr, sigmayerr):
    """
    Error on ellipticity.
    """
    e = _ellipticityFromGaussian(sigmax, sigmay)
    err = e * np.sqrt((sigmaxerr / e)**2 + (sigmayerr / e)**2)
    return err


def _R2FromGaussian(sigmax, sigmay, pixel=0.1):
    """
    R2.
    """
    return (sigmax * pixel)**2 + (sigmay * pixel)**2


def _R2err(sigmax, sigmay, sigmaxerr, sigmayerr):
    """
    Error on R2.
    """
    err = np.sqrt((2 * _R2FromGaussian(sigmax, sigmay))**2 * sigmaxerr**2 +
                  (2 * _R2FromGaussian(sigmax, sigmay))**2 * sigmayerr**2)
    return err


def _amplitudeFromPeak(peak, x, y, radius, x_0=10, y_0=10):
    """
    This function can be used to estimate an Airy disc amplitude from the peak pixel, centroid and radius.
    """
    rz = jn_zeros(1, 1)[0] / np.pi
    r = np.sqrt((x - x_0) ** 2 + (y - y_0) ** 2) / (radius / rz)
    if r == 0.:
        return peak
    rt = np.pi * r
    z = (2.0 * j1(rt) / rt)**2
    amp = peak / z
    return amp


def _peakFromTruth(theta, size=21):
    """
    Derive the peak value from the parameters used for simulations.
    """
    amplitude, center_x, center_y, radius, focus, width_x, width_y = theta
    x = np.arange(0, size)
    y = np.arange(0, size)
    x, y = np.meshgrid(x, y)
    airy = models.AiryDisk2D(amplitude, center_x, center_y, radius)
    adata = airy.evaluate(x, y, amplitude, center_x, center_y, radius)
    return adata.max()
