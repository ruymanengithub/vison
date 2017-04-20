"""
Measuring a shape of an object
==============================

Simple class to measure quadrupole moments and ellipticity of an object.

:requires: NumPy
:requres: PyFITS

:author: Sami-Matias Niemi, Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""
import os, datetime, unittest
import numpy as np
import pyfits as pf

from pdb import set_trace as stop
#from astropy.io import fits as fts


class Shapemeter():
    """
    Provides methods to measure the shape of an object.

    :param data: stamp to be analysed.
    :type data: np.ndarray
    :param log: logger
    :type log: instance
    :param kwargs: additional keyword arguments
    :type kwargs: dict

    Settings dictionary contains all parameter values needed.
    """

    def __init__(self, data, log=None, **kwargs):
        """
        :param data: stamp to be analysed.
        :type data: ndarray
        :param log: logger
        :type log: instance
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        Settings dictionary contains all parameter values needed.
        """
        self.data = data.copy()
        self.log = log

        NY, NX = self.data.shape
        self.NX = NX
        self.NY = NY

        self.shsettings = dict(iterations=4,
                             sampling=1.0,
                             platescale=120.0,
                             pixelSize=12.0,
                             sigma=0.75,
                             weighted=True,
                             conservePeak=True,
                             debug=False,
                             fixedPosition=False,
                             fixedX=None,
                             fixedY=None)
       
        self.shsettings.update(kwargs)
        

        for key, value in self.shsettings.iteritems():
            if self.log is not None: self.log.info('%s = %s' % (key, value))


    def quadrupoles(self):
        """
        Derive quadrupole moments and ellipticity from the input image.

        :param img: input image data
        :type img: ndarray

        :return: quadrupoles, centroid, and ellipticity (also the projected components e1, e2)
        :rtype: dict
        """
        if self.log is not None: self.log.info('Deriving quadrupole moments')
        image = self.data.copy()

        #normalization factor
        imsum = float(np.sum(image))

        #generate a mesh coordinate grid
        sizeY, sizeX = image.shape
        Xvector = np.arange(0, sizeX)
        Yvector = np.arange(0, sizeY)
        Xmesh, Ymesh = np.meshgrid(Xvector, Yvector)

        # No centroid given, take from data and weighting with input image
        Xcentre = np.sum(Xmesh.copy() * image.copy()) / imsum
        Ycentre = np.sum(Ymesh.copy() * image.copy()) / imsum

        #coordinate array
        Xarray = Xcentre * np.ones([sizeY, sizeX])
        Yarray = Ycentre * np.ones([sizeY, sizeX])

        #centroided positions
        Xpos = Xmesh - Xarray
        Ypos = Ymesh - Yarray

        #squared and cross term
        Xpos2 = Xpos * Xpos
        Ypos2 = Ypos * Ypos
        XYpos = Ypos * Xpos

        #integrand
        Qyyint = Ypos2 * image.copy()
        Qxxint = Xpos2 * image.copy()
        Qxyint = XYpos * image.copy()

        #sum over and normalize to get the quadrupole moments
        Qyy = np.sum(Qyyint) / imsum
        Qxx = np.sum(Qxxint) / imsum
        Qxy = np.sum(Qxyint) / imsum

        if self.log is not None: self.log.info('(Qxx, Qyy, Qxy) = (%f, %f, %f)' % (Qxx, Qyy, Qxy))

        #derive projections and ellipticity
        denom = Qxx + Qyy
        e1 = (Qxx - Qyy) / denom
        e2 = 2. * Qxy / denom
        ellipticity = np.sqrt(e1*e1 + e2*e2)

        #also a and b
        a = np.sqrt(.5 * (Qxx + Qyy + np.sqrt((Qxx - Qyy)**2 + 4.*Qxy*Qxy)))
        b = np.sqrt(.5 * (Qxx + Qyy - np.sqrt((Qxx - Qyy)**2 + 4.*Qxy*Qxy)))

        #check that ellipticity is reasonable
        if ellipticity > 1.0:
            if self.log is not None: self.log.error('Ellipticity greater than 1 derived, will set it to unity!')
            ellipticity = 1.0

        if self.log is not None: self.log.info('Centroiding (x, y) = (%f, %f) and ellipticity = %.4f (%.4f, %.4f)' %
                      (Ycentre+1, Xcentre+1, ellipticity, e1, e2))

        res = dict(ellipticity=ellipticity, e1=e1, e2=e2, Qxx=Qxx, Qyy=Qyy, Qxy=Qxy,
                   centreY=Ycentre, centreX=Xcentre,
                   a=a, b=b)
        return res


    def circular2DGaussian(self, x, y, sigma):
        """
        Create a circular symmetric Gaussian centered on x, y.

        :param x: x coordinate of the centre
        :type x: float
        :param y: y coordinate of the centre
        :type y: float
        :param sigma: standard deviation of the Gaussian, note that sigma_x = sigma_y = sigma
        :type sigma: float

        :return: circular Gaussian 2D profile and x and y mesh grid
        :rtype: dict
        """
        if self.log is not None: self.log.info('Creating a circular symmetric 2D Gaussian with sigma=%.3f centered on (x, y) = (%f, %f)' % (sigma, x, y))

        #x and y coordinate vectors
        Gyvect = np.arange(1, self.NY + 1)
        Gxvect = np.arange(1, self.NX + 1)

        #meshgrid
        Gxmesh, Gymesh = np.meshgrid(Gxvect, Gyvect)

        #normalizers
        sigmax = 1. / (2. * sigma**2)
        sigmay = sigmax #same sigma in both directions, thus same normalizer

        #gaussian
        exponent = (sigmax * (Gxmesh - x)**2 + sigmay * (Gymesh - y)**2)
        try: Gaussian = np.exp(-exponent) / (2. * np.pi * sigma*sigma)
        except: stop()

        if self.shsettings['conservePeak']:
            #normalize to unity
            Gaussian /= np.max(Gaussian)

        output = dict(GaussianXmesh=Gxmesh, GaussianYmesh=Gymesh, Gaussian=Gaussian)

        return output


    def ellip2DGaussian(self, x, y, sigmax, sigmay):
        """
        Create a two-dimensional Gaussian centered on x, y.

        :param x: x coordinate of the centre
        :type x: float
        :param y: y coordinate of the centre
        :type y: float
        :param sigmax: standard deviation of the Gaussian in x-direction
        :type sigmax: float
        :param sigmay: standard deviation of the Gaussian in y-direction
        :type sigmay: float


        :return: circular Gaussian 2D profile and x and y mesh grid
        :rtype: dict
        """
        if self.log is not None: self.log.info('Creating a 2D Gaussian with sigmax=%.3f and sigmay=%.3f centered on (x, y) = (%f, %f)' %
                      (sigmax, sigmay, x, y))

        #x and y coordinate vectors
        Gyvect = np.arange(1, self.NY + 1)
        Gxvect = np.arange(1, self.NX + 1)

        #meshgrid
        Gxmesh, Gymesh = np.meshgrid(Gxvect, Gyvect)

        #normalizers
        sigx = 1. / (2. * sigmax**2)
        sigy = 1. / (2. * sigmay**2)

        #gaussian
        exponent = (sigx * (Gxmesh - x)**2 + sigy * (Gymesh - y)**2)
        Gaussian = np.exp(-exponent) / (2. * np.pi * sigmax*sigmay)

        if self.shsettings['conservePeak']:
            #normalize to unity
            Gaussian /= np.max(Gaussian)

        output = dict(GaussianXmesh=Gxmesh, GaussianYmesh=Gymesh, Gaussian=Gaussian)

        return output


    def measureRefinedEllipticity(self):
        """
        Derive a refined iterated polarisability/ellipticity measurement for a given object.

        By default polarisability/ellipticity is defined in terms of the Gaussian weighted quadrupole moments.
        If self.shsettings['weighted'] is False then no weighting scheme is used.

        The number of iterations is defined in self.shsettings['iterations'].

        :return: centroids [indexing stars from 1], ellipticity (including projected e1 and e2), and R2
        :rtype: dict
        """
        self.shsettings['sampleSigma'] = self.shsettings['sigma'] / self.shsettings['pixelSize'] * \
                                       self.shsettings['platescale'] / self.shsettings['sampling']

        if self.log is not None: 
            self.log.info('Sample sigma used for weighting = %f' % self.shsettings['sampleSigma'])

        if self.shsettings['fixedPosition']:
            if self.log is not None: self.log.info('Using a fixed ')
            quad = dict(centreX=self.shsettings['fixedX'], centreY=self.shsettings['fixedY'])
        else:
            if self.log is not None: self.log.info('The initial estimate for the mean values are taken from the unweighted quadrupole moments.')
            quad = self.quadrupoles(self.data.copy())
            

        for x in range(self.shsettings['iterations']):
            if self.shsettings['weighted']:
                if self.log is not None: self.log.info('Iteration %i with circular symmetric Gaussian weights' % x)
                if self.shsettings['fixedPosition']:
                    gaussian = self.circular2DGaussian(self.shsettings['fixedX'],
                                                       self.shsettings['fixedY'],
                                                       self.shsettings['sampleSigma'])
                else:
                    gaussian = self.circular2DGaussian(quad['centreX']+1.,
                                                       quad['centreY']+1.,
                                                       self.shsettings['sampleSigma'])

                GaussianWeighted = self.data.copy() * gaussian['Gaussian'].copy()
            else:
                if self.log is not None: self.log.info('Iteration %i with no weighting' % x)
                GaussianWeighted = self.data.copy()

            quad = self.quadrupoles(GaussianWeighted.copy())
        
        # The squared radius R2 in um2
        R2 = quad['Qxx'] * self.shsettings['sampling']**2 + quad['Qyy'] * self.shsettings['sampling']**2.
        R2arcsec = R2 * (self.shsettings['pixelSize'] / self.shsettings['platescale'])**2.

        if self.shsettings['debug']:
            self.writeFITS(gaussian['Gaussian'], 'GaussianWeightingFunction.fits')
            self.writeFITS(GaussianWeighted, 'GaussianWeighted.fits')

        res = dict(centreX=quad['centreX']+1, centreY=quad['centreY']+1,
                   e1=quad['e1'], e2=quad['e2'],
                   ellipticity=quad['ellipticity'],
                   R2=R2,
                   R2arcsec=R2arcsec,
                   GaussianWeighted=GaussianWeighted,
                   a=quad['a'], b=quad['b'])
        return res


    def writeFITS(self, data, output):
        """
        Write out a FITS file using PyFITS.

        :param data: data to write to a FITS file
        :type data: ndarray
        :param output: name of the output file
        :type output: string

        :return: None
        """
        if os.path.isfile(output):
            os.remove(output)

        #create a new FITS file, using HDUList instance
        ofd = pf.HDUList(pf.PrimaryHDU())

        #new image HDU
        hdu = pf.ImageHDU(data=data)

        #update and verify the header
        hdu.header.add_history('If questions, please contact Ruyman Azzollini (r.azzollini at ucl.ac.uk).')
        hdu.header.add_history('This file has been created with the vison package at %s' \
                               % datetime.datetime.isoformat(datetime.datetime.now()))
        hdu.verify('fix')

        ofd.append(hdu)

        #write the actual file
        ofd.writeto(output)
        if self.log is not None: self.log.info('Wrote %s' % output)


class TestShape(unittest.TestCase):
    """
    Unit tests for the shape class.
    """
    def setUp(self):
        from vissim.support import logger as lg
        if self.log is not None: self.log = lg.setUpLogger('shapeTesting.log')

        self.psffile12x = '../data/psf12x.fits'
        self.psffile = '../data/psf1x.fits'
        self.tolerance = 1.e-7
        self.sigma = 40.0
        self.sigmax = 67.25
        self.sigmay = 24.15
        self.sigmax2 = 77.12343
        self.sigmay2 = 42.34543
        self.xcent = 500.
        self.ycent = 500.

        #create 2D Gaussians that will be used for testing
        self.GaussianCirc = Shapemeter(np.zeros((1000, 1000)), self.log).circular2DGaussian(self.xcent,
                                                                                                  self.ycent,
                                                                                                  self.sigma)['Gaussian']
        self.Gaussian = Shapemeter(np.zeros((1000, 1000)), self.log).ellip2DGaussian(self.xcent,
                                                                                      self.ycent,
                                                                                      self.sigmax,
                                                                                      self.sigmay)['Gaussian']
        self.Gaussian2 = Shapemeter(np.zeros((1000, 1000)), self.log).ellip2DGaussian(self.xcent,
                                                                                       self.ycent,
                                                                                       self.sigmax2,
                                                                                       self.sigmay2)['Gaussian']

    def test_ellipticity_noweighting_circular_Gaussian(self):
        expected = 0.0
        settings = dict(weighted=False)
        actual = Shapemeter(self.GaussianCirc, self.log, **settings).measureRefinedEllipticity()['ellipticity']
        self.assertAlmostEqual(expected, actual, msg='exp=%f, got=%f' % (expected, actual), delta=self.tolerance)


    def test_noweighting_Gaussian(self):
        expected = np.abs((self.sigmax**2 - self.sigmay**2) / (self.sigmax**2 + self.sigmay**2))
        settings = dict(weighted=False)
        actual = Shapemeter(self.Gaussian, self.log, **settings).measureRefinedEllipticity()
        ae = actual['ellipticity']
        ae1 = actual['e1']
        ae2 = actual['e2']
        R2 = actual['R2']
        R2exp = self.sigmax**2 + self.sigmay**2
        self.assertAlmostEqual(expected, ae, msg='exp=%f, got=%f' % (expected, ae), delta=self.tolerance)
        self.assertAlmostEqual(expected, ae1, msg='exp=%f, got=%f' % (expected, ae1), delta=self.tolerance)
        self.assertAlmostEqual(0.0, ae2, msg='exp=%f, got=%f' % (expected, ae2), delta=self.tolerance)
        self.assertAlmostEqual(R2exp, R2, msg='exp=%f, got=%f' % (R2exp, R2), delta=self.tolerance)


    def test_noweighting_Gaussian2(self):
        expected = np.abs((self.sigmax2**2 - self.sigmay2**2) / (self.sigmax2**2 + self.sigmay2**2))
        settings = dict(weighted=False, iterations=40)
        actual = Shapemeter(self.Gaussian2, self.log, **settings).measureRefinedEllipticity()
        ae = actual['ellipticity']
        ae1 = actual['e1']
        ae2 = actual['e2']
        R2 = actual['R2']
        R2exp = self.sigmax2**2 + self.sigmay2**2
        self.assertAlmostEqual(expected, ae, msg='exp=%f, got=%f' % (expected, ae), delta=self.tolerance)
        self.assertAlmostEqual(expected, ae1, msg='exp=%f, got=%f' % (expected, ae1), delta=self.tolerance)
        self.assertAlmostEqual(0.0, ae2, msg='exp=%f, got=%f' % (expected, ae2), delta=self.tolerance)
        self.assertAlmostEqual(R2exp, R2, msg='exp=%f, got=%f' % (R2exp, R2), delta=1e-4)


    def test_ellipticity_Gaussian(self):
        expected = np.abs((self.sigmax**2 - self.sigmay**2) / (self.sigmax**2 + self.sigmay**2))
        settings = dict(sigma=3000., iterations=40)
        actual = Shapemeter(self.Gaussian, self.log, **settings).measureRefinedEllipticity()['ellipticity']
        self.assertAlmostEqual(expected, actual, msg='exp=%f, got=%f' % (expected, actual), delta=1e-5)


    def test_centroiding_weighting_Gaussian(self):
        expected = self.xcent, self.ycent
        actual = Shapemeter(self.Gaussian, self.log).measureRefinedEllipticity()
        self.assertAlmostEqual(expected[0], actual['centreX'],
                               msg='exp=%f, got=%f' % (expected[0], actual['centreX']), delta=self.tolerance)
        self.assertAlmostEqual(expected[1], actual['centreY'],
                               msg='exp=%f, got=%f' % (expected[1], actual['centreY']), delta=self.tolerance)


    def test_R2_noweighting_circular_Gaussian(self):
        expected = 2 * self.sigma**2
        settings = dict(weighted=False)
        actual = Shapemeter(self.GaussianCirc, self.log, **settings).measureRefinedEllipticity()['R2']
        self.assertAlmostEqual(expected, actual, msg='exp=%f, got=%f' % (expected, actual), delta=self.tolerance)


    def test_PSF12x(self):
        expected = 0.045536
        R2exp = 5.010087
        data = pf.getdata(self.psffile12x)
        settings = dict(sampling=1/12.0)
        actual = Shapemeter(data, self.log, **settings).measureRefinedEllipticity()
        self.assertAlmostEqual(expected, actual['ellipticity'],
                               msg='exp=%f, got=%f' % (expected, actual['ellipticity']), delta=10*self.tolerance)
        self.assertAlmostEqual(R2exp, actual['R2'],
                               msg='exp=%f, got=%f' % (R2exp, actual['R2']), delta=10*self.tolerance)


    def test_PSF(self):
        expected = 0.045437
        R2exp = 4.959904
        data = pf.getdata(self.psffile)
        actual = Shapemeter(data, self.log).measureRefinedEllipticity()
        self.assertAlmostEqual(expected, actual['ellipticity'],
                               msg='exp=%f, got=%f' % (expected, actual['ellipticity']), delta=10*self.tolerance)
        self.assertAlmostEqual(R2exp, actual['R2'],
                               msg='exp=%f, got=%f' % (R2exp, actual['R2']), delta=10*self.tolerance)


    def test_gaussian_weighting(self):
        settings = dict(debug=True)
        data = pf.getdata(self.psffile)
        data /= np.max(data)
        _ = Shapemeter(data, self.log, **settings).measureRefinedEllipticity()
        gweighted = pf.getdata('GaussianWeighted.fits')
        self.assertAlmostEqual(np.max(data), np.max(gweighted))


if __name__ == '__main__':
    #testing section
    suite = unittest.TestLoader().loadTestsFromTestCase(TestShape)
    unittest.TextTestRunner(verbosity=3).run(suite)
