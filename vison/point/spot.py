# -*- coding: utf-8 -*-
"""

Spot Stamp Class.

Created on Thu Apr 20 15:35:08 2017


:author: Ruyman Azzollini


"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np

from shape import Shapemeter
from photom import Photometer
from gauss import Gaussmeter
from vison.datamodel.ccd import gain
# END IMPORT


class Spot(Shapemeter, Photometer, Gaussmeter):
    """
    Provides methods to do point-source analysis on a stamp.
    Aimed at basic analysis:
     - Photometry
     - Quadrupole Moments
     - Gaussian Fit

    :param data: stamp to be analysed.
    :type data: np.ndarray
    :param log: logger
    :type log: instance
    :param kwargs: additional keyword arguments
    :type kwargs: dict

    Settings dictionary contains all parameter values needed.
    """

    def __init__(self, data, log=None, verbose=False, lowerleft=(None,),
                 **kwargs):
        """
        :param data: stamp to be analysed.
        :type data: ndarray
        :param log: logger
        :type log: instance
        :param verbose: verbosity switch
        :type verbose: bool
        :param lowerleft: lower left corner coordinates (x,y)
        :type lowerleft: tuple
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        Settings dictionary contains all parameter values needed.
        """

        #super(Spot, self).__init__(data,log=log,**kwargs)

        Shapemeter.__init__(self, data, log, verbose, **kwargs)
        Photometer.__init__(self, data, log, verbose, **kwargs)
        Gaussmeter.__init__(self, data, log, verbose, **kwargs)

        self.data = data.copy()
        self.log = log

        NX, NY = self.data.shape
        self.NX = NX
        self.NY = NY

        self.xcen = NX / 2
        self.ycen = NY / 2

        if lowerleft[0] is not None:
            self.x0 = lowerleft[0]
            self.y0 = lowerleft[1]

    def get_photom(self):
        """

        measurements:'apflu','eapflu','bgd','ebgd'
        """
        raise NotImplementedError

        res = dict()
        return res

    def get_shape_Gauss(self):
        """

        :return: res = dict(i0,ei0,x,ex,y,ey,
                    sigma_x,esigma_x,sigmay,esigma_y,
                    fwhm_x,efwhm_x,
                    fwhm_y,efwhm_y,
                    fluence,efluence)

        """

        Gpars, eGpars = self.fit_Gauss()

        res = dict(i0=Gpars[0], ei0=eGpars[0],
                   x=Gpars[1], ex=eGpars[1], y=Gpars[2], ey=eGpars[2],
                   sigma_x=Gpars[3], esigma_x=eGpars[3],
                   sigma_y=Gpars[4], esigma_y=eGpars[4])
        res['fwhm_x'] = res['sigma_x'] * 2.355
        res['efwhm_x'] = res['esigma_x'] * 2.355
        res['fwhm_y'] = res['sigma_y'] * 2.355
        res['efwhm_y'] = res['esigma_y'] * 2.355
        sigma_maj = np.max([res['sigma_x'], res['sigma_y']])
        sigma_min = np.min([res['sigma_x'], res['sigma_y']])
        q = sigma_min / sigma_maj

        res['fluence'] = 2. * np.pi * res['i0'] * q * sigma_maj**2.
        res['efluence'] = 2. * np.pi * res['ei0'] * q * sigma_maj**2.

        return res

    def get_shape_Moments(self):
        """

        :return: res = dict(x,y,ellip,e1,e2,a,b)

        """

        rawres = self.measureRefinedEllipticity()
        res = dict(x=rawres['centreX'] - 1., y=rawres['centreY'] - 1.,
                   ellip=rawres['ellipticity'], e1=rawres['e1'],
                   e2=rawres['e2'], R2=rawres['R2'],
                   a=rawres['a'], b=rawres['b'])

        return res

    def get_shape_easy(self, method='G', debug=False):
        """

        """

        if debug:
            return dict()

        if method == 'G':
            res = self.get_shape_Gauss()

        if method == 'M':
            res = self.get_shape_Moments()

        return res

    def measure_basic(self, rap=10, rin=15, rout=-1, gain=gain, debug=False):
        """
        # TODO:
        #   get basic statistics, measure and subtract background
        #   update centroid
        #   do aperture photometry
        #   pack-up results and return

        :parameter rap: source extraction aperture radius.
        :parameter rin: inner radius of background annulus.
        :parameter rout: outer radius of background annulus (-1 to set bound by image area).
        :parameter gain: image gain (e-/ADU).

        """

        if rout < 0:
            rout = max(self.data.shape)
        if debug:
            return dict(bgd=0., peak=1., fluence=1., efluence=0.,
                        x=1., y=1., fwhmx=1., fwhmy=1.)

        bgd, sbgd = self.measure_bgd(rin, rout)
        self.sub_bgd(rin, rout)
        peak = self.data.max()
        xcen, ycen, fwhmx, fwhmy = self.get_centroid(full=True)
        centre = self.xcen, self.ycen

        flu, eflu = self.doap_photom(centre, rap, rin, rout, gain=gain,
                                     doErrors=True,
                                     subbgd=True)

        x, y = self.xcen, self.ycen
        x_ccd = x + self.x0
        y_ccd = y + self.y0

        res = dict(bgd=bgd, peak=peak, fluence=flu, efluence=eflu,
                   x=x, y=y, x_ccd=x_ccd, y_ccd=y_ccd,
                   fwhmx=fwhmx, fwhmy=fwhmy)

        return res
