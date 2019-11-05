# -*- coding: utf-8 -*-
"""
Gaussian Model of Point-like Sources
====================================

Simple class to do Gaussian Fitting to a spot.

:requires: NumPy, astropy

Created on Thu Apr 20 16:42:47 2017


:author: Ruyman Azzollini



"""

# IMPORT STUFF
from pdb import set_trace as stop

import numpy as np
from astropy.modeling import models, fitting

from vison.point.basis import SpotBase
# END IMPORT


class Gaussmeter(SpotBase):
    """
    Provides methods to measure the shape of an object using
    a 2D Gaussian Model.

    :param data: stamp to be analysed.
    :type data: np.ndarray
    :param log: logger
    :type log: instance
    :param kwargs: additional keyword arguments
    :type kwargs: dict

    Settings dictionary contains all parameter values needed.
    """

    def __init__(self, data, log=None, verbose=False, **kwargs):
        """
        :param data: stamp to be analysed.
        :type data: ndarray
        :param log: logger
        :type log: instance
        :param verbose: bool switch
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        Settings dictionary contains all parameter values needed.
        """

        super(Gaussmeter, self).__init__(data, log)

        self.gasettings = dict()

        self.gasettings.update(kwargs)

        for key, value in self.gasettings.iteritems():
            if self.log is not None:
                self.log.info('%s = %s' % (key, value))

    def fit_Gauss(self):
        """ """

        i00 = self.data.max()

        gaus = models.Gaussian2D(
            i00, self.xcen, self.ycen, x_stddev=0.5, y_stddev=0.5)
        gaus.theta.fixed = True  # fix angle
        p_init = gaus
        fit_p = fitting.LevMarLSQFitter()

        XX, YY = np.meshgrid(np.arange(self.NX),
                             np.arange(self.NY), indexing='xy')
        rawres = fit_p(p_init, XX, YY, self.data)

        params = rawres.parameters[0:-1]  # theta is fixed
        try:
            eparams = np.diag(fit_p.fit_info['param_cov'])**0.5
        except BaseException:
            Emsg = 'Gaussmeter.fit_Gauss did not converge'
            if self.log is not None:
                self.log.info(Emsg)
            else:
                print Emsg
            eparams = np.zeros_like(params) + np.nan

        return params, eparams
