# -*- coding: utf-8 -*-
"""
Doing Aperture Photometry of an object
======================================

Simple class to do Gaussian Fitting to a spot.

:requires: NumPy, astropy

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

Created on Thu Apr 20 16:42:47 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop

import numpy as np
from astropy.modeling import models, fitting
# END IMPORT

class Gaussmeter():
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
        
        self.gasettings = dict(x=NX/2,y=NY/2)

        self.gasettings.update(kwargs)
        

        for key, value in self.settings.iteritems():
            if self.log is not None: self.log.info('%s = %s' % (key, value))
    
    
    def fit_Gauss(self):
        """ """
        
        i00 = self.data.max()
        
        gaus = models.Gaussian2D(i00, self.NX, self.NY, x_stddev=0.5, y_stddev=0.5)
        gaus.theta.fixed = True  #fix angle
        p_init = gaus
        fit_p = fitting.LevMarLSQFitter()
        
        XX, YY = np.meshgrid(np.arange(self.NX), np.arange(self.NY),indexing='xy')
        p = fit_p(p_init, XX, YY, self.data)
        print p
        
        stop()