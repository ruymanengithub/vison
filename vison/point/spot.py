# -*- coding: utf-8 -*-
"""



Created on Thu Apr 20 15:35:08 2017

@author: raf
"""

# IMPORT STUFF
from shape import Shapemeter
from photom import Photometer
from gauss import Gaussmeter
# END IMPORT

class Spot(Shapemeter,Photometer,Gaussmeter):
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

        sizeY, sizeX = self.data.shape
        self.settings = dict(sizeX=sizeX,
                             sizeY=sizeY)

        self.settings.update(kwargs)
        

        for key, value in self.settings.iteritems():
            if self.log is not None: self.log.info('%s = %s' % (key, value))
        
    