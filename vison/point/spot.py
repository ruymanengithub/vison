# -*- coding: utf-8 -*-
"""




:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

Created on Thu Apr 20 15:35:08 2017

"""

# IMPORT STUFF
from pdb import set_trace as stop
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
        
        #super(Spot, self).__init__(data,log=log,**kwargs)
        
        Shapemeter.__init__(self,data,log,**kwargs)
        Photometer.__init__(self,data,log,**kwargs)
        Gaussmeter.__init__(self,data,log,**kwargs)
        
        
        self.data = data.copy()
        self.log = log

        NX, NY = self.data.shape
        self.NX = NX
        self.NY = NY
        
        self.xcen = NX/2
        self.ycen = NY/2
