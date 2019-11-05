# -*- coding: utf-8 -*-
"""


:author: Ruyman Azzollini

Created on Thu Apr 20 18:56:40 2017

"""


class SpotBase(object):
    """ """

    def __init__(self, data, log=None, verbose=False):
        """
        :param data: stamp to be analysed.
        :type data: ndarray
        :param log: logger
        :type log: instance
        :param verbose: bool
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        Settings dictionary contains all parameter values needed.
        """

        self.data = data.copy()
        self.log = log
        self.verbose = verbose

        NX, NY = self.data.shape
        self.NX = NX
        self.NY = NY

        self.xcen = NX / 2
        self.ycen = NY / 2
