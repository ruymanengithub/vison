#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Unit-testing for CCDPile class.


Created on Mon May  7 09:47:07 2018

:author: Ruyman Azzollini

"""
# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import unittest
from vison.datamodel import ccd as ccdmod
# END IMPORT


class TestCCDPileClass(unittest.TestCase):

    def setUp(self):
        self.NAXIS1 = 4238
        #self.wQ = self.NAXIS1/2
        self.NAXIS2 = 4172
        #self.hQ = self.NAXIS2/2
        self.offset = 1000.
        # bad pixel on each Quad, canonical orientation
        #self.badpixel = (100, 100)
        self.prescan = ccdmod.prescan
        self.overscan = ccdmod.overscan
        self.voverscan = ccdmod.voverscan
        self.tolerance = 1.e-7
        self.gain = 3.5
        self.RON = 4.5 / self.gain
        self.StackSize = 5

        baseimg = np.ones((self.NAXIS1, self.NAXIS2), dtype='float32')
        eimg = np.ones((self.NAXIS1, self.NAXIS2), dtype='float32') * self.RON

        self.ccdobjs = []

        for i in range(1, self.StackSize + 1):

            self.ccdobjs.append(ccdmod.CCD(withpover=True))
            img = baseimg * i
            img[:i, :] = np.nan

            self.ccdobjs[i - 1].add_extension(data=img, label='IMAGE')
            self.ccdobjs[i - 1].add_extension(data=eimg, label='UNCERTAINTY')

            mask = np.isnan(img)
            self.ccdobjs[i - 1].get_mask(mask)

        self.ccdpile = ccdmod.CCDPile(ccdobjList=self.ccdobjs, extension=0,
                                      withpover=True)

    #@unittest.skip("REDTAG")
    #@unittest.expectedFailure
    def test_stack_mean(self):

        stacked, stdimg = self.ccdpile.stack(method='mean', dostd=True)

        residuals_stack = stacked[0:6, 0].data -\
            np.array([0., 1., 1.5, 2., 2.5, 3.])
        residuals_std = stdimg[0:6, 0].data - \
            np.array([0., 0., 0.5, 0.81649661, 1.11803401,
                      1.41421354])
        self.assertAlmostEqual(residuals_stack.sum(), 0.)
        self.assertAlmostEqual(residuals_std.sum(), 0.)

    #@unittest.expectedFailure
    def test_stack_median(self):

        stacked, stdimg = self.ccdpile.stack(method='median', dostd=True)

        residuals_stack = stacked[0:6, 0].data -\
            np.array([0., 1., 1.5, 2., 2.5, 3.])
        residuals_std = stdimg[0:6, 0].data - \
            np.array([0., 0., 0.5, 0.81649661, 1.11803401,
                      1.41421354])
        self.assertAlmostEqual(residuals_stack.sum(), 0.)
        self.assertAlmostEqual(residuals_std.sum(), 0.)


if __name__ == '__main__':
    unittest.main()
