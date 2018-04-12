#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

EUCLID-VIS Ground Calibration Campaign

Model of the calibration OGSE




Created on Fri Sep  8 12:11:55 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# OVERHEADS

overheads = dict(ff_readout=72.,
                 writing=25.,
                 flat2point=40.,
                 delta_focus=0.,
                 fw=2.)

# saturation exposure times for flat-source, s
# updated on 3rd FEB 2018 after Cold Commissioning of Fac#1

tFWC_flat = dict(nm590=0.82,
                 nm640=0.82,
                 nm730=0.94,
                 nm800=1.87,
                 nm880=9.4,
                 nm0=1.6)  # after DRY-RUN...

# saturation exposure times for point-source, ms, orientative

# updated on 3rd FEB 2018 after Cold Commissioning of Fac#1

tFWC_point = dict(nm590=0.17,
                  nm640=0.22,
                  nm730=0.21,
                  nm800=0.39,
                  nm880=1.95,
                  nm0=3.)  # ms

# FWHM (lambda), in pixels, GUESS

fwhm_lambda = dict(nm590=1.3,
                   nm640=1.4,
                   nm730=1.5,
                   nm800=1.5,
                   nm880=1.6,
                   nm0=1.5)  # ms

# Filter Wheel Positions

FW = dict(F1=590, F2=640,
          F3=730, F4=800,
          F5=880, F6=0)


def get_FW_ID(wavelength):
    """returns FW key corresponding to input wavelength.
    :param wavelength: integer, wavelength.
    """
    return [key for key in FW if FW[key] == wavelength][0]
