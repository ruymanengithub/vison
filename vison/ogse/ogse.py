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

overheads = dict(ff_readout = 72.,
                 writing = 25.,
                 flat2point = 40.,
                 fw=2.)

# saturation exposure times for flat-source, s, place-holders

tFWC_flat = dict(nm590=2.,
                nm640=2.,
                nm730=2.,
                nm800=2.,
                nm890=2.,
                nm0=200.) 

# saturation exposure times for point-source, ms, orientative


tFWC_point = dict(nm590=2.,
                nm640=2.,
                nm730=2.,
                nm800=2.,
                nm890=2.,
                nm0=200.) # ms

# FWHM (lambda), in pixels, orientative
                  
fwhm_lambda = dict(nm590=1.3,
                nm640=1.4,
                nm730=1.5,
                nm800=1.5,
                nm890=1.6,
                nm0=1.5) # ms

# Filter Wheel Positions

FW = dict(F1=590,F2=640,
          F3=730,F4=800,
          F5=890,F6=0)



def get_FW_ID(wavelength):
    """returns FW key corresponding to input wavelength.
    :param wavelength: integer, wavelength.
    """
    return [key for key in FW if FW[key] == wavelength][0] 


                   
