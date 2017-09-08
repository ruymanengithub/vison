#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

EUCLID-VIS Ground Calibration Campaign

Model of the calibration OGSE




Created on Fri Sep  8 12:11:55 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# saturation exposure times for flat-source, ms, orientative

tFWC_flat = dict(nm590=2.E3,
                nm640=2.E3,
                nm730=2.E3,
                nm800=2.E3,
                nm880=2.E3,
                ND=200.*1.E3) 

# saturation exposure times for point-source, ms, orientative


tFWC_point = dict(nm590=2.E3,
                nm640=2.E3,
                nm730=2.E3,
                nm800=2.E3,
                nm880=2.E3,
                ND=200.*1.E3) # ms

# FWHM (lambda), in pixels, orientative
                  
fwhm_lambda = dict(nm590=1.3,
                nm640=1.4,
                nm730=1.5,
                nm800=1.5,
                nm880=1.6,
                ND=1.5) # ms

FW = dict(Filter1=590,Filter2=640,
          Filter3=730,Filter4=800,
          Filter5=880,Filter6=0)



def get_FW_ID(wavelength):
    """returns FW key corresponding to input wavelength.
    :param wavelength: integer, wavelength.
    """
    return [key for key in FW if FW[key] == wavelength][0] 


                   
