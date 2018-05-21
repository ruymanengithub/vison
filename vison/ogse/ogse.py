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


# FOCUS 

# TODO: NEEDS to be addressed by CHAMBER-ID
mirror_nom = dict(F1=86., F2=86., F3=86., F4=86., F5=86.,
                  F6=86.)  # nominal focus positions of mirror


# SATURATION TIMES

# saturation exposure times for flat-source, s
# updated on 3rd FEB 2018 after Cold Commissioning of Fac#1

#tFWC_flat = dict(nm590=0.82,
#                 nm640=0.82,
#                 nm730=0.94,
#                 nm800=1.87,
#                 nm880=9.4,
#                 nm0=1.6)  # after DRY-RUN...

# saturation exposure times for flat-source, s
# updated on 17 MAY 2018 before DRY-RUN-JUN18

#tFWC_flat = dict(nm590=1.02,
#                 nm640=0.78,
#                 nm730=0.94,
#                 nm800=0.53,
#                 nm880=2.02,
#                 nm0=1.6)  # GUESSING from tFWC_point...       

tFWC_flat = dict(nm590=100.,
                 nm640=100.,
                 nm730=100.,
                 nm800=200.,
                 nm880=200.,
                 nm0=200.)  # Using diffuser... 21/MAY/2018   
                 
# saturation exposure times for point-source, s, orientative

# updated on 3rd FEB 2018 after Cold Commissioning of Fac#1

#tFWC_point = dict(nm590=0.17,
#                  nm640=0.22,
#                  nm730=0.21,
#                  nm800=0.39,
#                  nm880=1.95,
#                  nm0=3.)  # s

# updated on 17th May 2018 before DRY-RUN-JUN18
tFWC_point = dict(nm590=0.21,
                  nm640=0.21,
                  nm730=0.21,
                  nm800=0.11,
                  nm880=0.82,
                  nm0=3.)  # s

# FWHM (lambda), in pixels, GUESS

fwhm_lambda = dict(nm590=1.3,
                   nm640=1.4,
                   nm730=1.5,
                   nm800=1.5,
                   nm880=1.6,
                   nm0=1.5)  # s


# Filter Wheel Positions

FW = dict(F1=590, F2=640,
          F3=730, F4=800,
          F5=880, F6=0)


def get_FW_ID(wavelength):
    """returns FW key corresponding to input wavelength.
    :param wavelength: integer, wavelength.
    """
    return [key for key in FW if FW[key] == wavelength][0]
