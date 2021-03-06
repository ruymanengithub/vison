#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Model of the calibration OGSE

Created on Fri Sep  8 12:11:55 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
from collections import OrderedDict
import copy

from vison.support import vjson
from vison import ogse_profiles
from vison.point import startracker as strmod
# END IMPORT

_path = os.path.abspath(ogse_profiles.__file__)
ogsepath = os.path.dirname(_path)

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

tFWC_flat = dict(nm590=0.82,
                 nm640=0.82,
                 nm730=0.94,
                 nm800=1.87,
                 nm880=9.4,
                 nm0=1.6)  # after DRY-RUN-FEB18

# saturation exposure times for flat-source, s
# updated on 17 MAY 2018 before DRY-RUN-JUN18

# tFWC_flat = dict(nm590=100.,
#                 nm640=100.,
#                 nm730=100.,
#                 nm800=200.,
#                 nm880=200.,
#                 nm0=200.)  # Using diffuser... 21/MAY/2018

# saturation exposure times for point-source, s, orientative

# updated on 3rd FEB 2018 after Cold Commissioning of Fac#1

tFWC_point = dict(nm590=0.17,
                  nm640=0.22,
                  nm730=0.21,
                  nm800=0.39,
                  nm880=1.95,
                  nm0=3.)  # s

# updated on 17th May 2018 before DRY-RUN-JUN18
# tFWC_point = dict(nm590=0.21,
#                  nm640=0.21,
#                  nm730=0.21,
#                  nm800=0.11,
#                  nm880=0.82,
#                  nm0=3.)  # s

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


default_profile = dict(
    tFWC_flat=dict(nm590=0.82,
                   nm640=0.82,
                   nm730=0.94,
                   nm800=1.87,
                   nm880=9.4,
                   nm0=1.6),  # after DRY-RUN-FEB18
    tFWC_point=dict(nm590=0.17,
                    nm640=0.22,
                    nm730=0.21,
                    nm800=0.39,
                    nm880=1.95,
                    nm0=3.),  # s, after DRY-RUN-FEB18

    mirror_nom=dict(F1=86., F2=86., F3=86., F4=86.,
                    F5=86., F6=86.),  # nominal focus positions of mirror

    FW=dict(F1=590, F2=640,
            F3=730, F4=800,
            F5=880, F6=0),

    fwhm_lambda=dict(nm590=1.3,
                     nm640=1.4,
                     nm730=1.5,
                     nm800=1.5,
                     nm880=1.6,
                     nm0=1.5)  # s

)


def get_FW_ID(wavelength, FW=FW):  # REDTAG: DEFAULT FW ON TESTS ONLY
    """returns FW key corresponding to input wavelength.
    :param wavelength: integer, wavelength.
    """
    return [key for key in FW if FW[key] == wavelength][0]


def get_wavelength(FW_ID, FW=FW):
    return FW['F%i' % FW_ID]


class Ogse(object):

    def __init__(self, CHAMBER=None, withpover=True, labels=None):

        self.CCDs = [1, 2, 3]
        if labels is None:
            self.labels = []
        else:
            assert isinstance(labels, list)
            self.labels = labels
        self.CHAMBER = CHAMBER
        self.profile = OrderedDict()
        self.load_ogse_profile(CHAMBER)
        self.startrackers = OrderedDict()
        self.load_startrackers(withpover=withpover)

    def get_profile_path(self, CHAMBER):
        profilef = 'CHAMBER_%s_profile.json' % CHAMBER
        fpath = os.path.join(ogsepath, profilef)
        return fpath

    def load_ogse_profile(self, CHAMBER=None):
        """ """

        if CHAMBER is None:
            self.profile = default_profile.copy()
            return

        fpath = self.get_profile_path(CHAMBER)
        profile = vjson.load_jsonfile(fpath)

        self.profile = profile.copy()

    def get_pspattern_paths(self):
        pspattern_paths = OrderedDict()
        for CCD in self.CCDs:
            pspatternf = str(self.profile['PATTERN']['CCD%i' % CCD])
            fpath = os.path.join(ogsepath, pspatternf)
            pspattern_paths['CCD%i' % CCD] = fpath
        return pspattern_paths

    def load_startrackers(self, withpover=True):
        """ """

        pspattern_paths = self.get_pspattern_paths()

        for CCD in self.CCDs:

            if len(self.labels) > 0:

                for label in self.labels:

                    istartracker = strmod.StarTracker(
                        CCD, withpover, patt_files=pspattern_paths)
                    self.startrackers['CCD%i' % CCD][label] = copy.deepcopy(istartracker)

            else:

                istartracker = strmod.StarTracker(
                    CCD, withpover, patt_files=pspattern_paths)
                self.startrackers['CCD%i' % CCD] = copy.deepcopy(istartracker)

    def save_ogse_profile(self, jsonpath):
        """ """
        profile = self.profile
        vjson.save_jsonfile(profile, jsonpath)

    def get_FW_ID(self, wavelength):
        """ """
        return get_FW_ID(wavelength, self.profile['FW'])

    def get_wavelength(self, FW_ID):
        return get_wavelength(FW_ID, self.profile['FW'])
