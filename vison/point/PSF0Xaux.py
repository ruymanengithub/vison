# -*- coding: utf-8 -*-
"""
Script to analyze Test PSF0X

PSF vs. Fluence, and Wavelength
   PSF01 - nominal temperature
   PSF02 - alternative temperatures

Tasks:

    - select exposures, get file names, extract data.
    - subtract offset level.
    - divide by Flat-field.
    - crop stamps of the sources on each CCD/Quadrant.
       - save snapshot figures of sources.
    - for each source:
       - measure shape using weighted moments
       - measure shape using Gaussian Fit
       - Forward Model the optomechanic+detector PSF

Created on Fri Nov 25 19:14:12 2016

@author: raf
"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from collections import OrderedDict

from vison.plot import baseclasses as plbaseclasses

# END IMPORT

class genplot_stats_beam:
    """ """



PSF0Xfigs = dict()
PSF0Xfigs['P0Xchecks_offsets'] = genplot_stats_beam
PSF0Xfigs['P0Xchecks_stds'] = genplot_stats_beam
PSF0Xfigs['BlueScreen'] = plbaseclasses.BlueScreen



            
            