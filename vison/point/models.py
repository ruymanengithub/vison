# -*- coding: utf-8 -*-
"""

FM-Calib. Campaign.

Library module with models for processing of point-source imaging data.

Created on Wed Apr 19 11:47:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
# END IMPORT


def fgauss2D(x,y,p):
    """
    A gaussian fitting function where
     p[0] = amplitude
     p[1] = x0
     p[2] = y0
     p[3] = sigmax
     p[4] = sigmay
     p[5] = floor
    """
    return p[0] * np.exp(-( 
                          (x-p[1])**2./(2.0*p[3]**2.) +
                          (y-p[2])**2./(2.*p[4]**2.)
                        )) + p[5]
    


