#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Common Values which are used by functions and classes throughout pipeline.


Created on Tue Jan 16 10:53:40 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from vison.datamodel.ccd import CCD
# END IMPORT

elvis = '6.5.X'
Quads = ['E', 'F', 'G', 'H']

ccdobj = CCD()
prescan = ccdobj.prescan
overscan = ccdobj.overscan
imgheight = ccdobj.NAXIS2/2
quad_width = ccdobj.NAXIS1/2
imgwidth = quad_width - prescan - overscan

#sumwell = dict(fwd_bas=[9.3, 4.725], fwd_e2v=[7.95, 5.825],
#               rwd_bas_v=[6.5, 0.175], rwd_bas_s=[6.5, 0.175],
#               rwd_bas_vs=[6.5, 0.175])
# Update to sumwell in 17 OCT 2018
sumwell = dict(fwd_bas=[8.4, 5.625], 
            fwd_e2v=[7.15, 6.725],
            rwd_bas_v=[8.4, 5.625],
            rwd_bas_s=[6.5, 0.175],
            rwd_bas_vs=[6.5, 0.175])

