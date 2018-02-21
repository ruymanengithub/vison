#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary functions to BIAS01 and DARK01 test scripts.

Created on Tue Nov 14 13:25:01 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop


# END IMPORT



def get_darkdefects_mask(ccdobj,threshold,subbgd=True,extension=-1):
    """ """
    
    Quads = ccdobj.Quads
    
    if subbgd:
        
        for Quad in Quads:
            
            Qimgmodel = ccdobj.get_region2Dmodel(Quad,area='img',kind='poly2D',
                pdegree=5,doFilter=False,filtsize=1,
                vstart=0,vend=2066,canonical=False,extension=extension)
            
            
        ccdobj.sub_bias()
    