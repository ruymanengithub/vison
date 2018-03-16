#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""


Created on Fri Mar 16 14:38:51 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import copy

from vison.datamodel import ccd as ccdmod

# END IMPORT

def get_injprofile_tpnorm(ccdobj,vstart,vend):
    """Produces a 2D Map of charge injection to be used in trap-pumping analysis,
    to obtain dipole maps."""
    
    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    NAXIS1 = ccdobj.NAXIS1
    NAXIS2 = ccdobj.NAXIS2
    
    ccdobj_dummy = copy.deepcopy(ccdobj)
    
    for Q in Quads:
        
        qimg = ccdobj.extract_region(Q,'img',vstart=vstart,vend=vend,canonical=True).copy()
        qimgmod = qimg.mean(axis=1).reshape(qimg.shape[0],1).repeat(qimg.shape[1],axis=1)
        
        qmod = np.ones((NAXIS1/2,NAXIS2/2),dtype='float32')
        
        qmod[prescan:-overscan,vstart:vend] = qimgmod.copy()
        
        ccdobj_dummy.set_quad(qmod,Q,canonical=True,extension=-1)
    
    injprofile_tpnorm = ccdobj_dummy.extensions[-1].data.copy()
    
    return injprofile_tpnorm

def gen_raw_dpmap_vtpump(ccdobj,vstart=0,vend=ccdmod.NrowsCCD,extension=-1):
    """ """
    
    Quads = ccdobj.Quads
    prescan = ccdobj.prescan
    overscan = ccdobj.overscan
    
    for Q in Quads:
        qdata = ccdobj.get_quad(Q,canonical=True,extension=extension).copy()
        avgrow = qdata[prescan:-overscan,vstart:vend].mean(axis=1)
        dipmap = np.ones_like(qdata)
        dipmap[prescan:-overscan] = qdata[prescan:-overscan]/avgrow[:,None]
        ccdobj.set_quad(dipmap,Q,canonical=True,extension=extension)
        
    return ccdobj

    