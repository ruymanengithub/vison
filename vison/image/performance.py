#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Performance parameters of the ROE+CCDs.
Compilation of CCD offsets, offset gradients, RONs... used for 
checks.

Created on Wed Nov  1 09:57:44 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

offsets = dict(CCD1=2000,CCD2=2000,CCD3=2000) # ADU
offsets_margins = [-200,200] # ADU
offsets_gradients = dict(CCD1=[[0,0],[5,10],[5,10]],CCD2=[[0,0],[0,10],[0,10]],
                         CCD3=[[0,0],[0,10],[0,10]]) # ADU
RONs = dict(CCD1=1.4,CCD2=1.4,CCD3=1.4) # ADU, STD
RONs_margins = [-0.2,0.2] # ADU, STD

perf_rdout = dict(offsets=offsets,offsets_margins=offsets_margins,
                  offsets_gradients=offsets_gradients,
                  RONs=RONs,RONs_margins=RONs_margins)

