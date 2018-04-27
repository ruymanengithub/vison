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

offsets = dict(CCD1=2000, CCD2=2000, CCD3=2000)  # ADU
offsets_margins = [-200, 200]  # ADU
offsets_lims = dict()
for iCCD in [1, 2, 3]:
    CCDkey = 'CCD%i' % iCCD
    offsets_lims[CCDkey] = [offsets[CCDkey] +
                            offsets_margins[0], offsets[CCDkey]+offsets_margins[1]]

# offsets_gradients = dict(CCDX=[[prescan-prescan+-margin],[image-prescan+-margin],[overscan-prescan+-margin]])
offsets_gradients = dict(CCD1=dict(pre=[0, 0], img=[5, 10], ove=[5, 10]), 
                         CCD2=dict(pre=[0, 0], img=[5, 10], ove=[5, 10]),
                         CCD3=dict(pre=[0, 0], img=[5, 10], ove=[5, 10]))  # ADU
RONs = dict(CCD1=1.4, CCD2=1.4, CCD3=1.4)  # ADU, STD
RONs_margins = [-0.2, 0.2]  # ADU, STD
RONs_lims = dict()
for CCD in [1, 2, 3]:
    CCDkey = 'CCD%i' % CCD
    RONs_lims[CCDkey] = [RONs[CCDkey] +
                         RONs_margins[0], RONs[CCDkey]+RONs_margins[1]]


perf_rdout = dict(offsets=offsets, offsets_margins=offsets_margins,
                  offsets_lims=offsets_lims,
                  offsets_gradients=offsets_gradients,
                  RONs=RONs, RONs_margins=RONs_margins,
                  RONs_lims=RONs_lims)

gains = dict(CCD1=3.3, CCD2=3.3, CCD3=3.3)
