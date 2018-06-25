#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

DS9 Regions tool.

Created on Fri May 18 15:02:07 2018

:author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np

# END IMPORT


def _replace_default(A, Nregs, default):
    if A is not None:
        assert len(A) == Nregs
    else:
        A = np.ones(Nregs, dtype='float32') * default
    return A


def get_body_circles(X, Y, R=None, radius=6.0):
    """ """

    Nregs = len(X)
    assert len(Y) == Nregs
    R = _replace_default(R, Nregs, default=radius)

    body = []
    for i in range(Nregs):
        body.append('circle(%.1f,%.1f,%.1f)' %
                    (X[i]+1., Y[i]+1., R[i]))

    return body


def get_body_ellipses(X, Y, A=None, B=None, THETA=None):
    """ """
    radius = 6.0

    Nregs = len(X)
    assert len(Y) == Nregs

    A = _replace_default(A, Nregs, default=radius)
    B = _replace_default(B, Nregs, default=radius)
    THETA = _replace_default(THETA, Nregs, default=0.)

    body = []

    for i in range(Nregs):
        body.append('ellipse(%.1f,%.1f,%.1f,%.1f,%.1f)' %
                    (X[i]+1., Y[i]+1., A[i], B[i], THETA[i]))

    return body


def save_spots_as_ds9regs(data, regfilename=None, regfile=None, regtype='circle', clobber=True):
    """ """
    color = 'green'
    width = 2

    assert (regfilename is not None) ^ (regfile is not None)

    fbody_dict = dict(circle=get_body_circles,
                      ellipse=get_body_ellipses)

    hdr_main = 'global color=%s dashlist=8 3 width=%i ' +\
        'font="helvetica 10 normal roman" select=1 highlite=1' +\
        ' dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'

    hdr = ['# Region file format: DS9 version 4.1',
           hdr_main % (color, width),
           'physical']

    body = fbody_dict[regtype](**data)

    def myprinter(f, hdr, body):
        for line in hdr:
            print >> f, line
        for line in body:
            print >>f, line

    if regfilename is not None:

        if clobber and os.path.exists(regfilename):
            os.system('rm %s' % regfilename)

        with open(regfilename, 'wa') as f:
            myprinter(f, hdr, body)
            f.close()
    elif regfile is not None:
        myprinter(regfile, hdr, body)
        regfile.close()
