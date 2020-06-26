#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Methods to simulate data. Used by ccd.CCD class.

:History:
Created on Wed Apr  4 11:13:30 2018

:author: Ruyman Azzollini


"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np

# END IMPORT


def simadd_flatilum(ccdobj, levels=None, extension=-1):
    """ """
    if levels is None:
        levels = dict(E=0., F=0., G=0., H=0.)

    for Q in ccdobj.Quads:
        quaddata = ccdobj.get_quad(Q, canonical=True, extension=extension)
        quaddata[ccdobj.prescan:-ccdobj.overscan,
                 0:ccdobj.NrowsCCD] += levels[Q]
        ccdobj.set_quad(quaddata, Q, canonical=True, extension=extension)


def simadd_points(ccdobj, flux, fwhm, CCDID='CCD1', dx=0, dy=0, extension=-1):
    """ """
    from vison.point.lib import Point_CooNom
    from vison.point.models import fgauss2D

    sigma = fwhm / 2.355
    i0 = flux / (2. * np.pi * sigma**2.)

    nx = 15
    ny = 15

    for Q in ccdobj.Quads:
        quaddata = ccdobj.get_quad(
            Q, canonical=False, extension=extension).copy()

        point_keys = Point_CooNom[CCDID][Q].keys()

        B = ccdobj.QuadBound[Q]

        for pkey in point_keys:
            xp, yp = Point_CooNom[CCDID][Q][pkey]

            x0 = xp - B[0] + dx
            y0 = yp - B[2] + dy

            xmin = int(np.round(x0)) - nx / 2
            xmax = xmin + nx
            ymin = int(np.round(y0)) - ny / 2
            ymax = ymin + ny

            x = np.arange(xmin, xmax)
            y = np.arange(ymin, ymax)

            xx, yy = np.meshgrid(x, y, indexing='xy')

            p = [i0, x0, y0, sigma, sigma, 0.]
            point_stamp = fgauss2D(xx, yy, p)

            quaddata[xmin:xmax, ymin:ymax] += point_stamp.copy()

        ccdobj.set_quad(quaddata, Q, canonical=False, extension=extension)


def simadd_bias(ccdobj, levels=None, extension=-1):

    if levels is None:
        levels = dict(E=2000, F=2000, G=2000, H=2000)

    for Q in ccdobj.Quads:
        B = ccdobj.QuadBound[Q]
        ccdobj.extensions[extension].data[B[0]:B[1], B[2]:B[3]] += levels[Q]


def simadd_ron(ccdobj, extension=-1):
    """ """

    for Q in ccdobj.Quads:

        B = ccdobj.QuadBound[Q]

        Qnoise = np.random.normal(
            loc=0.0, scale=ccdobj.rn[Q], size=(ccdobj.wQ, ccdobj.hQ))
        ccdobj.extensions[extension].data[B[0]:B[1], B[2]:B[3]] += Qnoise


def simadd_poisson(ccdobj, extension=-1):
    """ """

    gain = ccdobj.gain['E']

    image_e = ccdobj.extensions[extension].data.copy() * gain

    rounded = np.rint(image_e)
    # ugly workaround for multiple rounding operations...
    residual = image_e - rounded
    rounded[rounded < 0.0] = 0.0
    image_e = np.random.poisson(rounded).astype(np.float64)
    image_e += residual

    image = image_e / gain

    ccdobj.extensions[extension].data = image.copy()


def sim_window(ccdobj, vstart, vend, extension=-1):
    """ """

    mask = np.ones((ccdobj.NAXIS1 / 2, ccdobj.NAXIS2 / 2), dtype='int8')
    mask[:, vstart:vend] = 0

    for Q in ccdobj.Quads:
        qdata = ccdobj.get_quad(Q, canonical=True, extension=extension)
        qdata[np.where(mask)] = 0
        ccdobj.set_quad(qdata, Q, canonical=True, extension=extension)


def simadd_injection(ccdobj, levels, on=1E6, off=0, extension=-1):

    on = min(on, ccdobj.NrowsCCD)

    ccdobj.simadd_flatilum(levels=levels, extension=-1)

    # ON/OFF masking

    NX, NY = ccdobj.NAXIS1 / 2, ccdobj.NAXIS2 / 2

    mask_onoff = np.zeros((NX, NY), dtype='int8')

    y0 = 0
    Fahrtig = False

    while not Fahrtig:
        y1 = min(y0 + on, NY)
        mask_onoff[:, y0:y1] = 1
        y0 = y1 + off
        if (y0 > NY - 1):
            break

    #from astropy.io import fits as fts
    # fts.writeto('mask.fits',mask_onoff.transpose().astype('int32'),overwrite=True)

    for Q in ccdobj.Quads:
        qdata = ccdobj.get_quad(Q, canonical=True, extension=extension)
        qdata[np.where(mask_onoff == 0)] = 0.
        ccdobj.set_quad(qdata, Q, canonical=True, extension=extension)
