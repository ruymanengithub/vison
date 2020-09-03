# -*- coding: utf-8 -*-
"""

Aperture Photometry of point-like objects
=========================================

Simple class to do aperture photometry on a stamp of a point-source.

:requires: NumPy

Created on Thu Apr 20 14:37:46 2017


:author: Ruyman Azzollini


"""

# IMPORT STUFF
from pdb import set_trace as stop
from scipy import ndimage as nd

import numpy as np
from vison.analysis.ellipse import dist_superellipse
from vison.point import algorithms as poalg
from vison.point.basis import SpotBase
# END IMPORT


class Photometer(SpotBase):
    """
    Provides methods to measure the shape of an object.

    :param data: stamp to be analysed.
    :type data: np.ndarray
    :param log: logger
    :type log: instance
    :param kwargs: additional keyword arguments
    :type kwargs: dict

    Settings dictionary contains all parameter values needed.
    """

    def __init__(self, data, log=None, verbose=False, **kwargs):
        """
        :param data: stamp to be analysed.
        :type data: ndarray
        :param log: logger
        :type log: instance
        :param verbose: bool switch
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        Settings dictionary contains all parameter values needed.
        """

        super(Photometer, self).__init__(data, log, verbose)

        self.phsettings = dict()
        self.phsettings.update(kwargs)

        self._gen_rad_mask()

        for key, value in list(self.phsettings.items()):
            if self.log is not None:
                self.log.info('%s = %s' % (key, value))

    def _gen_rad_mask(self):
        """ """
        n = self.NX, self.NY
        x = self.xcen
        y = self.ycen

        self.radmask = dist_superellipse(n, (x, y), q=1, pos_ang=0., c=0.)

    def measure_bgd(self, rin, rout):
        """ """

        _statops = dict(masked=dict(
            median=np.ma.median,
            std=np.ma.std),
            nonmasked=dict(
                median=np.nanmedian,
                std=np.nanstd))

        bgd_mask = ((self.radmask >= rin) & (self.radmask <= rout))

        if np.ma.isMaskedArray(self.data):
            _oper = _statops['masked']
        else:
            _oper = _statops['nonmasked']

        bgd = _oper['median'](self.data[bgd_mask])
        sbgd = _oper['std'](self.data[bgd_mask])

        return bgd, sbgd

    def sub_bgd(self, rin, rout):
        """ """
        bgd, sbgd = self.measure_bgd(rin, rout)
        self.data -= bgd
        return bgd, sbgd

    def get_centroid(self, rap=None, full=False):
        """
        TODO:
            add aperture masking
        """
        #xcen0 = self.data.shape[1]/2
        #ycen0 = self.data.shape[0]/2
        medimg = nd.median_filter(self.data, [5, 5])
        ycen0, xcen0 = np.unravel_index(np.argmax(medimg), medimg.shape)
        xcen, ycen, xcen2, ycen2, xcen3, ycen3 = poalg.fwcentroid(self.data, checkbox=1, 
                    maxiterations=10,
                    threshold=1e-3, halfwidth=6, verbose=False,
                    full=full, CEN0=(xcen0, ycen0))
        self.xcen = xcen
        self.ycen = ycen

        fwhmx = np.sqrt(xcen2 - xcen**2.) * 2.355  # Gaussian Approx.
        fwhmy = np.sqrt(ycen2 - ycen**2.) * 2.355  # Gaussian Approx.

        if full:
            return xcen, ycen, fwhmx, fwhmy
        return xcen, ycen

    def doap_photom(self, centre, rap, rin=-1., rout=-1., gain=3.5, doErrors=True,
                    subbgd=False):
        """ """

        img = self.data.copy()
        x, y = centre

        if doErrors:
            ring_mask = self._get_ring_mask(x, y, rin, rout)
            Nbgd = len(np.where(ring_mask)[0])
            bgd, sbgd = self.measure_bgd(rin, rout)

        if subbgd:
            if not doErrors:
                bgd, sbgd = self.measure_bgd(rin, rout)
            img -= bgd
        else:
            bgd = 0.

        apmask = self._get_circ_mask(x, y, rap)
        Nap = len(np.where(apmask)[0])

        apflu = np.sum(self.data[apmask]) - bgd
        sapflu = np.std(self.data[apmask])

        # safebag = dict(img=img,
        #               ampmask=apmask,
        #               rin=rin,
        #               rout=rout,
        #               centre=centre,
        #               bgd=bgd,
        #               ring_mask=ring_mask,
        #               apflu=apflu,
        #               sapflu=sapflu)

        if doErrors:
            ebgd = (sbgd * gain) / np.sqrt(Nbgd)
            sapflu = np.sqrt(apflu * gain + (sbgd * gain)**2. * Nap + ebgd**2.) / gain
            return apflu, sapflu  # , safebag

        return apflu

    def _get_circ_mask(self, x, y, radius):
        """ """
        mask = self.radmask <= radius
        return mask

    def _get_ring_mask(self, x, y, rin, rout):
        """ """
        mask = (self.radmask >= rin) & (self.radmask <= rout)
        return mask
