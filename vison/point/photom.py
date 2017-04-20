# -*- coding: utf-8 -*-
"""

Doing Aperture Photometry of an object
======================================

Simple class to do aperture photometry on a stamp of a point-source.

:requires: NumPy

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk


Created on Thu Apr 20 14:37:46 2017

:author: Ruyman Azzollini
:contact: r.azzollini _at_ ucl.ac.uk
"""

# IMPORT STUFF
import numpy as np
from vison.analysis.ellipse import dist_superellipse

from basis import SpotBase
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

    def __init__(self, data, log=None, **kwargs):
        """
        :param data: stamp to be analysed.
        :type data: ndarray
        :param log: logger
        :type log: instance
        :param kwargs: additional keyword arguments
        :type kwargs: dict

        Settings dictionary contains all parameter values needed.
        """
        
        super(Photometer,self).__init__(data,log)
        
        self.phsettings = dict()
        self.phsettings.update(kwargs)
        
        self._gen_rad_mask()
        

        for key, value in self.phsettings.iteritems():
            if self.log is not None: self.log.info('%s = %s' % (key, value))
                
    
    def _gen_rad_mask(self):
        """ """
        n = self.NX,self.NY
        x = self.xcen
        y = self.ycen
        
        self.radmask = dist_superellipse(n,(x,y),q=1,pos_ang=0.,c=0.)
    
    
    def measure_bgd(self,rin,rout):
        """ """
        bgd_mask = ((self.radmask >= rin) & (self.radmask <= rout))
        bgd = np.median(self.data[bgd_mask])
        sbgd = np.std(self.data[bgd_mask])
        
        return bgd,sbgd
    
    def sub_bgd(self,rin,rout):
        """ """
        bgd,sbgd = self.measure_bgd(rin,rout,dostd=False)
        self.data -= bgd
        return bgd,sbgd
    
    def doap_photom(self,centre,rap,rin=-1.,rout=-1.,subbgd=False):
        """ """
        
        x,y = centre
        
        if subbgd:
            _ = self.sub_bgd(rin,rout)
        
        apmask = self._get_circ_mask(x,y,rap)
        
        apflu = np.sum(self.data[apmask])
        sapflu = np.std(self.data[apmask])
        
        return apflu,sapflu
    
    
    def _get_circ_mask(self,x,y,radius):
        """ """
        mask = self.radmask <= radius        
        return mask
    
    def _get_ring_mask(self,x,y,rin,rout):
        """ """
        mask = (self.radmask >= rin) & (self.radmask <= rout)        
        return mask