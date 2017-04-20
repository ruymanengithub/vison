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
# END IMPORT


class Photometer():
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
        self.data = data.copy()
        self.log = log

        NY, NX = self.data.shape
        self.NX = NX
        self.NY = NY
        
        
        self.phsettings = dict(x=NX/2,y=NY/2)

        self.settings.update(kwargs)
        
        self._gen_rad_mask()
        

        for key, value in self.settings.iteritems():
            if self.log is not None: self.log.info('%s = %s' % (key, value))
                
    
    def _gen_rad_mask(self):
        """ """
        n = self.NX,self.NY
        x = self.phsettings['x']
        y = self.phsettings['y']
        
        self.radmask = dist_superellipse(n,(x,y),q=1,pos_ang=0.,c=0.)
    
    
    def measure_bgd(self,rin,rout):
        """ """
        bgd_mask = ((self.radmask >= rin) & (self.radmask <= rout))
        bgd = np.sum(self.data[bgd_mask])
        return bgd
        
    
    def sub_bgd(self):
        """ """
        bgd = self.measure_bgd()
        self.data -= bgd
    
    def doaperture_photom(self,centre,rap,rin=-1.,rout=1.,subbgd=False):
        """ """
        
        img = self.data.copy()
        x,y = centre
        
        if subbgd:
            self.sub_bgd()
        
        apmask = self._get_circ_mask(x,y,rap)
        
        apflu = np.sum(img[apmask])
        
        return apflu
        
    
    def _get_circ_mask(self,x,y,radius):
        """ """
        mask = self.radmask <= radius        
        return mask
    
    def _get_ring_mask(self,x,y,rin,rout):
        """ """
        mask = (self.radmask >= rin) & (self.radmask <= rout)        
        return mask