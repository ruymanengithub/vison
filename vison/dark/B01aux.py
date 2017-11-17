#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Auxiliary Functions and resources to BIAS01.

Created on Tue Nov 14 13:54:34 2017

:author: Ruyman Azzollini
:contact: r.azzollini__at__ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os

from vison.plot.classes import Fig, CCD2DPlot
# END IMPORT


class plB01offsets(Fig):
    
    def __init__(self):
        super(plB01offsets,self).__init__()
        
        self.CCD = 0
        self.figname = 'BIAS01_offset_vs_time_CCD%i.png' 
        self.caption = 'CCD%i: Some Caption' 
        self.texfraction = 0.8
        self.data = dict()
    
    def configure(self,**kwargs):
        """ """
        
        defaults = dict(path='./',CCD=0)
        defaults.update(kwargs)
        
        self.CCD = kwargs['CCD']
        self.figname = self.figname % self.CCD
        self.caption = self.caption % self.CCD
        
        path = kwargs['path']
        self.figname = os.path.join(path,self.figname)
        
    def build_data(self,parent):
        """ """
        
        dd = parent.dd
        
        indices = parent.dd.indices
        CCDs = indices[indices.names.index('CCD')].vals
        Quads = indices[indices.names.index('Quad')].vals
        
        data = dict()
        for Q in ['E','F','G','H']:
            ixQ = Quads.index(Q)
            ixCCD = CCDs.index(self.CCD)
            data[Q] = dict()
            data[Q]['x'] = dict()
            data[Q]['y'] = dict()
            for sec in ['pre','img','ove']:
                data[Q]['x'][sec] = dd.mx['time'][:,ixCCD].copy()
                data[Q]['y'][sec] = dd.mx['offset_%s' % sec][:,ixCCD,ixQ].copy()

        self.data = data
        
        
    def plot(self,**kwargs):
        """ """
        meta = dict(suptitle='CCD %i' % self.CCD,doNiceXDate=True)
        plotobj = CCD2DPlot(self.data,meta=meta)
        plotobj.render(self.figname)


B01figs = dict()
for CCD in [1,2,3]: B01figs['B01offsets_CCD%i' % CCD] = plB01offsets