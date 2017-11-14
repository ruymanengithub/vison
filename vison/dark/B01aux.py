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

from vison.plot.classes import Fig, CCD2DPlot
# END IMPORT


class B01offsets(Fig):
    
    def __init__(self,CCD):
        super(B01offsets,self).__init__()
        self.CCD = CCD
        self.figname = 'BIAS01_offset_vs_time_CCD%i.png' % self.CCD
        self.caption = 'CCD%i: Some Caption' % self.CCD
        self.texfraction = 0.8
        
    def plot(self,parent,**kwargs):
        
        figname = kwargs['figname']
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
        
        meta = dict(suptitle='CCD %i' % self.CCD,doNiceXDate=True)
        plotobj = CCD2DPlot(data,meta=meta)
        plotobj.render(figname)
        self.figname = figname
    

B01figs = dict(B01offsets_CCD1=B01offsets(1),
                      B01offsets_CCD2=B01offsets(2),
                      B01offsets_CCD3=B01offsets(3))