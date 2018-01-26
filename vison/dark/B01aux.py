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
from collections import OrderedDict

from vison.plot import baseclasses as plbaseclasses

# END IMPORT

class plB01check(plbaseclasses.Fig):
    
    def __init__(self):
        super(plB01check,self).__init__()
        
        self.figname = 'BIAS01_%s_vs_time.png' 
        self.caption = 'BIAS01: %s vs. time.' 
        self.texfraction = 1.1
        self.stat = ''
        self.suptitle = ''
        self.data = dict()
    
    def configure(self,**kwargs):
        """ """
        defaults = dict(path='./',stat='offset')
        defaults.update(kwargs)
        self.stat = defaults['stat']
        self.figname = self.figname % self.stat
        self.caption = self.caption % self.stat
        path = defaults['path']
        self.figname = os.path.join(path,self.figname)
        self.suptitle = 'BIAS01-checks: %s' % self.stat
        
    def build_data(self,parent):
        """ """
        
        dd = parent.dd
        indices = parent.dd.indices
        CCDs = indices[indices.names.index('CCD')].vals
        Quads = indices[indices.names.index('Quad')].vals
        
        data = OrderedDict()
        for ixCCD,CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            data[CCDkey] = OrderedDict()
            for iQ,Q in enumerate(Quads):
                data[CCDkey][Q] = OrderedDict()
                data[CCDkey][Q]['x'] = OrderedDict()
                data[CCDkey][Q]['y'] = OrderedDict()
                for sec in ['pre','img','ove']:
                    data[CCDkey][Q]['x'][sec] = dd.mx['time'][:,ixCCD].copy()
                    #data[CCDkey][Q]['x'][sec] = np.arange(len(dd.mx['time'][:,ixCCD])) # dd.mx['time'][:,ixCCD].copy()
                    data[CCDkey][Q]['y'][sec] = dd.mx['%s_%s' % (self.stat,sec)][:,ixCCD,iQ].copy()
        
        self.data = data
        
        
    def plot(self,**kwargs):
        """ """
        meta = dict(suptitle=self.suptitle,
                    doNiceXDate=True,doLegend=True)
        kwargs.update(meta)
        plotobj = plbaseclasses.Beam2DPlot(self.data,**kwargs)
        plotobj.render(self.figname)
        

B01figs = dict()
B01figs['B01checks_offsets'] = [plB01check,dict(stat='offset')]
B01figs['B01checks_stds'] = [plB01check,dict(stat='std')]
B01figs['BlueScreen'] = [plbaseclasses.BlueScreen,dict()]