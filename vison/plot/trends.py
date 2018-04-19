#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Plotting classes shared across tasks/sub-tasks and derived from plots.baseclasses.
They have in common that they show trends with time of some variables / stats.

Created on Fri Jan 26 16:18:43 2018

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy

from vison.plot import figclasses
from vison.plot import baseplotclasses
# END IMPORT


class Fig_Basic_Checkstat(figclasses.Fig_XY_fromDD):

    def __init__(self, *args, **kwargs):
        super(Fig_Basic_Checkstat, self).__init__(*args,**kwargs)
        #self.figname = figname

        self.stats = []
        self.trendaxis = ''
        self.suptitle = ''
        self.caption = ''
        self.texfraction = 1.1
        self.data = dict()
        self.plotclass = baseplotclasses.BeamPlotYvX
        
    def build_data(self,parent,**kwargs):
        """ """
        
        indicesList = ['CCD','Quad']
        
        dd = parent.dd
        indices = parent.dd.indices
        indicesDict = dict()
        for key in indicesList:
            indicesDict[key] = indices[indices.names.index(key)].vals        
        
        ddmx = dict()
        ddmx[self.trendaxis] = copy.deepcopy(dd[self.trendaxis][:])
        for label in self.stats:
            ddmx[label] = copy.deepcopy(dd[label][:])
        
        self.data = self._build_data_yvsx(ddmx,indicesDict,indicesList,self.stats,
                                     self.trendaxis)
        
    def configure(self, **kwargs):
        """ """
        defaults = dict(path='./', stats=[], trendaxis='time')
        defaults.update(kwargs)
        self.stats = defaults['stats']
        self.trendaxis = defaults['trendaxis']
        
        super(Fig_Basic_Checkstat,self).configure(**defaults)
        
    def plot(self, **kwargs):
        """ """
        meta = dict(suptitle=self.suptitle,
                    doNiceXDate=True, doLegend=True)
        meta.update(kwargs)
        
        super(Fig_Basic_Checkstat,self).plot(**meta)


