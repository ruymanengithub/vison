#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Plotting classes shared across tasks/sub-tasks and derived from plots.baseclasses.
They have in common that they show trends with time of some variables / stats.

Created on Fri Jan 26 16:18:43 2018

@author: raf
"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from collections import OrderedDict

import baseclasses
# END IMPORT


class pl_basic_checkstat(baseclasses.Fig):

    def __init__(self):
        super(pl_basic_checkstat, self).__init__()
        #self.figname = figname

        self.stats = []
        self.trendaxis = ''
        self.suptitle = ''
        self.caption = ''
        self.texfraction = 1.1
        self.data = dict()

    def configure(self, **kwargs):
        """ """
        defaults = dict(path='./', stats=[], trendaxis='time')
        defaults.update(kwargs)
        self.stats = defaults['stats']
        self.trendaxis = defaults['trendaxis']
        if 'figname' in defaults:
            self.figname = defaults['figname']
        if 'caption' in defaults:
            self.caption = defaults['caption']
        path = defaults['path']
        self.figname = os.path.join(path, self.figname)
        if 'suptitle' in defaults:
            self.suptitle = defaults['suptitle']

    def build_data(self, parent):
        """ """

        dd = parent.dd
        indices = parent.dd.indices
        CCDs = indices[indices.names.index('CCD')].vals
        Quads = indices[indices.names.index('Quad')].vals

        data = OrderedDict()
        for ixCCD, CCD in enumerate(CCDs):
            CCDkey = 'CCD%i' % CCD
            data[CCDkey] = OrderedDict()
            for iQ, Q in enumerate(Quads):
                data[CCDkey][Q] = OrderedDict()
                data[CCDkey][Q]['x'] = OrderedDict()
                data[CCDkey][Q]['y'] = OrderedDict()

                for stat in self.stats:
                    data[CCDkey][Q]['x'][stat] = dd.mx[self.trendaxis][:, ixCCD].copy()
                    data[CCDkey][Q]['y'][stat] = dd.mx[stat][:, ixCCD, iQ].copy()

        self.data = data.copy()

    def plot(self, **kwargs):
        """ """
        meta = dict(suptitle=self.suptitle,
                    doNiceXDate=True, doLegend=True)
        meta.update(kwargs)

        plotobj = baseclasses.Beam2DPlot(self.data, **meta)
        plotobj.render(self.figname)
