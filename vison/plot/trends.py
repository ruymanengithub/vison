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
import baseclasses
# END IMPORT

class pl_checkstat(baseclasses.Fig):
    
    def __init__(self,stats,figname='',suptitle='',caption=''):
        super(pl_checkstat,self).__init__()
        
        self.figname = figname
        self.caption = caption
        self.texfraction = 1.1
        self.stats = stats
        self.suptitle = suptitle
        self.data = dict()
    
#    def configure(self,**kwargs):
#        """ """
#        defaults = dict(path='./',stat='offset')
#        defaults.update(kwargs)
#        self.stat = defaults['stat']
#        self.figname = self.figname % self.stat
#        self.caption = self.caption % self.stat
#        path = defaults['path']
#        self.figname = os.path.join(path,self.figname)
#        self.suptitle = 'BIAS01-checks: %s' % self.stat
#        
#    def build_data(self,parent):
#        """ """
#        
#        dd = parent.dd
#        indices = parent.dd.indices
#        CCDs = indices[indices.names.index('CCD')].vals
#        Quads = indices[indices.names.index('Quad')].vals
#        
#        data = OrderedDict()
#        for ixCCD,CCD in enumerate(CCDs):
#            CCDkey = 'CCD%i' % CCD
#            data[CCDkey] = OrderedDict()
#            for iQ,Q in enumerate(Quads):
#                data[CCDkey][Q] = OrderedDict()
#                data[CCDkey][Q]['x'] = OrderedDict()
#                data[CCDkey][Q]['y'] = OrderedDict()
#                for sec in ['pre','img','ove']:
#                    data[CCDkey][Q]['x'][sec] = dd.mx['time'][:,ixCCD].copy()
#                    #data[CCDkey][Q]['x'][sec] = np.arange(len(dd.mx['time'][:,ixCCD])) # dd.mx['time'][:,ixCCD].copy()
#                    data[CCDkey][Q]['y'][sec] = dd.mx['%s_%s' % (self.stat,sec)][:,ixCCD,iQ].copy()
#        
#        self.data = data
#        
#        
#    def plot(self,**kwargs):
#        """ """
#        meta = dict(suptitle=self.suptitle,
#                    doNiceXDate=True,doLegend=True)
#        kwargs.update(meta)
#        plotobj = plbaseclasses.Beam2DPlot(self.data,**kwargs)
#        plotobj.render(self.figname)
