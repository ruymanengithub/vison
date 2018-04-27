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
import numpy as np

from vison.plot import figclasses
from vison.plot import baseplotclasses
from vison.datamodel.core import vColumn
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
        
    
    def build_data(self,dd,**kwargs):
        """ """
        
        index2margin = 'ix'
        indices2param = ['CCD','Quad']
        stats = self.stats
        
        subdd = dict()
        subdd[self.trendaxis] = copy.deepcopy(dd.mx[self.trendaxis])
        trendaxes = subdd[self.trendaxis].indices.get_names()
        
        assert index2margin in trendaxes
        assert np.all([ixkey in [index2margin]+indices2param for ixkey in trendaxes])
        
        refstatindices = dd.mx[self.stats[0]].indices.get_names() # all stats must have same indices
        indices2tag = [index for index in refstatindices if index not in [index2margin]+indices2param]
        
        assert len(indices2tag)<=1
                
        if len(stats)>1:
            assert len(indices2tag)==0
        
        if len(indices2tag)==1:
            index2tag = indices2tag[0]
        else:
            index2tag = None
        
        
        if len(stats)>1 and index2tag is None:
            
            tags = stats
            
            for tag in tags:
                subdd[tag] = copy.deepcopy(dd.mx[tag])
                _statindices =  subdd[tag].indices.get_names()
                assert _statindices == refstatindices
                
        elif len(stats)==1 and index2tag is not None:
            
            tags = dd.mx[stats[0]].indices.get_vals(index2tag)
            stat = stats[0]
            
            for tag in tags:
                subdd[tag] = vColumn(dd.mx[stat][...,0],dd.mx[stat].name,dd.mx[stat].indices[:-1])
            
        elif len(stats)==1 and index2tag is None:
            
            tags = stats
            subdd[tags[0]] = copy.deepcopy(dd.mx[tags[0]])
        
        else:
            raise AssertionError
        
        
        self.data = self._build_data_yvsx(subdd,tags,
                self.trendaxis,index2margin,indices2param)
        
        
    
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


