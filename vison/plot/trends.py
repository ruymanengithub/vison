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
import os
from collections import OrderedDict
import numpy as np
import copy

from vison.plot import baseclasses
# END IMPORT

from functools import reduce
import operator
from itertools import product
from collections import defaultdict


def getFromDict(dataDict,mapList):
    return reduce(operator.getitem,mapList,dataDict)

def setInDict(dataDict,mapList,value):
    getFromDict(dataDict,mapList[:-1])[mapList[-1]] = value

class pl_basic_checkstat(baseclasses.Fig):

    def __init__(self, *args, **kwargs):
        super(pl_basic_checkstat, self).__init__(*args,**kwargs)
        #self.figname = figname

        self.stats = []
        self.trendaxis = ''
        self.suptitle = ''
        self.caption = ''
        self.texfraction = 1.1
        self.data = dict()
        self.plotclass = baseclasses.Beam2DPlot
        
        
    def _build_data(self,ddmx,indicesDict,indicesList,labels,trendaxis):
        
        plnames = ['x','y']
    
        iterable_keys = []
        for indname in indicesList:
            iterable_keys.append(indicesDict[indname])
        iterable_keys += [plnames] + [labels]
        
        iterable_ix = []
        for indname in indicesList:
            iterable_ix.append(np.arange(len(indicesDict[indname])).tolist())
        iterable_ix += [[None]*len(plnames)] + [[None]*len(labels)]
        
        
        # nested dictionary based on recursion and collections.defaultdict
        nested_dict = lambda: defaultdict(nested_dict)
        data = nested_dict()
        
        indices_tuples = [item for item in product(*iterable_keys)]
        ix_tuples = [item for item in product(*iterable_ix)]
        
        # assignation of values
        
        for i in range(len(indices_tuples)):
            indices_tuple = indices_tuples[i]
            ix_tuple = ix_tuples[i]
            
            axis = indices_tuple[-2]
            stat = indices_tuple[-1]
            
            if axis == 'x':
                value = ddmx[trendaxis][ix_tuple[:-2]]
            elif axis == 'y':
                value = ddmx[stat][ix_tuple[:-2]]
                
            setInDict(data,indices_tuple,value)
        
        # Convert defaultdict to dict
        
        for indices_tuple in indices_tuples:
            for j in range(1,len(indicesList)+2):
                value = getFromDict(data,indices_tuple[:-j])
                if isinstance(value,defaultdict):
                    setInDict(data,indices_tuple[:-j],dict(value))
        data = dict(data)
        
        return data
    
    def build_data(self,parent):
        """ """
        dd = parent.dd
        indices = parent.dd.indices
        indicesList = ['CCD','Quad']
        indicesDict = dict()
        for key in indicesList:
            indicesDict[key] = indices[indices.names.index(key)].vals        
        
        ddmx = dict()
        ddmx[self.trendaxis] = copy.deepcopy(dd[self.trendaxis][:])
        for label in self.stats:
            ddmx[label] = copy.deepcopy(dd[label][:])
        
        self.data = self._build_data(ddmx,indicesDict,indicesList,self.stats,
                                     self.trendaxis)
        
#    def build_data(self, parent):
#        """ """
#
#        dd = parent.dd
#        indices = parent.dd.indices
#        CCDs = indices[indices.names.index('CCD')].vals
#        Quads = indices[indices.names.index('Quad')].vals
#
#        data = OrderedDict()
#        for ixCCD, CCD in enumerate(CCDs):
#            CCDkey = 'CCD%i' % CCD
#            data[CCDkey] = OrderedDict()
#            for iQ, Q in enumerate(Quads):
#                data[CCDkey][Q] = OrderedDict()
#                data[CCDkey][Q]['x'] = OrderedDict()
#                data[CCDkey][Q]['y'] = OrderedDict()
#                for stat in self.stats:
#                    data[CCDkey][Q]['x'][stat] = dd.mx[self.trendaxis][:, ixCCD].copy()
#                    data[CCDkey][Q]['y'][stat] = dd.mx[stat][:, ixCCD, iQ].copy()
#        #print self.stats,self.trendaxis        
#        self.data = data.copy()

    def configure(self, **kwargs):
        """ """
        defaults = dict(path='./', stats=[], trendaxis='time')
        defaults.update(kwargs)
        self.stats = defaults['stats']
        self.trendaxis = defaults['trendaxis']
        
        super(pl_basic_checkstat,self).configure(**defaults)
        
    def plot(self, **kwargs):
        """ """
        meta = dict(suptitle=self.suptitle,
                    doNiceXDate=True, doLegend=True)
        meta.update(kwargs)
        
        super(pl_basic_checkstat,self).plot(**meta)
