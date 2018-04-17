#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 16:17:13 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import os
import operator
import numpy as np
from collections import defaultdict
from itertools import product
from functools import reduce
import baseplotclasses

from matplotlib import pyplot as plt
# END IMPORT

def getFromDict(dataDict,mapList):
    return reduce(operator.getitem,mapList,dataDict)

def setInDict(dataDict,mapList,value):
    getFromDict(dataDict,mapList[:-1])[mapList[-1]] = value


class Fig(object):
    
    plotclass = None

    def __init__(self, figname=''):

        self.figname = figname
        self.texfraction = 1.0
        self.caption = ''
        self.suptitle = ''
        self.data = dict()

    def configure(self, **kwargs):
        """ """
        defaults = dict(suptitle='')
        defaults.update(kwargs)
        if 'figname' in defaults:
            self.figname = defaults['figname']
        if 'caption' in defaults:
            self.caption = defaults['caption']
        path = defaults['path']
        self.figname = os.path.join(path, self.figname)
        if 'suptitle' in defaults:
            self.suptitle = defaults['suptitle']

    def plot(self,**kwargs):
        """ """
        plotobj = self.plotclass(self.data,**kwargs)
        plotobj.render(self.figname)

class Fig_Beam2DPlot(Fig):
    plotclass = baseplotclasses.BeamPlotXvY

class Fig_Beam1DHist(Fig):
    plotclass = baseplotclasses.Beam1DHist
    

class Fig_XY_fromDD(Fig):

    
    def _build_data_xvsy(self,ddmx,indicesDict,indicesList,labels,trendaxis):
        
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
        
        # on-the-fly creation and assignation of values
        
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
        
        # Convert defaultdict 'back' to dict()
        
        for indices_tuple in indices_tuples:
            for j in range(1,len(indicesList)+2):
                value = getFromDict(data,indices_tuple[:-j])
                if isinstance(value,defaultdict):
                    setInDict(data,indices_tuple[:-j],dict(value))
        data = dict(data)
        
        return data

class BlueScreen(Fig):
    """ """

    def __init__(self):
        from vison import data as vdata

        super(BlueScreen, self).__init__()

        self.rootfigure = os.path.join(
            vdata.__path__[0], 'BlueScreenErrorWindows.png')
        self.figname = 'BlueScreen%s.png'
        self.caption = ''
        self.texfraction = 0.7
        self.plotclass = baseplotclasses.ImgShow

    def configure(self, **kwargs):
        """ """

        defaults = dict(path='./',
                        caption='', tag='',
                        meta=dict(title=''))
        defaults.update(kwargs)
        self.caption = defaults['caption']
        path = defaults['path']
        self.figname = os.path.join(path, self.figname % defaults['tag'])

    def build_data(self, parent):
        """ """
        self.data = plt.imread(self.rootfigure)

