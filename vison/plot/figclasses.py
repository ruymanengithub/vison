#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 16:17:13 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
from collections import defaultdict
from itertools import product
from vison.plot import baseplotclasses
from vison.support import utils

from matplotlib import pyplot as plt
# END IMPORT

# def getFromDict(dataDict,mapList):
#    return reduce(operator.getitem,mapList,dataDict)
#
# def setInDict(dataDict,mapList,value):
#    getFromDict(dataDict,mapList[:-1])[mapList[-1]] = value
#


class Fig(object):

    def plotclass(x): return None

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

    def plot(self, **kwargs):
        """ """
        plotobj = self.plotclass(self.data, **kwargs)
        plotobj.render(self.figname)

    def build_data(self, *args, **kwargs):
        raise NotImplementedError("child class implements abstract method")


class Fig_Beam2DPlot(Fig):
    plotclass = baseplotclasses.BeamPlotYvX


class Fig_Beam1DHist(Fig):
    plotclass = baseplotclasses.Beam1DHist

class Fig_BeamImgShow(Fig):
    plotclass = baseplotclasses.BeamImgShow


class Fig_XY_fromDD(Fig):

    def _build_data_yvsx(self, subdd, stats, trendaxis, index2margin, indices2param):

        assert isinstance(subdd, dict)
        assert isinstance(stats, list)
        assert isinstance(trendaxis, str)


        iterable_data_keys = []
        for indname in indices2param:
            iterable_data_keys.append(
                subdd[stats[0]].indices.get_vals(indname))


        data_keys_tuples = [item for item in product(*iterable_data_keys)]

        trend_indices = subdd[trendaxis].indices
        # all stats SHOULD have same indices
        stats_indices = subdd[stats[0]].indices

        def get_ix_tuples(indobj, indexlist, margin=None):

            iterable_ix = []
            if margin is not None:
                iterable_ix += [[slice(None)]]

            for _ixname in indexlist:
                iterable_ix.append(np.arange(indobj.get_len(_ixname)))
            return [item for item in product(*iterable_ix)]

        ix_tuples = get_ix_tuples(
            stats_indices, indices2param, margin=index2margin)

        # nested dictionary based on recursion and collections.defaultdict
        def nested_dict(): return defaultdict(nested_dict)
        data = nested_dict()

        NdimsTrend = len(trend_indices.shape)
        NdimsStats = len(stats_indices.shape)

        # broadcasting trends axis matrix, if needed

        sup_trendmx = subdd[trendaxis][:].copy()

        def add_axis(arr, repeats):
            return arr[..., np.newaxis].repeat(repeats, axis=-1).copy()

        if NdimsStats > NdimsTrend:
            stats_shape = stats_indices.shape
            for i in stats_shape[NdimsTrend:]:
                sup_trendmx = add_axis(sup_trendmx, i)

        # on-the-fly creation and assignation of values

        assert len(data_keys_tuples) == len(ix_tuples)

        full_keys_tuples = []

        for i in range(len(data_keys_tuples)):

            keys_tuple = data_keys_tuples[i]
            ix_tuple = ix_tuples[i]

            for stat in stats:
                _keysx = tuple(list(keys_tuple)+['x']+[stat])

                valuex = sup_trendmx[ix_tuple].copy()

                utils.setInDict(data, _keysx, valuex)

                _keysy = tuple(list(keys_tuple)+['y']+[stat])

                valuey = subdd[stat][ix_tuple].copy()
                utils.setInDict(data, _keysy, valuey)

                full_keys_tuples += [_keysx, _keysy]

        # Convert defaultdict 'back' to dict()

        Nlevels = len(full_keys_tuples[0])

        for keys_tuple in full_keys_tuples:
            for j in range(1, Nlevels):
                value = utils.getFromDict(data, keys_tuple[:-j])
                if isinstance(value, defaultdict):
                    utils.setInDict(data, keys_tuple[:-j], dict(value))
        data = dict(data)
        
        data['labelkeys'] = stats

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

    def build_data(self, *args, **kwargs):
        """ """
        self.data = plt.imread(self.rootfigure)
