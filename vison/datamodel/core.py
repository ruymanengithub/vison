#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

DataDict Class : holds data and results across sub-tasks of a "task" (Test).
This is the CORE data-structure used to do analysis and report results.


Created on Thu Sep 21 16:47:09 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
from collections import OrderedDict
from copy import deepcopy
import string as st
import astropy as ast
import itertools
import tempfile

from vison.datamodel import EXPLOGtools as ELtools
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.support import flags
# END IMPORT


class vIndex(object):
    """ """

    def __init__(self, name, vals=None, N=0):
        """ """
        self.name = name
        if vals is None:
            vals = []
        
        if len(vals) != 0:
            self.vals = vals
            self.len = len(self.vals)
        elif (len(vals) == 0) and N != 0:
            self.vals = np.arange(N).tolist()
            self.len = N
        else:
            raise TypeError

    def get_vals(self):
        return self.vals

    def get_len(self):
        return self.len

    def __str__(self):
        return '("%s": %i)' % (self.name, self.len)

    def __repr__(self):
        return '("%s": %i)' % (self.name, self.len)


class vMultiIndex(list, object):

    def __init__(self, IndexList=None):
        
        if IndexList is None:
            IndexList = []

        assert (isinstance(IndexList, list) or isinstance(IndexList, vMultiIndex) or
                (isinstance(IndexList, vIndex)))

        try:
            self += IndexList
        except TypeError:
            self += [IndexList]
        self.names = self.get_names()
        self.shape = self.get_shape()

    def get_names(self):
        """ """
        names = []
        for item in self:
            names.append(item.name)
        return names

    def find(self, indexname):
        return self.names.index(indexname)

    def get_vals(self, indexname):
        """ """
        return self[self.find(indexname)].get_vals()

    def get_len(self, indexname):
        return self[self.find(indexname)].get_len()

    def get_shape(self):
        shape = []
        for item in self:
            shape.append(item.len)
        return tuple(shape)

    def update_names(self):
        self.names = self.get_names()

    def update_shape(self):
        self.shape = self.get_shape()

    def append(self, *args):
        super(vMultiIndex, self).append(*args)
        self.update_names()
        self.update_shape()

    def __add__(self, *args):
        super(vMultiIndex, self).__add__(*args)
        self.update_names()
        self.update_shape()

    def pop(self, *args):
        super(vMultiIndex, self).pop(*args)
        self.update_names()
        self.update_shape()

    def __getslice__(self, i, j):
        return vMultiIndex(super(vMultiIndex, self).__getslice__(i, j))

    def __str__(self):
        inside = st.join(['%s' % item.__str__() for item in self], ',')
        return '[%s]' % inside


class vColumn(object):
    """ """

    def __init__(self, array, name, indices):
        """ """

        self.array = array.copy()
        self.shape = self.array.shape
        assert ((isinstance(indices, list)) or (
            isinstance(indices, vMultiIndex)))
        for index in indices:
            assert isinstance(index, vIndex)

        assert len(self.shape) == len(indices)

        for i in range(len(self.shape)):
            assert self.shape[i] == indices[i].len, stop()

        if isinstance(indices, list):
            self.indices = vMultiIndex(indices)
        else:
            self.indices = indices
        self.name = name

    def name_indices(self):
        """ """
        print self.indices.names

    def __repr__(self):
        return self.array.__repr__()

    def __call__(self):
        return self.array

    def __getslice__(self, i, j):
        return self.array.__getslice__(i, j)

    def __setslice__(self, i, j, y):
        return self.array.__setslice__(i, j, y)

    def __setitem__(self, i, y):
        return self.array.__setitem__(i, y)

    def __getitem__(self, i):
        return self.array.__getitem__(i)

    def __str__(self):
        return '%s: %s' % (self.name, self.array.__str__())


class DataDict(object):
    """ """

    def __init__(self, meta=None):
        """ """
        
        if meta is None:
            meta = OrderedDict()

        self.meta = meta
        self.mx = OrderedDict()
        self.colnames = []
        self.indices = vMultiIndex()
        self.products = dict()  # data products
        self.flags = flags.Flags()
        self.compliances = OrderedDict()

    def loadExpLog(self, explog):
        """ """
        # CCDs = [1,2,3]
        CCDs = ['CCD1', 'CCD2', 'CCD3']

        ObsID = explog['ObsID'].data
        uObsID = np.unique(ObsID)
        Nobs = len(uObsID)

        ccdcol = explog['CCD']

        ObsIndex = vIndex('ix', N=Nobs)
        commIndices = vMultiIndex([ObsIndex, vIndex('CCD', CCDs)])

        self.addColumn(uObsID, 'ObsID', [ObsIndex])

        for key in explog.colnames:

            if key == 'ObsID':
                continue

            arrlist = []
            for CCDkey in CCDs:
                #CCDkey = 'CCD%i' % CCDindex
                arrlist.append(explog[key].data[np.where(ccdcol == CCDkey)])
            array = np.array(arrlist).transpose()
            self.addColumn(array, key, commIndices)

        return None

    def initColumn(self, name, indices, dtype='float32', valini=0.):
        """ """

        assert isinstance(indices, vMultiIndex)

        shape = indices.shape
        if dtype[0] != 'S':
            array = np.zeros(shape, dtype=dtype) + valini
        else:
            array = np.zeros(shape, dtype=dtype)
            array[:] = valini

        if name in self.colnames:
            self.dropColumn(name)
            
        
        self.addColumn(array, name, indices)

    def addColumn(self, array, name, indices, ix=-1):
        """ """

        column = vColumn(array, name, indices)

        self.mx[column.name] = column
        if ix==-1:
            self.colnames.append(column.name)
        else:
            self.colnames.insert(ix,column.name)
        colindices = column.indices

        selfindnames = self.indices.names

        for ic, index in enumerate(colindices):
            if index.name in selfindnames:
                ix = selfindnames.index(index.name)
                assert np.all(index.vals == self.indices[ix].vals)
                assert ic == ix
            else:
                self.indices.append(index)

    def name_indices(self):
        """ """
        print self.indices.names

    def col_has_index(self, colname, indexname):
        """ """
        assert colname in self.colnames
        assert indexname in self.indices.names
        if indexname in self.mx[colname].indices.names:
            return True
        return False

    def dropColumn(self, colname):
        """ """

        assert colname in self.colnames

        colindices = self.mx[colname].indices
        colnames = self.colnames

        indexescounter = []
        for colindex in colindices:
            indexescounter.append(len([1 for icolname in colnames if
                                       colindex.name in self.mx[icolname].indices.names]))

        ixunique = np.where(np.array(indexescounter) == 1)
        for ixix in ixunique[0]:
            self.indices.pop(ixix)

        self.mx.pop(colname)
        self.colnames.pop(self.colnames.index(colname))
    
    
    def flattentoTable(self):
        """ """

        t = ast.table.Table()

        for col in self.mx:
            cindex = self.mx[col].indices
            carray = self.mx[col]().copy()
            if len(cindex) == 1:
                t[col] = ast.table.Column(carray)
            elif len(cindex) > 1:
                _names = cindex[1:].names
                _vals = [item.vals for item in cindex[1:]]
                _ixix = [range(item.len) for item in cindex[1:]]
                prod_vals = list(itertools.product(*_vals))
                prod_ixix = list(itertools.product(*_ixix))
                tmp_suf = st.join([('_%s' % item)+'%s' for item in _names], '')

                nix0 = cindex[0].len

                for i in range(len(prod_vals)):
                    suf = tmp_suf % prod_vals[i]
                    key = col + suf
                    jxx = [tuple([jx] + list(prod_ixix[i]))
                           for jx in range(nix0)]

                    carray = np.array([self.mx[col][_jxx] for _jxx in jxx])

                    t[key] = ast.table.Column(carray)
        
        return t
        

    def saveToFile(self, outfile, format='ascii.commented_header'):
        """ """
        
        t = self.flattentoTable()
        t.write(outfile, format=format, overwrite=True)

    def __cmp__(self, other):
        return self.__dict__ == other.__dict__


def useCases():
    """ 

    #TODO:

        # create a DataDict object from an exposure log.
        # add a column indexed by ObsID, CCD and Quad
        # drop a column
        # create a column from an operation on several columns with different dimensions
        # save to a text / excel file
        # save to a pickle file

    """

    dpath = '/home/raf/WORK/EUCLID/CALIBRATION/PipelineDevel/TEST_DATA/24_Feb_80/'
    explogf = os.path.join(dpath, 'EXP_LOG_240280.txt')
    elvis = '6.5.X'
    explog = ELtools.loadExpLog(explogf, elvis=elvis)

    dd = DataDict()

    dd.loadExpLog(explog)

    print 'dd.indices'
    print dd.indices

    ans1 = dd.col_has_index('ObsID', 'ix')
    print 'ix in ObsID-col: %s' % ans1
    ans2 = dd.col_has_index('ObsID', 'CCD')
    print 'CCD in ObsID-col: %s' % ans2

    Xindices = deepcopy(dd.indices)
    Xindices.append(vIndex('Quad', vals=['E', 'F', 'G', 'H']))

    ArrShape = []
    for index in Xindices:
        ArrShape.append(index.len)
    ArrShape = tuple(ArrShape)

    print 'Adding spot_fluence'

    spot_fluence = np.zeros(ArrShape, dtype='float32')

    dd.addColumn(spot_fluence, 'spot_fluence', Xindices)

    print 'dd indices = '
    print dd.indices

    ans3 = dd.col_has_index('spot_fluence', 'Quad')
    print 'Quad in spot_fluence: %s' % ans3

    def all_checker(key): return dd.col_has_index(key, 'Quad')

    ans4 = np.any(map(all_checker, [
                  colname for colname in dd.colnames if colname != 'spot_fluence']))
    print 'Quad in colnames != spot_fluence: %s' % ans4

    print 'Dropping spot_fluence'
    dd.dropColumn('spot_fluence')

    print 'Columns: ', dd.colnames
    print 'Indices: ', dd.indices

    OCQindices = Xindices

    OCindices = Xindices[0:2]

    dd.initColumn('gain', OCindices, dtype='float32', valini=1.)
    dd.mx['gain'][:, 0] = 1.
    dd.mx['gain'][:, 1] = 2.
    dd.mx['gain'][:, 2] = 3.

    dd.initColumn('spot_fwhm', OCQindices, dtype='float32', valini=0.)

    nObs = dd.indices.get_len('ix')
    nCCD = dd.indices.get_len('CCD')
    nQuad = dd.indices.get_len('Quad')

    for iObs in range(nObs):
        for jCCD in range(nCCD):
            for kQ in range(nQuad):

                dd.mx['spot_fwhm'][iObs, jCCD, kQ] = (dd.mx['ObsID'][iObs] +
                                                      dd.mx['gain'][iObs, jCCD]*1.e4) * (-1.)**kQ

    # SAVING
    cwd = os.path.abspath(os.path.curdir)

    # PICKLING

    f1 = tempfile.NamedTemporaryFile(
        mode='w+b', prefix='tmp', dir=cwd, delete=True)

    with f1:

        pickf = f1.name

        data = dict(dd=dd)
        cPickleDumpDictionary(data, output=pickf)

        res = cPickleRead(pickf)

        f1.close()

    ddret = res['dd']

    print 'Pickled == Original : %s' % (dd == ddret,)

    outf = 'test_DD.txt'

    f2 = tempfile.NamedTemporaryFile(mode='w+b', prefix='tmp', dir=cwd,
                                     delete=True)
    with f2:
        try:
            dd.saveToFile(outf, format='ascii')
            print 'saved catalog to temp file: %s' % f2.name
            print 'erased %s...' % f2.name
        except:
            print 'failed to save catalog to temp file'


if __name__ == '__main__':

    useCases()
