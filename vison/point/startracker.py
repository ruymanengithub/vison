#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 16 11:12:02 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import unittest
from collections import OrderedDict
import numpy as np
#import collections
from astropy.io import ascii
from astropy.table import Table
import astroalign as aa
import skimage
import copy
import os

from vison.datamodel import ccd as ccdmod
from vison import ogse_profiles
# END IMPORT


_path = os.path.abspath(ogse_profiles.__file__)
ogsepath = os.path.dirname(_path)

default_patt_file = 'Pattern_CHAMBER_A_JUN18.txt'


starnames = ['ALPHA', 'BRAVO', 'CHARLIE', 'DELTA', 'ECHO']


class StarTracker(object):

    def __init__(self, CCD, withpover=True, patt_file=None):

        self.Quads = ['E', 'F', 'G', 'H']
        self.starnames = starnames
        if patt_file is None:
            self.patt_file = copy.deepcopy(default_patt_file)
        else:
            self.patt_file = copy.deepcopy(patt_file)

        self.CCD = CCD
        self.withpover = withpover  # parallel overscans

        self.cobj = ccdmod.CCD(withpover=self.withpover)

        self.Pattern = OrderedDict(x=np.array([]), y=np.array([]),
                                   ID=np.array([]))

        self.init_Pattern(CCD)

    def _coodict_2_nested(self, coodict):
        nested = OrderedDict()
        for Q in self.Quads:
            nested[Q] = OrderedDict()
            for star in self.starnames:
                ID = '%s%s' % (Q, star[0])
                ixID = np.where(coodict['ID'] == ID)
                nested[Q][star] = \
                    (coodict['X'][ixID][0], coodict['Y'][ixID][0])
        return nested

    def _nested_2_coodict(self, nested):
        """ """
        coodict = OrderedDict()
        X = []
        Y = []
        ID = []
        for Q in self.Quads:
            for star in self.starnames:
                X.append(nested[Q][star][0])
                Y.append(nested[Q][star][1])
                ID.append('%s%s' % (Q, star[0]))
        coodict['X'] = np.array(X).copy()
        coodict['Y'] = np.array(Y).copy()
        coodict['ID'] = np.array(ID).copy()

        return coodict

    def init_Pattern(self, CCD):
        """ """
        pattfile = self.patt_files['CCD%i' % CCD]
        self.load_Patt_fromfile(pattfile)

    def _load_gencoos_fromfile(self, coosfile):
        """ """
        if coosfile is None:
            return

        coodict = OrderedDict()
        coostab = ascii.read(coosfile)
        coodict['ID'] = coostab['ID'].data.copy()
        coodict['X'] = coostab['X'].data.copy()
        coodict['Y'] = coostab['Y'].data.copy()

        return coodict

    def _save_gencoos_tofile(self, datadict, outfile, overwrite=False):
        ID = datadict['ID'].copy()
        X = datadict['X'].copy()
        Y = datadict['Y'].copy()
        t = Table([ID, X, Y], names=('ID', 'X', 'Y'))
        ascii.write(t, output=outfile, format='commented_header',
                    delimiter='\t',
                    formats=dict(X='%.2f', Y='%.2f'),
                    overwrite=overwrite)

    def load_Patt_fromdict(self, Pattdict):
        assert isinstance(Pattdict, (dict, OrderedDict))
        assert 'ID' in Pattdict
        assert 'X' in Pattdict
        assert 'Y' in Pattdict
        assert len(Pattdict['X']) == len(Pattdict['Y'])
        assert len(Pattdict['X']) == len(Pattdict['ID'])

        self.Pattern = copy.deepcopy(Pattdict)

    def load_Patt_fromfile(self, Pattfile):
        self.Pattern = self._load_gencoos_fromfile(Pattfile)

    def save_Patt_tofile(self, Pattfile, overwrite=False):
        self._save_gencoos_tofile(self.Pattern, Pattfile, overwrite=overwrite)

    def convert_Phys_2_CCD(self, X, Y):
        """ """
        return self.cobj.cooconv_Phys_2_CCD(X, Y)

    def convert_CCD_2_Phys(self, X, Y):
        return self.cobj.cooconv_CCD_2_Phys(X, Y)

    def get_allCCDcoos(self, nested=False):
        """ """
        X = self.Pattern['X'].copy()
        Y = self.Pattern['Y'].copy()

        Xc, Yc = self.convert_Phys_2_CCD(X, Y)

        coodict = OrderedDict()
        coodict['X'] = Xc.copy()
        coodict['Y'] = Yc.copy()
        coodict['ID'] = self.Pattern['ID'].copy()

        if nested:
            nesteddict = self._coodict_2_nested(coodict)
            return nesteddict
        else:
            return coodict

    def get_CCDcoos(self, Quad, Star):
        ID = '%s%s' % (Quad, Star[0])
        ixID = np.where(self.Pattern['ID'] == ID)
        pt = self.Pattern  # alias
        Xp, Yp = pt['X'][ixID][0], pt['Y'][ixID][0]
        Xc, Yc = self.convert_Phys_2_CCD(Xp, Yp)
        return (Xc, Yc)

    def get_similaritymx(scale, rotation_deg, translation):
        """ 

        :param scale: float, scale factor.
        :param rotation_deg: float, rotation angel in ccw direction in degrees.
        :param translation: tuple/list/array, x,y translation parameters.

        """

        ss = skimage.transform.SimilarityTransform(scale=scale,
                                                   rotation=np.radians(rotation_deg), translation=translation)
        return ss.params

    def find_patt_transform(self, X, Y):
        """ """
        source = zip(self.Pattern['X'], self.Pattern['Y'])
        target = zip(X, Y)
        transf, (s_list, t_list) = aa.find_transform(source, target)
        return transf

    def _apply_transform(self, X, Y, similaritymx):

        def f(coo): return np.dot(similaritymx, coo)
        Xp, Yp, shouldbe1s = zip(*map(f, zip(X, Y, np.ones_like(X))))
        Xp = np.array(Xp)
        Yp = np.array(Yp)

        return Xp, Yp

    def apply_patt_transform(self, similaritymx):
        """ """
        X = self.Pattern['X'].copy()
        Y = self.Pattern['Y'].copy()
        Xp, Yp = self._apply_transform(X, Y, similaritymx)
        self.Pattern['X'] = Xp.copy()
        self.Pattern['Y'] = Yp.copy()


class TestStarTracker(unittest.TestCase):
    """
    Unit tests for the point-source location tools.
    """

    def setUp(self):
        """ """
        self.stracker = StarTracker(CCD=2, withpover=True)

    @unittest.skip("REDTAG")
    def test_something(self):
        """ """
        pass


def testit():
    suite = unittest.TestLoader().loadTestsFromTestCase(TestStarTracker)
    unittest.TextTestRunner(verbosity=3).run(suite)


if __name__ == '__main__':

    testit()
