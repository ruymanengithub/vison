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

default_patt_files = OrderedDict(
    CCD1=os.path.join(ogsepath, 'Pattern_CCD1_CHAMBER_A_JUN18.txt'),
    CCD2=os.path.join(ogsepath, 'Pattern_CCD2_CHAMBER_A_JUN18.txt'),
    CCD3=os.path.join(ogsepath, 'Pattern_CCD3_CHAMBER_A_JUN18.txt'),)


starnames = ['ALPHA', 'BRAVO', 'CHARLIE', 'DELTA', 'ECHO']


def sort_coordinate_pairs(coopairs):
    """ """
    return sorted(coopairs, key=lambda k: [k[0], k[1]])


class StarTracker(object):

    def __init__(self, CCD, withpover=True, patt_files=None):

        self.Quads = ['E', 'F', 'G', 'H']
        self.starnames = starnames
        if patt_files is None:
            self.patt_files = default_patt_files.copy()
        else:
            self.patt_files = patt_files.copy()

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

    def get_similaritymx(self, scale, rotation_rad, translation):
        """

        :param scale: float, scale factor.
        :param rotation_deg: float, rotation angel in ccw direction in degrees.
        :param translation: tuple/list/array, x,y translation parameters.

        """

        # ss = skimage.transform.SimilarityTransform(scale=scale,
        #            rotation=np.radians(rotation_deg), translation=translation)
        ss = skimage.transform.SimilarityTransform(scale=scale,
                                                   rotation=rotation_rad, translation=translation)
        return ss.params

    def find_patt_transform(self, Xt, Yt, Full=False, discardQ=None, debug=False):
        """ """
        Xs = self.Pattern['X']
        Ys = self.Pattern['Y']
        if discardQ is not None:
            ID = self.Pattern['ID']
            sel = np.array([i for i in range(len(ID)) if ID[i][0] not in discardQ])
            Xs = Xs[sel].copy()
            Ys = Ys[sel].copy()

        source = list(zip(Xs, Ys))
        target = list(zip(Xt, Yt))

        source = sort_coordinate_pairs(source)
        target = sort_coordinate_pairs(target)

        xs, ys = tuple(zip(*source))
        xt, yt = tuple(zip(*target))

        try:
            transf, (s_list, t_list) = aa.find_transform(source, target)

        except BaseException:
            raise RuntimeError

        if debug:
            #from pylab import plot,show
            from matplotlib import pyplot as plt
            mx = self.get_similaritymx(transf.scale, transf.rotation, transf.translation)
            xtp, ytp = self._apply_transform(s_list[:, 0], s_list[:, 1], mx)
            fig = plt.figure()
            ax1 = fig.add_subplot(121)
            ax1.plot(xtp, ytp, 'ro')
            ax1.plot(xt, yt, 'k.')
            ax2 = fig.add_subplot(122)
            ax2.plot(xtp - t_list[:, 0], ytp - t_list[:, 1], 'bo')
            plt.show()

        if Full:
            return transf, (s_list, t_list)
        else:
            transf

    def _apply_transform(self, X, Y, similaritymx):

        def f(coo): return np.dot(similaritymx, coo)
        Xp, Yp, shouldbe1s = list(zip(*list(map(f, list(zip(X, Y, np.ones_like(X)))))))
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
