# -*- coding: utf-8 -*-
"""


Library module with useful data and functions for processing of point-source 
imaging data.

Created on Wed Apr  5 10:21:05 2017

:author: Ruyman Azzollini (except where indicated)

"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict

from vison.datamodel import ccd as ccdmod
#from vison.pipe import lib as plib
from vison.support import context
import copy
import numpy as np
from vison.point import spot as spotmod
from vison.datamodel import ccd
from vison.image import sextractor as sex
# END IMPORT


stampw = 75

# ALPHA BRAVO
#    CHARLIE
# DELTA ECHO
#
# Positions are CCD-based (not Q-based), and referenced to output node of H Quad (lower left)

Point_CooNom = {'names': ['ALPHA', 'BRAVO', 'CHARLIE', 'DELTA', 'ECHO'],
                'CCD1': {'E': {'ALPHA': (context.prescan+0.2*context.imgwidth, context.imgheight*(1.+0.8)),
                               'BRAVO': (context.prescan+0.8*context.imgwidth, context.imgheight*(1.+0.8)),
                               'CHARLIE': (context.prescan+0.5*context.imgwidth, context.imgheight*(1.+0.5)),
                               'DELTA': (context.prescan+0.2*context.imgwidth, context.imgheight*(1.+0.2)),
                               'ECHO': (context.prescan+0.8*context.imgwidth, context.imgheight*(1.+0.2))},
                         'F': {'ALPHA': (context.quad_width+context.prescan+0.2*context.imgwidth, context.imgheight*(1.+0.8)),
                               'BRAVO': (context.quad_width+context.prescan+0.8*context.imgwidth, context.imgheight*(1.+0.8)),
                               'CHARLIE': (context.quad_width+context.prescan+0.5*context.imgwidth, context.imgheight*(1.+0.5)),
                               'DELTA': (context.quad_width+context.prescan+0.2*context.imgwidth, context.imgheight*(1.+0.2)),
                               'ECHO': (context.quad_width+context.prescan+0.8*context.imgwidth, context.imgheight*(1.+0.2))},
                         'H': {'ALPHA': (context.prescan+0.2*context.imgwidth, context.imgheight*(0.8)),
                               'BRAVO': (context.prescan+0.8*context.imgwidth, context.imgheight*(0.8)),
                               'CHARLIE': (context.prescan+0.5*context.imgwidth, context.imgheight*(0.5)),
                               'DELTA': (context.prescan+0.2*context.imgwidth, context.imgheight*(0.2)),
                               'ECHO': (context.prescan+0.8*context.imgwidth, context.imgheight*(0.2))},
                         'G': {'ALPHA': (context.quad_width+context.prescan+0.2*context.imgwidth, context.imgheight*(0.8)),
                               'BRAVO': (context.quad_width+context.prescan+0.8*context.imgwidth, context.imgheight*(0.8)),
                               'CHARLIE': (context.quad_width+context.prescan+0.5*context.imgwidth, context.imgheight*(0.5)),
                               'DELTA': (context.quad_width+context.prescan+0.2*context.imgwidth, context.imgheight*(0.2)),
                               'ECHO': (context.quad_width+context.prescan+0.8*context.imgwidth, context.imgheight*(0.2))}
                         }}

for iCCD in range(2, 4):
    Point_CooNom['CCD%i' % iCCD] = dict()
    for Q in context.Quads:
        Point_CooNom['CCD%i' % iCCD][Q] = copy.deepcopy(
            Point_CooNom['CCD1'][Q])

# MEASUREMENTS

Point_CooNom['CCD2'] = OrderedDict(E=OrderedDict(ALPHA=(542.0, 1725.0),
                                                 BRAVO=(1716.0, 1700.0),
                                                 CHARLIE=(1126.0, 1124.0),
                                                 DELTA=(521.0, 551.0),
                                                 ECHO=(1695.0, 537.0)),
                                   F=OrderedDict(ALPHA=(578.0, 1686.0),
                                                 BRAVO=(1745.0, 1668.0),
                                                 CHARLIE=(1141.0, 1144.0),
                                                 DELTA=(553.0, 522.0),
                                                 ECHO=(1723.0, 496.0)),
                                   G=OrderedDict(ALPHA=(534.0, 1590.0),
                                                 BRAVO=(1702.0, 1571.0),
                                                 CHARLIE=(1139.0, 1033.0),
                                                 DELTA=(515.0, 415.0),
                                                 ECHO=(1685.0, 394.0)),
                                   H=OrderedDict(ALPHA=(554.0, 1626.0),
                                                 BRAVO=(1725.0, 1606.0),
                                                 CHARLIE=(1131.0, 1029.0),
                                                 DELTA=(531.0, 446.0),
                                                 ECHO=(1706.0, 435.0)))


def extract_spot(ccdobj, coo, Quad, log=None, stampw=25):
    """ """

    B = ccdobj.QuadBound[Quad]

    x0Q, y0Q = tuple([int(np.round(item)) for item in coo])
    x0Q -= B[0]
    y0Q -= B[2]

    corners = [x0Q-stampw/2, x0Q-stampw/2+stampw,
               y0Q-stampw/2, y0Q-stampw/2+stampw]

    # Cut-out stamp of the spot
    stamp = ccdobj.get_cutout(corners, Quad, canonical=False)
    
    lowerleft = (corners[0]+B[0], corners[2]+B[2])
    spot = spotmod.Spot(stamp, log=log, lowerleft=lowerleft)

    return spot


def gen_point_mask(Quad, width=stampw, sources='all', coodict=Point_CooNom):
    """ """

    fkccdobj = ccdmod.CCD()

    NAXIS1, NAXIS2 = fkccdobj.NAXIS1, fkccdobj.NAXIS2

    msk = np.zeros((NAXIS1, NAXIS2), dtype='int8')

    if sources == 'all':
        sourcelist = coodict['names']
    else:
        sourcelist = copy.deepcopy(sources)

    for source in sourcelist:

        coo = coodict[Quad][source]

        x0Q, y0Q = tuple([int(np.round(item)) for item in coo])

        corners = [x0Q-stampw/2, x0Q-stampw/2+stampw,
                   y0Q-stampw/2, y0Q-stampw/2+stampw]

        msk[corners[0]:corners[1], corners[2]:corners[3]] = 1

    return msk


class _Transf(object):
    rotation = 1.4*np.pi
    scale = 2.0
    translation = [2000.0,2000.0]


def _get_transformation(self,ffits, stracker, 
                dpath, _sexconfig, tag, devel=False):   
        
            
    #vstart = self.dd.mx['vstart'][iObs][jCCD]
    #vend = self.dd.mx['vend'][iObs][jCCD]
    
    SExCatroot = 'StarFinder_%s' % tag
    
    # DEFAULTS
    transf = _Transf()
    s_list = np.arange(5).tolist()
    t_list = copy.deepcopy(s_list)
    
    if not devel:
        
        ccdobj = ccd.CCD(ffits)
        
        img = ccdobj.extensions[-1].data.T.copy()
        
        SExCat = sex.easy_run_SEx(img, SExCatroot, 
                    sexconfig=_sexconfig,                                      
                    cleanafter=True)
        
        sel = np.where((SExCat['ELONGATION']<2.) &
                       (SExCat['A_IMAGE']<5.) &
                       (SExCat['A_IMAGE']>0.2) &
                       (SExCat['B_IMAGE']<5.) &
                       (SExCat['B_IMAGE']>0.2) &
                       (SExCat['FLUX_AUTO']>5000.) &
                       (SExCat['ISOAREA_IMAGE']>3.) &
                       (SExCat['ISOAREA_IMAGE']<100.)
                       )
        
        if len(sel[0])> 0:
            
            X_IMAGE = SExCat['X_IMAGE'][sel].copy() - 1.
            Y_IMAGE = SExCat['Y_IMAGE'][sel].copy() - 1.
            
            
            X_PHYS, Y_PHYS = self.ccdcalc.cooconv_CCD_2_Phys(X_IMAGE,Y_IMAGE)                    
            whatQ = self.ccdcalc.get_Q_of_CCDcoo(X_IMAGE,Y_IMAGE)[0]
            
            ixnonan = np.where(~(np.isnan(X_PHYS) | np.isnan(Y_PHYS)) &
                    (whatQ != 'G'))
            X_PHYS = X_PHYS[ixnonan]
            Y_PHYS = Y_PHYS[ixnonan]
            
            
            try:
                transf, (s_list, t_list) = stracker.find_patt_transform(X_PHYS,Y_PHYS,
                          Full=True, discardQ=['G'],debug=False)
            except:
            
                if self.log is not None:
                    self.log.info('Failed Locking Stars on %s' % tag)
        else:
            if self.log is not None:
                    self.log.info('Did not Find enough Stars to Lock-on on %s' % tag)
        
        return transf, s_list, t_list