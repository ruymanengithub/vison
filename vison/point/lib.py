# -*- coding: utf-8 -*-
"""

FM-Calib. Campaign.

Library module with useful data and functions for processing of point-source 
imaging data.

Created on Wed Apr  5 10:21:05 2017

:author: Ruyman Azzollini (except where indicated)
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop

from vison.datamodel import ccd as ccdmod
#from vison.pipe import lib as plib
from vison.support import context
import copy
import numpy as np
from vison.point import spot as spotmod
# END IMPORT



mirror_nom = dict(F1=70.,F2=70.,F3=70.,F4=70.,F5=70.,F6=70.) # nominal focus positions of mirror

stampw = 25

# ALPHA BRAVO
#    CHARLIE
# DELTA ECHO
#
# Positions are CCD-based (not Q-based), and referenced to output node of H Quad (lower left)

Point_CooNom = {'names' : ['ALPHA','BRAVO','CHARLIE','DELTA','ECHO'],
    'CCD1':{'E':{'ALPHA':(context.prescan+0.2*context.imgwidth,context.imgheight*(1.+0.8)),
                              'BRAVO':(context.prescan+0.8*context.imgwidth,context.imgheight*(1.+0.8)),
                             'CHARLIE':(context.prescan+0.5*context.imgwidth,context.imgheight*(1.+0.5)),
                              'DELTA':(context.prescan+0.2*context.imgwidth,context.imgheight*(1.+0.2)),
                              'ECHO':(context.prescan+0.8*context.imgwidth,context.imgheight*(1.+0.2))},
            'F':{'ALPHA':(context.quad_width+context.prescan+0.2*context.imgwidth,context.imgheight*(1.+0.8)),
                            'BRAVO':(context.quad_width+context.prescan+0.8*context.imgwidth,context.imgheight*(1.+0.8)),
                            'CHARLIE':(context.quad_width+context.prescan+0.5*context.imgwidth,context.imgheight*(1.+0.5)),
                            'DELTA':(context.quad_width+context.prescan+0.2*context.imgwidth,context.imgheight*(1.+0.2)),
                            'ECHO':(context.quad_width+context.prescan+0.8*context.imgwidth,context.imgheight*(1.+0.2))},
            'H':{'ALPHA':(context.prescan+0.2*context.imgwidth,context.imgheight*(0.8)),
                             'BRAVO':(context.prescan+0.8*context.imgwidth,context.imgheight*(0.8)),
                            'CHARLIE':(context.prescan+0.5*context.imgwidth,context.imgheight*(0.5)),
                            'DELTA':(context.prescan+0.2*context.imgwidth,context.imgheight*(0.2)),
                            'ECHO':(context.prescan+0.8*context.imgwidth,context.imgheight*(0.2))},
            'G':{'ALPHA':(context.quad_width+context.prescan+0.2*context.imgwidth,context.imgheight*(0.8)),
                            'BRAVO':(context.quad_width+context.prescan+0.8*context.imgwidth,context.imgheight*(0.8)),
                            'CHARLIE':(context.quad_width+context.prescan+0.5*context.imgwidth,context.imgheight*(0.5)),
                            'DELTA':(context.quad_width+context.prescan+0.2*context.imgwidth,context.imgheight*(0.2)),
                            'ECHO':(context.quad_width+context.prescan+0.8*context.imgwidth,context.imgheight*(0.2))}
                            }}

for iCCD in range(2,4):
    Point_CooNom['CCD%i' % iCCD] = dict()
    for Q in context.Quads:
        Point_CooNom['CCD%i' % iCCD][Q] = copy.deepcopy(Point_CooNom['CCD1'][Q])


def extract_spot(ccdobj,coo,Quad,log=None,stampw=25):
    """ """
    
    B = ccdobj.QuadBound[Quad]
    
    x0Q,y0Q = tuple([int(np.round(item)) for item in coo])
    x0Q -= B[0]
    y0Q -= B[2]
    
    corners = [x0Q-stampw/2,x0Q-stampw/2+stampw,
                y0Q-stampw/2,y0Q-stampw/2+stampw]
    
    # Cut-out stamp of the spot
    stamp = ccdobj.get_cutout(corners,Quad,canonical=False)
    
    spot = spotmod.Spot(stamp,log=log)
    
    return spot

def gen_point_mask(CCD,Quad,width=stampw,sources='all'):
    """ """
    
    fkccdobj = ccdmod.CCD()
    
    NAXIS1,NAXIS2 = fkccdobj.NAXIS1,fkccdobj.NAXIS2
    
    msk = np.zeros((NAXIS1,NAXIS2),dtype='int8')
    
    allsources = Point_CooNom['names']
    if sources == 'all':
        sourcelist = allsources
    else:
        sourcelist = [sources]
    
    for source in sourcelist:
        
        coo = Point_CooNom['CCD%i' % CCD][Quad][source]
        
        x0Q,y0Q = tuple([int(np.round(item)) for item in coo])
    
        corners = [x0Q-stampw/2,x0Q-stampw/2+stampw,
                y0Q-stampw/2,y0Q-stampw/2+stampw]
        
        msk[corners[0]:corners[1],corners[2]:corners[3]] = 1
    
    
    return msk

    


