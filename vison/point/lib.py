# -*- coding: utf-8 -*-
"""

FM-Calib. Campaign.

Library module with useful data and functions for processing of point-source 
imaging data.

Created on Wed Apr  5 10:21:05 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from vison.pipe import lib as plib
import copy
# END IMPORT


mirror_nom = dict(Filter1=70.,Filter2=70.,Filter3=70.,
                  Filter4=70.,Filter5=70.,Filter6=70.)# nominal focus positions of mirror

Point_CooNom = {'CCD1':{'E':{'ALPHA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(1.+0.8)),
                              'BRAVO':(plib.prescan+0.8*plib.imgwidth,plib.imgheight*(1.+0.8)),
                             'CHARLIE':(plib.prescan+0.5*plib.imgwidth,plib.imgheight*(1.+0.5)),
                              'DELTA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(1.+0.2)),
                              'ECHO':(plib.prescan+0.8*plib.imgwidth,plib.imgheight*(1.+0.2))},
                        'F':{'ALPHA':(plib.quad_width+plib.prescan+0.2*plib.imgwidth,plib.imgheight*(1.+0.8)),
                            'BRAVO':(plib.quad_width+plib.prescan+0.8*plib.imgwidth,plib.imgheight*(1.+0.8)),
                            'CHARLIE':(plib.quad_width+plib.prescan+0.5*plib.imgwidth,plib.imgheight*(1.+0.5)),
                            'DELTA':(plib.quad_width+plib.prescan+0.2*plib.imgwidth,plib.imgheight*(1.+0.2)),
                            'ECHO':(plib.quad_width+plib.prescan+0.8*plib.imgwidth,plib.imgheight*(1.+0.2))},
                        'H':{'ALPHA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(0.8)),
                             'BRAVO':(plib.prescan+0.8*plib.imgwidth,plib.imgheight*(0.8)),
                            'CHARLIE':(plib.prescan+0.5*plib.imgwidth,plib.imgheight*(0.5)),
                            'DELTA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(0.2)),
                            'ECHO':(plib.prescan+0.8*plib.imgwidth,plib.imgheight*(0.2))},
                        'G':{'ALPHA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(0.8)),
                            'BRAVO':(plib.prescan+0.8*plib.imgwidth,plib.imgheight*(0.8)),
                            'CHARLIE':(plib.prescan+0.5*plib.imgwidth,plib.imgheight*(0.5)),
                            'DELTA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(0.2)),
                            'ECHO':(plib.prescan+0.8*plib.imgwidth,plib.imgheight*(0.2))}
                            }}

for iCCD in range(2,4):
    Point_CooNom['CCD%i' % iCCD] = dict()
    for Q in plib.Quads:
        Point_CooNom['CCD%i' % iCCD][Q] = copy.deepcopy(Point_CooNom['CCD1'][Q])



