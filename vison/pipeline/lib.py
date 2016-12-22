# -*- coding: utf-8 -*-
"""

Support Script with Variables and Functions used for FM-Calib. analysis

Created on Wed Nov 30 11:11:27 2016


:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk
"""

# IMPORT STUFF
import copy
from vissim.datamodel.ccd import CCD
# END IMPORT



Quads =['E','F','G','H']
RON = 4.5 # e-
gain = 3.5 # e-/ADU

ccdobj = CCD()

prescan = ccdobj.prescan
overscan = ccdobj.overscan
imgheight = ccdobj.NAXIS2/2
quad_width = ccdobj.NAXIS1/2
imgwidth = quad_width - prescan - overscan

FW = dict(Filter1='570nm',Filter2='650nm',
          Filter3='700nm',Filter4='800nm',
          Filter5='900nm',Filter6='ND')

Point_CooNom = {'CCD1':{'E':{'ALPHA':(prescan+0.2*imgwidth,imgheight*(1.+0.8)),
                              'BRAVO':(prescan+0.8*imgwidth,imgheight*(1.+0.8)),
                             'CHARLIE':(prescan+0.5*imgwidth,imgheight*(1.+0.5)),
                              'DELTA':(prescan+0.2*imgwidth,imgheight*(1.+0.2)),
                              'ECHO':(prescan+0.8*imgwidth,imgheight*(1.+0.2))},
                        'F':{'ALPHA':(quad_width+prescan+0.2*imgwidth,imgheight*(1.+0.8)),
                            'BRAVO':(quad_width+prescan+0.8*imgwidth,imgheight*(1.+0.8)),
                            'CHARLIE':(quad_width+prescan+0.5*imgwidth,imgheight*(1.+0.5)),
                            'DELTA':(quad_width+prescan+0.2*imgwidth,imgheight*(1.+0.2)),
                            'ECHO':(quad_width+prescan+0.8*imgwidth,imgheight*(1.+0.2))},
                        'H':{'ALPHA':(prescan+0.2*imgwidth,imgheight*(0.8)),
                             'BRAVO':(prescan+0.8*imgwidth,imgheight*(0.8)),
                            'CHARLIE':(prescan+0.5*imgwidth,imgheight*(0.5)),
                            'DELTA':(prescan+0.2*imgwidth,imgheight*(0.2)),
                            'ECHO':(prescan+0.8*imgwidth,imgheight*(0.2))},
                        'G':{'ALPHA':(prescan+0.2*imgwidth,imgheight*(0.8)),
                            'BRAVO':(prescan+0.8*imgwidth,imgheight*(0.8)),
                            'CHARLIE':(prescan+0.5*imgwidth,imgheight*(0.5)),
                            'DELTA':(prescan+0.2*imgwidth,imgheight*(0.2)),
                            'ECHO':(prescan+0.8*imgwidth,imgheight*(0.2))}
                            }}

for iCCD in range(2,4):
    Point_CooNom['CCD%i' % iCCD] = dict()
    for Q in Quads:
        Point_CooNom['CCD%i' % iCCD][Q] = copy.deepcopy(Point_CooNom['CCD1'][Q])
