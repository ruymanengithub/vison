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
from pdb import set_trace as stop

from vison.pipe import lib as plib
import copy
import numpy as np
# END IMPORT


mirror_nom = dict(Filter1=70.,Filter2=70.,Filter3=70.,
                  Filter4=70.,Filter5=70.,Filter6=70.)# nominal focus positions of mirror

Point_CooNom = {'names' : ['ALPHA','BRAVO','CHARLIE','DELTA','ECHO'],
    'CCD1':{'E':{'ALPHA':(plib.prescan+0.2*plib.imgwidth,plib.imgheight*(1.+0.8)),
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



def get_photom_spot(stamp,x0,y0):
    """ 
    
    measurements:'apflu','eapflu','bgd','ebgd'
    """
    
    results = dict()    
    return results

def get_Gauss_spot():
    """ 

    measurements:'i0_G','ei0_G','flu_G','eflu_G',
                 'x0_G','ex0_G','y0_G','ey0_G',
                 'sx_G','esx_G','sy_G','esy_G',
                 'e1','e2','e','R2'

    
    """
    
    results = dict()
    return results
    
def get_Quadrupoles_spot():
    """ 
    
    measurements:'e1','e2','e','R2'                 
    
    """
    
    results = dict()
    return results

def fwcentroid(image, checkbox=1, maxiterations=30, threshold=1e-5, halfwidth=35, 
               verbose=False, full=False):
    """ Implement the Floating-window first moment centroid algorithm
        chosen for JWST target acquisition.

        See JWST-STScI-001117 and JWST-STScI-001134 for details.

        This code makes no attempt to vectorize or optimize for speed;
        it's pretty much just a straight verbatim implementation of the
        IDL-like pseudocode provided in JWST-STScI-001117


        Parameters
        ----------
        image : array_like
            image to centroid
        checkbox : int
            size of moving checkbox for initial peak pixel guess. Default 1
        maxiterations : int
            Max number of loops. Default 30
        halfwidth : int
            Half width of the centroid box size (less 1). Specify as a scalar, or a tuple Xhalfwidth, Yhalfwidth.
            Empirical tests suggest this parameter should be at *least* the PSF FWHM for convergence,
            preferably some small factor larger
        threshold : float
            Position threshold for convergence

        Returns
        --------
        (ycen, xcen) : float tuple
            Measured centroid position. Note that this is returned in Pythonic
            Y,X order for use as array indices, etc.


        -Marshall Perrin 2011-02-11
    """
    

    if hasattr(halfwidth, '__iter__'):
        XHW, YHW = halfwidth[0:2]
    else:
        XHW, YHW = halfwidth, halfwidth

    # Determine starting peak location
    if checkbox > 1:
        raise NotImplemented("Checkbox smoothing not done yet")
    else:
        # just use brightest pixel
        w = np.where(image == image.max())
        YPEAK, XPEAK = w[0][0], w[1][0]
        if verbose: print "Peak pixels are %d, %d" % (XPEAK, YPEAK)


    # Calculate centroid for first iteration
    SUM = 0.0
    XSUM = 0.0
    XSUM2 = 0.0
    XSUM3 = 0.0
    YSUM = 0.0
    YSUM2 = 0.0
    YSUM3 = 0.0
    CONVERGENCEFLAG = False


    for i in np.arange( 2*XHW+1)+ XPEAK-XHW :
        for j in np.arange( 2*YHW+1) +YPEAK-YHW :
            XLOC = i
            YLOC = j
            SUM += image[j,i]
            XSUM += XLOC * image[j,i]
            XSUM2 += XLOC**2 * image[j,i]
            XSUM3 += XLOC**3 * image[j,i]
            YSUM += YLOC * image[j,i]
            YSUM2 += YLOC**2 * image[j,i]
            YSUM3 += YLOC**3 * image[j,i]
    XCEN = XSUM / SUM
    XMOM2 = XSUM2 / SUM
    XMOM3 = XSUM3 / SUM
    YCEN = YSUM / SUM
    YMOM2 = YSUM2 / SUM
    YMOM3 = YSUM3 / SUM

    oldXCEN = XCEN
    oldYCEN = YCEN

    if verbose:
        print( "After initial calc, cent pos is  (%f, %f)" % (XCEN, YCEN))

    # Iteratively calculate centroid until solution converges,
    # use more neighboring pixels and apply weighting:
    for k in range(maxiterations):
        SUM = 0.0
        XSUM = 0.0
        XSUM2 = 0.0
        XSUM3 = 0.0
        YSUM = 0.0
        YSUM2 = 0.0
        YSUM3 = 0.0
        for i in np.arange( 2*(XHW+1)+1)+ int(oldXCEN)-(XHW+1) :
            for j in np.arange( 2*(YHW+1)+1) +int(oldYCEN)-(YHW+1) :
                #stop()
                #-- Calculate weights
                #Initialize weights to zero:
                XWEIGHT = 0
                YWEIGHT = 0
                #Adjust weights given distance from current centroid:
                XOFF = np.abs(i - oldXCEN)
                YOFF = np.abs(j - oldYCEN)
                #If within original centroid box, set the weight to one:
                if (XOFF <= XHW): XWEIGHT = 1
                elif (XOFF > XHW) and (XOFF < XHW+1):
                    #Else if on the border, then weight needs to be scaled:
                    XWEIGHT = XHW + 1 - XOFF
                #If within original centroid box, set the weight to one:
                if (YOFF <= YHW): YWEIGHT = 1
                elif (YOFF > YHW) and (YOFF < YHW+1):
                    #Else if on the border, then weight needs to be scaled:
                    YWEIGHT = YHW + 1 - YOFF
                WEIGHT = XWEIGHT * YWEIGHT

                #Centroid, second moment, and third moment calculations
                #XLOC = i - int(XCEN) + XHW + 2
                #YLOC = j - int(YCEN) + YHW + 2
                XLOC = i
                YLOC = j

                #print "pix (%d, %d) weight %f" % (i, j, WEIGHT)
                SUM = SUM + image[j,i] * WEIGHT
                XSUM = XSUM + XLOC * image[j,i] * WEIGHT
                XSUM2 = XSUM2 + XLOC**2 * image[j,i] * WEIGHT
                XSUM3 = XSUM3 + XLOC**3 * image[j,i] * WEIGHT
                YSUM = YSUM + YLOC * image[j,i] * WEIGHT
                YSUM2 = YSUM2 + YLOC**2 * image[j,i] * WEIGHT
                YSUM3 = YSUM3 + YLOC**3 * image[j,i] * WEIGHT
        XCEN = XSUM / SUM
        XMOM2 = XSUM2 / SUM
        XMOM3 = XSUM3 / SUM
        YCEN = YSUM / SUM
        YMOM2 = YSUM2 / SUM
        YMOM3 = YSUM3 / SUM
        
        if verbose:
            print( "After iter %d , cent pos is  (%f, %f)" % (k, XCEN, YCEN))

        #Check for convergence:
        if (np.abs(XCEN - oldXCEN) <= threshold and
            np.abs(YCEN - oldYCEN) <= threshold):
            CONVERGENCEFLAG = True
            break
        else:
            CONVERGENCEFLAG = False
            oldXCEN = XCEN
            oldYCEN = YCEN

    if not CONVERGENCEFLAG:
        print("Algorithm terminated at max iterations without convergence.")
        
    if full:
        return YCEN, XCEN, YMOM2, XMOM2, YMOM3, XMOM3
    return  YCEN, XCEN
