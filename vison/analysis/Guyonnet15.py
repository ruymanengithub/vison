# -*- coding: utf-8 -*-
"""

Library with functions that implement the algorithms described in Guyonnet+15.
"Evidence for self-interaction of charge distribution in CCDs"
Guyonnet, Astier, Antilogus, Regnault and Doherty 2015


Notes:
    
- I renamed "x" (pixel boundary index) to "b", to avoid confusion with
  cartesian "x".
- In paper, X belongsto [(0,1),(1,0),(0,-1),(-1,0)]. Here
  b is referred to as cardinal points "N","E","S","W". 
  It is linked to matrix index ib, running between 0 and 3.
     
Created on Thu Sep 22 11:38:24 2016

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop

import numpy as np
from scipy.special import expi as Ei
from scipy import optimize as optimize
import copy
from scipy.signal import convolve2d

from astropy.io import fits as fts
from matplotlib import pyplot as plt
# END IMPORT

pixpitch = 12. # um, pixel pitch

def frdist(i,j,ib):
    """ 
    Distance from the source charge to considered boundary "b"
    
    :param i: pixel coordinate i
    :param j: pixel coordinate j
    :param ib: boundary index [0, 1, 2, 3]
    
    :return: distance r(ijb)
    
    
    """
    x = i*1. + np.cos((1.-ib)*90*np.pi/180.)*0.5
    y = j*1. + np.sin((1.-ib)*90*np.pi/180.)*0.5
    
    return (x*x+y*y)**0.5

def ftheta_bound(i,j,ib):
    """
    "[theta_i,j^X is] the angle between the source-boundary vector and the
    normal to the boundary".
           
    :param i: pixel coordinate i
    :param j: pixel coordinate j
    :param ib: boundary index [0, 1, 2, 3]

    :return: theta_i,j^x
    
    """
    
    # source-boundary vector
    
    sb_v = np.array([i*1.+np.cos((1.-ib)*90*np.pi/180.)*0.5,
                     j*1.+np.sin((1.-ib)*90.*np.pi/180.)*0.5])
    
    # normal vector # unit vector pointing 'outwards' of pixel i,j
    n_v = np.array([np.cos((1.-ib)*90*np.pi/180.),
                    np.sin((1.-ib)*90.*np.pi/180.)])
    
    
    arg = np.dot(sb_v,n_v)
    modulus = np.sqrt(np.dot(sb_v,sb_v)) * 1.
    
    return np.arccos(arg/modulus)
    

def fpred_aijb(p,i,j,ib):
    """
    
    'The smoothing model assumes that a_{ij}^x coefficients are the
    product of a function of distance from the source charge to the
    considered boundary (r_{ij}) and that it also trivially depends on
    the angle between the source-boundary vector and the normal to
    the boundary (theta_{i,j}^x)'
    
    Eq. 18
    
    :param p: parameters of the radial function (list of 2)
    :param i: pixel coordinate i
    :param j: pixel coordinate j
    :param ib: boundary index [0, 1, 2, 3]
    
    :return: f(rij)cos(theta_ij^x)
    
    """
    
    p0,p1 = p
    
    rdist = frdist(i,j,ib)
    
    fr = p0 * (Ei(p1*rdist))
    
    theta_ijb = ftheta_bound(i,j,ib)
    
    aijb = fr * np.cos(theta_ijb)
    
    #if np.abs(aijb) > 3.: stop()

    return aijb


def fun_p(x,*p):
    """auxiliary function to 'solve_for_psmooth'
    """
    
    ii = x[0,:]
    jj = x[1,:]
    
    model = np.zeros(len(ii),dtype='float32')
    
    for ix in range(len(ii)):
        i = ii[ix]
        j = jj[ix]
        for ib in range(4):    
            model[ix] += fpred_aijb(p,i,j,ib)
            
    return model



def solve_for_psmooth(covij,var,mu,doplot=False):
    """Solving (p0,p1) parameters in Eq. 18 using covariance matrix and 
    measured covariance  matrix.
    
    :param covij: array, covariance matrix
    :param var: float, variance
    :param mu: float, expected value of pixel values ("mean" of flat-field)
    :param doplot: bool, if True, plot data and best fit model 
    
    :return: best-fit parameters, and errors: 2 tuples of 2 elements each
    
    """
    
    ni,nj = covij.shape # shape of covariance matrix
    
    assert ni == nj # only squared covariance matrices allowed in this  bar
    N = ni
    
    i = np.arange(N)
    j = np.arange(N)
    
    ii,jj = np.meshgrid(i,j,indexing='xy') # 2D index arrays
    
    ixsel = np.where((ii.flatten() >= 1) & (jj.flatten()>=1))
         # 'No constraints are applied on the displacement of the boundaries
         # of the pixels (1,0) and (0,1) because the analysis of correlation
         # properties presented in Sect. 4.3.1 indicates a different
         # behaviour for those than for the rest of the coefficients.
         
    nsel = len(ixsel[0])
    
    xdata = np.zeros((2,nsel),dtype='float32')
    xdata[0,:] = ii.flatten()[ixsel]
    xdata[1,:] = jj.flatten()[ixsel]
    
    r = (xdata[0,:]**2.+xdata[1,:]**2.)**0.5 * pixpitch # um
    ydata = (covij/(var*mu)).flatten()[ixsel]
    
    order = r.argsort()
    r = r[order]
    ydata = ydata[order].copy()
    xdata[0,:] = xdata[0,order].copy()
    xdata[1,:] = xdata[1,order].copy()
    
    p0 = (1.E-3,-1.) # initial values
         
    
    popt,pcov = optimize.curve_fit(fun_p,xdata,ydata,p0=p0,sigma=None,absolute_sigma=False)
    epopt = np.sqrt(np.diagonal(pcov)) # errors
    
    if doplot:
        print 'popt = #', popt
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(r,ydata,'ko:',label='data')
        ax.plot(r,fun_p(xdata,*p0),'r--',label='initial guess')
        ax.plot(r,fun_p(xdata,*popt),'b--',label='solution')
        ax.set_xlabel('Distance (um)')
        ax.set_ylabel('covij/(var*mu)')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles,labels,loc='upper right')
        plt.tight_layout()
        plt.show()
    
    return popt,epopt

def _build_aijb(X,indepBounds,psmooth):
    """ 
    
    Builds the matrix aijb, out of a matrix of independent coefficients "X",
    the (ijb) indices of such coefficients, and the parameters of the prediction
    function fpred(ijb).
    
    indexes of aijb
    
    i : 0, n-1
    j : 0, n-1
    x : 0, 3 (N,E,S,W) : (0, 90, 180, 270) degrees
    
    :param X: array, matrix N^2-1 with independent coefficients
    :param indepBounds: list of 3-tuples, indices of the independent boundary coefficients aijb
    :param psmooth: p parameters (2) of fpred(p,ijb)
    :return: matrix aijx
    
    """
    N = int((len(X)+1)**0.5)
    aijb = np.zeros((N,N,4),dtype='float32') + np.nan
    
    # First we asign the independent coefficients
    
    for zx, indepBound in enumerate(indepBounds):
        i,j,x = indepBound
        aijb[i,j,x] = X[zx]
    
    # "Uncompression" of the coefficients
    # It must be done in reverse order to the "compression"
    
    # First elements not on First Row or First Column
    
    # North of lower is like South of upper, if not at top floor
    
    for ip in range(N-1,1-1,-1):
        for jp in range(N-1,1-1,-1):
              
            #print ip, jp
                 
            if jp < N-1:
                # North of lower is like South of upper, if not at top (inverted)
                aijb[ip,jp,0] = aijb[ip,jp+1,2]*-1.
            elif jp == N-1:
                # if at top, use prediction
                aijb[ip,jp,0] = fpred_aijb(psmooth[0],ip,jp,0)

    # West is like South times factor
    
    for ip in range(N-1,1-1,-1):
        for jp in range(N-1,1-1,-1):
            
            #print ip, jp
            
            ratioWS = fpred_aijb(psmooth[0],ip,jp,3.)/fpred_aijb(psmooth[0],ip,jp,2.)#*\
                               #(ftheta_bound(ip,jp,2)/ftheta_bound(ip,jp,2))**2. # right?
                               
            aijb[ip,jp,3] = aijb[ip,jp,2] * ratioWS

    
    # East is like west of right
    
    for ip in range(N-1,1-1,-1):
        for jp in range(N-1,1-1,-1):
            
            #print ip,jp
            
            if ip<N-1:
                aijb[ip,jp,1] = aijb[ip+1,jp,3]*-1. # (inverted)
            elif ip == N-1:
                aijb[ip,jp,1] = fpred_aijb(psmooth[0],ip,jp,1.) # prediction
    
    # First Row: jp==0
    
    jp=0
    for ip in range(N-1,1-1,-1):
        
        #print ip,jp
        
        if ip < N-1:
            # East is like West of Right
            aijb[ip,jp,1] = aijb[ip+1,jp,3]*-1. # inverted
        elif ip == N-1:
            aijb[ip,jp,1] = fpred_aijb(psmooth[0],ip,jp,1.) # predict
        
        # North is like south of upper (inverted)
        aijb[ip,jp,0] = aijb[ip,jp+1,2] * -1.
        
        # South is like North on the horizontal
        aijb[ip,jp,2] = aijb[ip,jp,0]
        
        
    # First Column: ip == 0
    
    ip = 0
    
    for jp in range(N-1,-1,-1):
        
        #print ip,jp
        
        if (jp>=0) and (jp < N-1):
            # North is like South of Upper if not at top (inverted)
            aijb[ip,jp,0] = aijb[ip,jp+1,2]*-1.
        elif jp == N-1:
            aijb[ip,jp,0] = fpred_aijb(psmooth[0],ip,jp,0.)
        
        
        if jp ==0:
            # South is like North on the (0,0) pixel
            aijb[ip,jp,2] = aijb[ip,jp,0]
        
        #East is like West of Right (inverted)
        aijb[ip,jp,1] = aijb[ip+1,jp,3]*-1.
        
        # West is like East on the vertical
        aijb[ip,jp,3] = aijb[ip,jp,1]
        
    
    return aijb


def solve_for_A_linalg(covij,var=1.,mu=1.,doplot=False,psmooth=None,returnAll=False):
    """Function to retrieve the A matrix of pixel boundaries
    displacements, given a matrix of pixel covariances, variance, and mu.
    
    if var==1 and mu==1, it is understood that covij is
    the correlation matrix.
    
    See section 6.1 of G15.
    
    :param covij: array, squared matrix with pixel covariances.
    :param var: float, variance of the flat-field.
    :param mu: float, mean value of the flat-field.
    :param doplot: if True, plot the fit of the fpred(ijb) function
    :param psmooth: coefficients of the fpred(aijb) function (Eq. 18)
    :param returnAll: bool, controls return values
    
    :return: if returnAll == True, return (aijb, psmooth), otherwise return aijb only
    
    """
    assert isinstance(covij,np.ndarray)
    assert covij.shape[0] == covij.shape[1] # it is a squared matrix
    
    N = int((covij.size)**0.5) # side of cov matrix
    
    if psmooth is None:
        psmooth = solve_for_psmooth(covij,var,mu,doplot=doplot) # smoothing parameters, fpred(ijb)
    
    Rij = covij / (var*mu) # correlation matrix
    
    B = np.zeros(N**2-1,dtype='float32') # Independent terms
    Bini = B.copy()
    
    for i in range(N):
        for j in range(N):
            if i == 0 and j==0: continue

            ix = i*N+j - 1
            
            B[ix] = Rij[i,j]
            Bini[ix] = Rij[i,j]
    
    Afull = np.zeros((4*N**2,N**2-1),dtype='float32') # N-bounds x N-corr
    
    # Boundaries Names List
    
    boundsList = []
    for i in range(N):
        for j in range(N):
            for x in range(4):
                boundsList.append((i,j,x))
    
    # Full Coefficients Matrix, Without Reduction
    
    for i in range(N):
        for j in range(N):
            
            if i == 0 and j==0: continue
            
            jyy = i * N + j - 1 # correlation index
            
            for ip in range(N):
                for jp in range(N):
                    
                    for x in range(4):
                        
                        ixx = (ip * N + jp)*4 + x # boundary index
                        
                        if (i==ip) and (j==jp): Afull[ixx,jyy] = 1.
    
    print 'sum(A)=%.1f' % Afull.sum()

    #  REDUCING THE COEFFICIENTS MATRIX
    
    for i in range(N):
        for j in range(N):
            
            if i ==0 and j==0: continue
            
            jyy = i * N + j -1 # correlation index
            
            # 1st Column: ip == 0
            
            ip = 0
            for jp in range(N):
                
                ixx = (ip * N + jp)*4
                
                ixxN = ixx + 0 #
                ixxE = ixx + 1 #
                ixxS = ixx + 2 #
                ixxW = ixx + 3 #
                
                ixpN=(ip * N + jp+1)*4 #
                ixpE=((ip+1) * N + jp)*4 #
                #print ip, jp
                
                # East is like West on the vertical 

                Afull[ixxE,jyy] += Afull[ixxW,jyy]
                Afull[ixxW,jyy] = 0.
                
                # West of Right is like East 
                Afull[ixpE+3,jyy] = (Afull[ixpE+3,jyy]+Afull[ixxE,jyy])*1.
                Afull[ixxE,jyy] = 0.
                
                if jp ==0: #(0,0)
                    # North is like South on the (0,0) pixel
                    Afull[ixxN,jyy] += Afull[ixxS,jyy]
                    Afull[ixxS,jyy] = 0. 
                
                if (jp>=0) and (jp < N-1):
                    # South of Upper is like North if not at top 
                    Afull[ixpN+2,jyy] = (Afull[ixpN+2,jyy] + Afull[ixxN,jyy])*1.
                    Afull[ixxN,jyy] = 0.
                elif jp == N-1:
                    # if at top, add constant to B
                    Afull[ixxN,jyy] = 0.
                    if (i==ip) and (j==jp):
                        B[jyy] -= fpred_aijb(psmooth[0],ip,jp,0)                                 

            # First Row: jp==0 (excluding (0,0))
            
            jp = 0
            for ip in range(1,N):
                
                ixx = (ip * N + jp)*4
                
                ixxN = ixx + 0 #
                ixxE = ixx + 1 #             
                ixxS = ixx + 2 #
                
                ixpN=(ip * N + jp+1)*4 #
                ixpE=((ip+1) * N + jp)*4 #
                
                
                # North is like South on the horizontal
                Afull[ixxN,jyy] += Afull[ixxS,jyy]
                Afull[ixxS,jyy] = 0.
                    
                # South of Upper is like North if not at top 
                Afull[ixpN+2,jyy] = ( Afull[ixpN+2,jyy] + Afull[ixxN,jyy])*1.
                Afull[ixxN,jyy] = 0.
                
                
                if (ip < N-1):
                    # West of Right is like East
                    Afull[ixpE+3,jyy] = ( Afull[ixpE+3,jyy] + Afull[ixxE,jyy] )*1.
                    Afull[ixxE,jyy] = 0.
                elif ip == N-1:
                    Afull[ixxE,jyy] = 0.
                    if (i==ip) and (j==jp):
                        B[jyy] -= fpred_aijb(psmooth[0],ip,jp,1.)

            # The rest
            
            # West of Right is like East

            for ip in range(1,N):
                for jp in range(1,N):
                    
                    ixx = (ip * N + jp)*4 #
                    ixxE = ixx + 1 #
                    ixpE= ((ip+1) * N + jp)*4 #
            
                    if ip < N-1:
                        # West of Right is like East 
                        Afull[ixpE+3,jyy] = ( Afull[ixpE+3,jyy] + Afull[ixxE,jyy]) * 1.
                        Afull[ixxE,jyy] = 0.
                    elif ip == N-1:
                        Afull[ixxE,jyy] = 0.
                        if (i==ip) and (j==jp):
                            B[jyy] -= fpred_aijb(psmooth[0],ip,jp,1.) 

            # South is like West times factor
            
            for ip in range(1,N):
                for jp in range(1,N):
                    
                    ixx = (ip * N + jp)*4
                    ixxS = ixx + 2 #
                    ixxW = ixx + 3 #
                        
                    # South is like West times factor 
                    ratioWS = fpred_aijb(psmooth[0],ip,jp,3.)/fpred_aijb(psmooth[0],ip,jp,2.)#*\
                               #(ftheta_bound(ip,jp,2)/ftheta_bound(ip,jp,2))**2.

                    Afull[ixxS,jyy] += (Afull[ixxW,jyy]) * ratioWS  
                    Afull[ixxW,jyy] = 0 
                    
            # South of upper is like North if not at top floor

            for ip in range(1,N):
                for jp in range(1,N):
                    
                    ixx = (ip * N + jp)*4
                    ixxN = ixx + 0 #
                    
                    ixpN=(ip * N + jp+1)*4 #
                
                            
                    if jp < N-1:
                        # South of Upper is like North if not at top 
                        Afull[ixpN+2,jyy] = (Afull[ixpN+2,jyy] + Afull[ixxN,jyy])*1.
                        Afull[ixxN,jyy] = 0.
                    elif jp == N-1:
                        # if at top, add constant to B
                        Afull[ixxN,jyy] = 0.
                        if (i == ip) and (j==jp): 
                            B[jyy] -= fpred_aijb(psmooth[0],ip,jp,0.) 
    
    # Inverting (x-1) coefficients based on parity of number of displacements
    # between correlation and coefficient pixel
    
    for i in range(N):
        for j in range(N):
            
            if i ==0 and j==0: continue
            
            jyy = i * N + j -1 # correlation index
            
            for ip in range(0,N):
                for jp in range(0,N):
                    for ib in range(4):
                        ixb = (ip * N + jp)*4 + ib #
                        
                        nsteps = np.abs(ip-i)+np.abs(jp-j)
                        sign = (-1.)**float(nsteps % 2 != 0)
                        
                        if Afull[ixb,jyy] != 0.:
                            Afull[ixb,jyy] *= sign
                            
    
    # Prunning the A matrix
    
    A = np.zeros((N**2-1,N**2-1),dtype='float32')

    indepBounds = []
    
    ixp = 0
    Nraw = 4*N**2
    for ix in range(Nraw):
        if np.any(Afull[ix,:] != 0):
            A[ixp,:] = Afull[ix,:].copy()
            indepBounds.append(boundsList[ix])
            ixp += 1
    
    Nindep = len(indepBounds)
    
    print '\n Independent Coeffs. = %i / %i' % (Nindep,Nraw)
    print '\n A.shape = ',A.shape
    

    Atrans = A.transpose() 

    X = np.linalg.solve(Atrans,B)
    
    aijb = _build_aijb(X,indepBounds,psmooth)
    

    # verify matrix

    print '\nMatrix Verification'    
    res = Rij - aijb.sum(axis=2) # should be close to zero, except on (0,0)
    ii,jj = np.meshgrid(np.arange(N),np.arange(N),indexing='ij')
    sel = ((ii>0) | (jj>0))
    
    print 'mean res = %.1e' % res[sel].mean()
    print 'max res = %.1e' % res[sel].max()
    print 'min res = %.1e' % res[sel].min()
    
    if returnAll:
        return aijb,psmooth
    else:
        return aijb


def _build_aijb_synth(psmooth,N):
    """ 
    Builds aijb matrix using fpred(p,ijb)
    
    :param psmooth: 2-tuple, coefficients of fpred(p,ijb) function
    :param N: dimension of aijb matrix
    
    :return: aijb matrix
    
    """
    
    aijb = np.zeros((N,N,4),dtype='float32') + np.nan
    
    for ip in range(N):
        for jp in range(N):
            for x in range(4):
                aijb[ip,jp,x] = fpred_aijb(psmooth,ip,jp,x)
    
    return aijb


def get_kernel(aijb,writeFits=False):
    """ 
    
    'kernel' is an array (2N-1)x(2N-1)x4. Each plane kernel[:,:,b] is
    a 2D array with the displacement coefficients aijb, in all
    directions around a pixel at (0,0).
    
    :param aijb: array, matrix with displacements in 1st quadrant
    :param writeFits: save kernel to 4 FITS files
    
    :return: kernel matrix, (2N-1)x(2N-1)x4
    
    """
    from astropy.io import fits as fts

    N,Np,shouldbe4 = aijb.shape
    assert N == Np
    assert shouldbe4 == 4
    
    kernel = np.zeros((2*N-1,2*N-1,4),dtype='float32') + np.nan
    
    # North
    
    for i in range(N):
        for j in range(N):
            
            kernel[i+N-1,j+N-1,0] = aijb[i,j,2] # first quad
            if i>0 and j>0:
                kernel[i+N-1,N-j-1,0] = aijb[i,j,0] # second quad
                kernel[N-i-1,N-j-1,0] = aijb[i,j,0] # third quad                
                kernel[N-i-1,j+N-1,0] = aijb[i,j,2] # fourth quad
            elif i==0 and j>0:
                kernel[i+N-1,N-j-1,0] = aijb[i,j,0] # -y
            elif i>0 and j==0:
                kernel[N-i-1,j+N-1,0] = aijb[i,j,2] # -x
    
    # South
    
    for i in range(N):
        for j in range(N):
            
            kernel[i+N-1,j+N-1,2] = aijb[i,j,0] # first quad
            if i>0 and j>0:
                kernel[i+N-1,N-j-1,2] = aijb[i,j,2] # second quad
                kernel[N-i-1,N-j-1,2] = aijb[i,j,2] # third quad                
                kernel[N-i-1,j+N-1,2] = aijb[i,j,0] # fourth quad
            elif i==0 and j>0:
                kernel[i+N-1,N-j-1,2] = aijb[i,j,2] # -y
            elif i>0 and j==0:
                kernel[N-i-1,j+N-1,2] = aijb[i,j,0] # -x       

    # East
    
    for i in range(N):
        for j in range(N):
            
            kernel[i+N-1,j+N-1,1] = aijb[i,j,3] # first quad
            if i>0 and j>0:
                kernel[i+N-1,N-j-1,1] = aijb[i,j,3] # second quad
                kernel[N-i-1,N-j-1,1] = aijb[i,j,1] # third quad                
                kernel[N-i-1,j+N-1,1] = aijb[i,j,1] # fourth quad
            elif i==0 and j>0:
                kernel[i+N-1,N-j-1,1] = aijb[i,j,3] # -y
            elif i>0 and j==0:
                kernel[N-i-1,j+N-1,1] = aijb[i,j,1] # -x       

    # West

    for i in range(N):
        for j in range(N):
            
            kernel[i+N-1,j+N-1,3] = aijb[i,j,1] # first quad
            if i>0 and j>0:
                kernel[i+N-1,N-j-1,3] = aijb[i,j,1] # second quad
                kernel[N-i-1,N-j-1,3] = aijb[i,j,3] # third quad                
                kernel[N-i-1,j+N-1,3] = aijb[i,j,3] # fourth quad
            elif i==0 and j>0:
                kernel[i+N-1,N-j-1,3] =aijb[i,j,1] # -y
            elif i>0 and j==0:
                kernel[N-i-1,j+N-1,3] = aijb[i,j,3] # -x       
    
    if writeFits:
        fts.writeto('kernel_North.fits',kernel[:,:,0].transpose(),clobber=True)
        fts.writeto('kernel_East.fits',kernel[:,:,1].transpose(),clobber=True)
        fts.writeto('kernel_South.fits',kernel[:,:,2].transpose(),clobber=True)
        fts.writeto('kernel_West.fits',kernel[:,:,3].transpose(),clobber=True)
    
    kernel[:,:,:] = kernel[::-1,::-1,:] # surprise!
    
    # x - y permutation
    #kernel = kernel.transpose(1,0,2) # BEWARE of the Dog!
    
    return kernel
    

def get_Rdisp(img,aijb):
    """ Retrieves map of relative displacements of pixel boundaries, 
    for input img and Aijb matrix.
    

        See G15 - Eq. 6    
    
    :param img: image, 2D array
    :param aijb: aijb matrix, 3D array NxNx4
    
    :return: array, relative displacements all boundaries of pixels in img
    
    """
    
    kernel = get_kernel(aijb,writeFits=False)
    
    NX,NY = img.shape
    
    Rdisp = np.zeros((NX,NY,4),dtype='float32')
    
    #fts.writeto('kernel_North.fits',kernel[:,:,0],clobber=True)
    #fts.writeto('kernel_East.fits',kernel[:,:,1],clobber=True)
    #fts.writeto('kernel_South.fits',kernel[:,:,2],clobber=True)
    #fts.writeto('kernel_West.fits',kernel[:,:,3],clobber=True)
    
    
    for ixside in range(4):
        sidekernel = kernel[:,:,ixside].copy()
        Rdisp[:,:,ixside] =  1/2. * convolve2d(img,sidekernel,mode='same') # Eq. 6

    return Rdisp
    
    
    
def get_deltaQ(img,aijb,writeFits=False):
    """
    
    Retrieves deltaQ map for input image and aijb matrix.
    

    See G15 - Eq. 11
    
    :param img: image, 2D array
    :param aijb: Aijb matrix, 3D array
    :param writeFits: save FITS file with resulting dQ map (optional)
    
    :return: array, matrix with delta-Q for each pixel in img, given aijb
    
    """
    sign = 1. 
    
    kernel = get_kernel(aijb,writeFits)*sign
    
    dimg = np.zeros_like(img,dtype='float32')
    
    shift = [[1,0],[1,1],[-1,0],[-1,1]]
    
    for ixside in range(4):
        sidekernel = kernel[:,:,ixside].copy()
        dimg +=  (1/4. * (img + np.roll(img,shift[ixside][0],shift[ixside][1])) * \
             convolve2d(img,sidekernel,mode='same')) # Eq. 11

    return dimg
    

def degrade_estatic(img,aijb):
    """Degrades an image according to matrix of pixel-boundaries deformations.
    Follows on Eq. 11 of G15. Adds delta-Q.
    
    :param img: image, 2D array
    :param aijb: Aijb matrix, 3D array
    
    :return: array, img + delta-Q    
    
    """    
    dimg = get_deltaQ(img,aijb)    
    return img + dimg
    
def correct_estatic(img,aijb):
    """Corrects an image from pixel-boundaries deformation due to
    electrostatic forces. Subtracts delta-Q.
    
    :param img: image, 2D array
    :param aijb: Aijb matrix, 3D array
    
    :return: array, img - delta-Q    
    
    """
    dimg = get_deltaQ(img,aijb)
    return img - dimg
    


def plot_maps_ftheta(f,ii,jj,suptitle=''):
    """ """
    
    fig = plt.figure()
    
    i = ii.flatten().copy()
    j = jj.flatten().copy()
    
    f_N = np.zeros(len(i))
    f_E = np.zeros(len(i))
    f_S = np.zeros(len(i))
    f_W = np.zeros(len(i))
    
    for ix in range(len(i)):
        f_N[ix] = f(i[ix],j[ix],0)
        f_E[ix] = f(i[ix],j[ix],1)
        f_S[ix] = f(i[ix],j[ix],2)
        f_W[ix] = f(i[ix],j[ix],3)
        
    f_N = f_N.reshape(ii.shape)
    f_E = f_E.reshape(ii.shape)
    f_S = f_S.reshape(ii.shape)
    f_W = f_W.reshape(ii.shape)
    
    mini = i.min()
    maxi = i.max()
    minj = j.min()
    maxj = j.max()

    
    ax1 = fig.add_subplot(221)
    ax1.imshow(np.cos(f_N),cmap=plt.cm.Blues,interpolation='none',extent=[mini,maxi,minj,maxj],origin='lower left')
    ax1.set_title('N')
    
    ax2 = fig.add_subplot(222)
    ax2.imshow(np.cos(f_E),cmap=plt.cm.Greens,interpolation='none',extent=[mini,maxi,minj,maxj],origin='lower left')
    ax2.set_title('E')
    
    ax3 = fig.add_subplot(223)
    ax3.imshow(np.cos(f_S),cmap=plt.cm.Reds,interpolation='none',extent=[mini,maxi,minj,maxj],origin='lower left')
    ax3.set_title('S')
    
    ax4 = fig.add_subplot(224)
    ax4.imshow(np.cos(f_W),cmap=plt.cm.Oranges,interpolation='none',extent=[mini,maxi,minj,maxj],origin='lower left')
    ax4.set_title('W')
    
    plt.suptitle(suptitle)
    
    plt.show()
    
    plt.close()
    
    #stop()
    

def plot_map(z,ii,jj,title=''):
    """ """
    
    fig = plt.figure()
    
    i = ii.flatten().copy()
    j = jj.flatten().copy()
    
    mini = i.min()
    maxi = i.max()
    minj = j.min()
    maxj = j.max()

    
    ax1 = fig.add_subplot(111)
    cax1= ax1.imshow(z,cmap=plt.cm.Greys,interpolation='none',extent=[mini,maxi,minj,maxj],origin='lower left')
    cbar = fig.colorbar(cax1,orientation='vertical')
    ax1.set_title(title)
    plt.tight_layout()
    plt.show()
    plt.close()

def generate_GaussPSF(N,sigma):
    """ 
    Create a circular symmetric Gaussian centered on the centre of
    a NxN matrix/image.
    
    """
    
    x = (N-1.)/2.
    y = (N-1.)/2.
    
    #x and y coordinate vectors
    Gyvect = np.arange(0, N )
    Gxvect = np.arange(0, N )

    #meshgrid
    Gxmesh, Gymesh = np.meshgrid(Gxvect, Gyvect)

    #normalizers
    sigmax = 1. / (2. * sigma**2)
    sigmay = sigmax #same sigma in both directions, thus same normalizer

    #gaussian
    exponent = (sigmax * (Gxmesh - x)**2 + sigmay * (Gymesh - y)**2)
    Gaussian = np.exp(-exponent) / (2. * np.pi * sigma*sigma)
    

    #normalize to unity
    Gaussian /= np.max(Gaussian)

    return Gaussian


def show_disps_CCD273(aijb,stretch=5.,peak=1.E5/3.5,N=25,sigma=1.6,title='',figname=''):
    """ """
    
    sign = 1. 
    
    gaussPSF = generate_GaussPSF(N,sigma)
    gaussPSF = gaussPSF / gaussPSF.max() * peak
    
    
    Rdisp_Q = get_Rdisp(gaussPSF,aijb) / peak # ??
    
    #stop()
                       
    x0 = -(N-1)/2.*pixpitch
    x1 = (N-1)/2. * pixpitch
    y0 = -(N-1)/2.*pixpitch
    y1 = (N-1)/2.*pixpitch
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(gaussPSF,origin='lower left',extent=(x0,x1,y0,y1))
    ax.set_title(title)
    
    
    for i in range(N):
        for j in range(N):
            
            c = 'w'
            #if ((i % 2 == 0) and (j % 2 == 0)) or ((i % 2 !=0) and (j % 2 !=0)):
            #       c = 'g'
            
            d = Rdisp_Q[i,j,:] * stretch * pixpitch  # N,E,S,W
            
            # ll, ul, ur, lr
            # displacements are inside-outside oriented
            x = (i-(N-1.)/2.)*pixpitch + (np.array([-0.5,-0.5,0.5,0.5,-0.5]))*pixpitch +\
                sign * np.array([-d[3],-d[3],d[1],d[1],-d[3]])
                
            y = (j-(N-1.)/2.)*pixpitch + (np.array([-0.5,0.5,0.5,-0.5,-0.5]))*pixpitch +\
                sign * np.array([-d[2],d[0],d[0],-d[2],-d[2]])
            
            ax.plot(x,y,'%s-' % c)
    
    ax.set_xlim((-5.*pixpitch,5.*pixpitch))
    ax.set_ylim((-5.*pixpitch,5.*pixpitch))
    
    if figname != '':
        plt.savefig(figname)
    else:
        plt.show()
    

def test0():
    """ """
    
    N = 5
    sigma_covij = 0.5
    p = (1.E-3,1.E-3)
    
    ii,jj = np.meshgrid(np.arange(N),np.arange(N),indexing='xy')
    rr = np.sqrt(ii**2.+jj**2.)
    covij = np.exp(-(rr**2./(2*sigma_covij**2.))) + 0.0005*sigma_covij
    
    
    i = ii.flatten().copy()
    j = jj.flatten().copy()
    
    fpred_N = np.zeros(len(i))
    fpred_E = np.zeros(len(i))
    fpred_S = np.zeros(len(i))
    fpred_W = np.zeros(len(i))
    
    for ix in range(len(i)):
        fpred_N[ix] = fpred_aijb(p,i[ix],j[ix],0)
        fpred_E[ix] = fpred_aijb(p,i[ix],j[ix],1)
        fpred_S[ix] = fpred_aijb(p,i[ix],j[ix],2)
        fpred_W[ix] = fpred_aijb(p,i[ix],j[ix],3)
        
    fpred_N = fpred_N.reshape(ii.shape)
    fpred_E = fpred_E.reshape(ii.shape)
    fpred_S = fpred_S.reshape(ii.shape)
    fpred_W = fpred_W.reshape(ii.shape)
    
    
    plot_maps_ftheta(ftheta_bound,ii,jj,r'$\theta_B$')
    plot_maps_ftheta(frdist,ii,jj,'f(r)')
    plot_map(np.log10(covij),ii,jj,title=r'$Cov_{ij}$')
    
    plot_map(fpred_N,ii,jj,title=r'$fpred_N$')
    plot_map(fpred_E,ii,jj,title=r'$fpred_E$')
    plot_map(fpred_S,ii,jj,title=r'$fpred_S$')
    plot_map(fpred_W,ii,jj,title=r'$fpred_W$')


def test_solve():
    """ """
    N=5
    sigma_covij = 0.5
    var = sigma_covij**2.
    mu = 1.
    
    ii,jj = np.meshgrid(np.arange(N),np.arange(N),indexing='xy')
    rr = np.sqrt(ii**2.+jj**2.)
    covij = sigma_covij**2. * mu * np.exp(-(rr**2./(2*(1.3)**2.))) #+ sigma_covij * 0.0005

    psmooth,epsmooth = solve_for_psmooth(covij,var,mu)
    
    Abest = solve_for_A_linalg(covij,var=1.,mu=1.)




def test_selfconsist():
    """ """
    N=3
    var = 1.
    mu = 1.
    
    #pin = [0.00237274,  0.13478362]
    pin = [0.0237274, -1.2]
    
    #ii,jj = np.meshgrid(np.arange(N),np.arange(N),indexing='xy')
    
    covij = np.zeros((N,N,4),dtype='float32')
    
    for i in range(N):
        for j in range(N):
            covij[i,j,0] = fpred_aijb(pin,i,j,0.)
            covij[i,j,1] = fpred_aijb(pin,i,j,1.)
            covij[i,j,2] = fpred_aijb(pin,i,j,2.)
            covij[i,j,3] = fpred_aijb(pin,i,j,3.)
    
    covij = covij.sum(axis=2)
    
    showcovij = False
    
    if showcovij:
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(covij,origin='lower left')
        cax.set_cmap('spectral')
        cbar = fig.colorbar(cax,orientation='vertical')
        ax.set_title('Cov_ij')
        ax.set_xlabel('i')
        ax.set_ylabel('j')
        plt.show()
    

    Asynth = _build_aijb_synth(pin,N)
    #Abest = solve_for_A_linalg(covij,var=var,mu=mu,doplot=False)
    Abest = solve_for_A_linalg(covij,var=var,mu=mu,doplot=False) #,psmooth=[pin])
    
    
    show_disps_CCD273(Asynth,stretch=30.,peak=1.E5/3.5,N=25,sigma=1.6,title='Input')
    show_disps_CCD273(Abest,stretch=30.,peak=1.E5/3.5,N=25,sigma=1.6,title='Output')
    
    singlepixmap = np.zeros((101,101),dtype='float32') + 0.01
    singlepixmap[50,50] = 1.
    
    k = get_kernel(Asynth,writeFits=True)
    
    
    kernel = degrade_estatic(singlepixmap,Abest)
    deltaQ = get_deltaQ(singlepixmap,Abest)
    kernelsynth = degrade_estatic(singlepixmap,Asynth)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax1.imshow(kernel)
    ax1.set_title('KERNEL')
    ax2 = fig.add_subplot(122)
    ax2.imshow(kernelsynth)
    ax2.set_title('KERNEL-SYNTH')
    plt.show()    
    
    
    fts.writeto('kernel_selfconsist.fits',kernel.transpose(),clobber=True)
    fts.writeto('deltaQ_selfconsist.fits',deltaQ.transpose(),clobber=True)
    fts.writeto('kernel_selfconsist_synth.fits',kernelsynth.transpose(),clobber=True)
    
    stop()

def test_getkernel():
    """ """
    
    xx,yy = np.meshgrid(np.arange(5),np.arange(5),indexing='xy')
    
    aijb = np.zeros((5,5,4),dtype='float32')
    
    for i in range(4):
        aijb[:,:,i] = (xx**2.+yy**2.)**0.5
    
    
    kernel = get_kernel(aijb,writeFits=True)
    
    
    
if __name__ == '__main__':
    
    #test_getkernel() 
    
    #test0()
    
    #test_solve()
    test_selfconsist()
