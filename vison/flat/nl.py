#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

NEEDSREVISION

Module with tools used in NL analysis.

Created on Mon Feb 5 15:51:00 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

#IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import datetime
#END IMPORT

FullDynRange = 2.**16
NLdeg=6

def get_exptime_atmiddynrange(flu1D,exp1D):
    """ """
    mod1d_fit = np.polyfit(flu1D/FullDynRange,exp1D,1) # a linear approx. is fine
    mod1d_pol = np.poly1d(mod1d_fit)
    t50 = mod1d_pol(0.5)
    return t50
    

def fitNL(fluences,exptimes):
    """ """
    
    assert fluences.shape[0] == exptimes.shape[0]
    assert fluences.ndim <= 2
    assert exptimes.ndim == 1
    
    Nexp = len(exptimes)

    if fluences.ndim == 2:
        Nsec = fluences.shape[1]
        #_exptimes = np.repeat(exptimes.reshape(Nexp,1),Nsec,axis=1)
        
        t50 = np.zeros(Nsec,dtype='float32') + np.nan
                      
        for i in range(Nsec):            
            t50[i] = get_exptime_atmiddynrange(fluences[:,i],exptimes)
            
        t50 = np.repeat(t50.reshape(1,Nsec),Nexp,axis=0)
    
        exptimes_bc = np.repeat(exptimes.reshape(Nexp,1),Nsec,axis=1)
        
    else:
        
        t50 = get_exptime_atmiddynrange(fluences[:,i],exptimes)
        
        exptimes_bc = exptimes.copy()
        
    YL = exptimes_bc/t50*FullDynRange/2.
    
    Z = 100.*(YL/fluences-1.)
    
    NLfit = np.polyfit(YL.flatten(),Z.flatten(),deg=NLdeg,full=False)
    
    NLpol = np.poly1d(NLfit)
    
    xNL = np.arange(1,FullDynRange,dtype='float32') 
    # array with NL fluences (1 to 2**16, in steps of ADU)
    
    ixmax = NLpol(xNL).argmax()
    maxNLpc = NLpol[ixmax]
    flu_maxNLpc = xNL[ixmax]
    
    fitresults = dict(coeffs=NLfit,NLdeg=NLdeg,maxNLpc=maxNLpc,
                      flu_maxNLpc=flu_maxNLpc)
    return fitresults


def wrap_fitNL(raw_data,exptimes,col_labels,times=np.array([]),TrackFlux=True,subBgd=True):
    """ """
    # col1 == BGD
    # colEVEN = STAB
    # colODD = Fluences != 0
    
    NObsIDs, Nsecs = raw_data.shape
    
    dtimes = np.array([(times[i]-times[0]).seconds for i in range(NObsIDs)],dtype='float32')
    
    col_numbers = np.array([int(item[3:]) for item in col_labels])
    
    ixboo_bgd = col_numbers == 1
    ixboo_stab = (col_numbers % 2 != 0) & (col_numbers>1) # ODD, > 1
    ixboo_fluences = col_numbers % 2 == 0 # EVEN
    
    if subBgd:
        bgd = np.nanmean(raw_data[ixboo_bgd,:],axis=0)
        bgd = np.repeat(bgd.reshape(1,Nsecs),NObsIDs,axis=0).copy()
        raw_data -= bgd
    
    if TrackFlux and len(times) == NObsIDs:
        
        st_dtimes = dtimes[ixboo_stab].copy()
        st_fluences = np.nanmean(raw_data[ixboo_stab,:],axis=1).copy()
        
        track_fit = np.polyfit(st_dtimes,st_fluences,2,full=False,cov=False)
        track_pol = np.poly1d(track_fit)
        
        track = track_pol(dtimes[ixboo_fluences])
        track /= np.median(track)
        track_bc = np.repeat(track.reshape(len(track),1),Nsecs,axis=1)
        
        raw_data[ixboo_fluences,:] /= track_bc
                
    
    fitresults = fitNL(raw_data[ixboo_fluences,:],exptimes[ixboo_fluences],)
    
    
    
    return fitresults
    

def test_wrap_fitNL():
    """ """
    
    col_labels = np.array(['col1', 'col1', 'col1', 'col1', 
       'col2', 'col2', 'col2', 'col2','col2', 
       'col3', 
       'col4', 'col4', 'col4', 'col4', 'col4', 
       'col5',
       'col6', 'col6', 'col6', 'col6', 'col6', 
       'col7', 
       'col8', 'col8','col8', 'col8', 'col8', 
       'col9', 
       'col10', 'col10', 'col10', 'col10','col10', 
       'col13', 
       'col14', 'col14', 'col14', 'col14', 'col14',
       'col15', 
       'col16', 'col16', 'col16', 'col16', 'col16', 
       'col17',
       'col18', 'col18', 'col18', 'col18', 'col18', 
       'col19', 
       'col20','col20', 'col20', 'col20', 'col20', 
       'col21', 
       'col22', 'col22','col22', 'col22', 'col22', 
       'col23', 
       'col24', 'col24', 'col24','col24', 'col24', 
       'col25'], 
      dtype='|S12')

    exptimes = np.array([  0. ,   0. ,   0. ,   0. ,   1.5,   1.5,   1.5,   1.5,   1.5,
        15. ,   3. ,   3. ,   3. ,   3. ,   3. ,  15. ,   6. ,   6. ,
         6. ,   6. ,   6. ,  15. ,   9. ,   9. ,   9. ,   9. ,   9. ,
        15. ,  15. ,  15. ,  15. ,  15. ,  15. ,  15. ,  21. ,  21. ,
        21. ,  21. ,  21. ,  15. ,  24. ,  24. ,  24. ,  24. ,  24. ,
        15. ,  27. ,  27. ,  27. ,  27. ,  27. ,  15. ,  30. ,  30. ,
        30. ,  30. ,  30. ,  15. ,  33. ,  33. ,  33. ,  33. ,  33. ,
        15. ,  36. ,  36. ,  36. ,  36. ,  36. ,  15. ])
    
    raw_data = np.repeat((exptimes/exptimes.max()*2.**16).reshape((70,1)),49,axis=1)
    dtobjs = np.array([datetime.datetime(2018, 2, 8, 3, 59, 17),
       datetime.datetime(2018, 2, 8, 4, 0, 55),datetime.datetime(2018, 2, 8, 4, 2, 28),
       datetime.datetime(2018, 2, 8, 4, 4, 7),datetime.datetime(2018, 2, 8, 4, 5, 43),
       datetime.datetime(2018, 2, 8, 4, 7, 19),datetime.datetime(2018, 2, 8, 4, 8, 54),
       datetime.datetime(2018, 2, 8, 4, 10, 29),datetime.datetime(2018, 2, 8, 4, 12, 5),
       datetime.datetime(2018, 2, 8, 4, 14, 2),datetime.datetime(2018, 2, 8, 4, 15, 40),
       datetime.datetime(2018, 2, 8, 4, 17, 18),datetime.datetime(2018, 2, 8, 4, 18, 55),
       datetime.datetime(2018, 2, 8, 4, 20, 38),datetime.datetime(2018, 2, 8, 4, 22, 15),
       datetime.datetime(2018, 2, 8, 4, 24, 5),datetime.datetime(2018, 2, 8, 4, 25, 46),
       datetime.datetime(2018, 2, 8, 4, 27, 27),datetime.datetime(2018, 2, 8, 4, 29, 6),
       datetime.datetime(2018, 2, 8, 4, 30, 47),datetime.datetime(2018, 2, 8, 4, 32, 28),
       datetime.datetime(2018, 2, 8, 4, 34, 18),datetime.datetime(2018, 2, 8, 4, 36, 2),
       datetime.datetime(2018, 2, 8, 4, 37, 46),datetime.datetime(2018, 2, 8, 4, 39, 28),
       datetime.datetime(2018, 2, 8, 4, 41, 19),datetime.datetime(2018, 2, 8, 4, 43, 11),
       datetime.datetime(2018, 2, 8, 4, 45),datetime.datetime(2018, 2, 8, 4, 46, 57),
       datetime.datetime(2018, 2, 8, 4, 48, 46),datetime.datetime(2018, 2, 8, 4, 50, 34),
       datetime.datetime(2018, 2, 8, 4, 52, 24),datetime.datetime(2018, 2, 8, 4, 54, 13),
       datetime.datetime(2018, 2, 8, 4, 56, 10),datetime.datetime(2018, 2, 8, 4, 58, 6),
       datetime.datetime(2018, 2, 8, 5, 0, 5),datetime.datetime(2018, 2, 8, 5, 2),
       datetime.datetime(2018, 2, 8, 5, 4, 2),datetime.datetime(2018, 2, 8, 5, 5, 58),
       datetime.datetime(2018, 2, 8, 5, 7, 48),datetime.datetime(2018, 2, 8, 5, 9, 48),
       datetime.datetime(2018, 2, 8, 5, 11, 46),datetime.datetime(2018, 2, 8, 5, 13, 43),
       datetime.datetime(2018, 2, 8, 5, 15, 42),datetime.datetime(2018, 2, 8, 5, 17, 40),
       datetime.datetime(2018, 2, 8, 5, 19, 30),datetime.datetime(2018, 2, 8, 5, 21, 40),
       datetime.datetime(2018, 2, 8, 5, 23, 46),datetime.datetime(2018, 2, 8, 5, 25, 46),
       datetime.datetime(2018, 2, 8, 5, 27, 47),datetime.datetime(2018, 2, 8, 5, 29, 49),
       datetime.datetime(2018, 2, 8, 5, 31, 39),datetime.datetime(2018, 2, 8, 5, 33, 44),
       datetime.datetime(2018, 2, 8, 5, 35, 49),datetime.datetime(2018, 2, 8, 5, 37, 52),
       datetime.datetime(2018, 2, 8, 5, 39, 56),datetime.datetime(2018, 2, 8, 5, 42, 8),
       datetime.datetime(2018, 2, 8, 5, 43, 58),datetime.datetime(2018, 2, 8, 5, 46, 7),
       datetime.datetime(2018, 2, 8, 5, 48, 14),datetime.datetime(2018, 2, 8, 5, 50, 26),
       datetime.datetime(2018, 2, 8, 5, 52, 34),datetime.datetime(2018, 2, 8, 5, 54, 43),
       datetime.datetime(2018, 2, 8, 5, 56, 42),datetime.datetime(2018, 2, 8, 5, 58, 58),
       datetime.datetime(2018, 2, 8, 6, 1, 8),datetime.datetime(2018, 2, 8, 6, 3, 18),
       datetime.datetime(2018, 2, 8, 6, 5, 29),datetime.datetime(2018, 2, 8, 6, 7, 40),
       datetime.datetime(2018, 2, 8, 6, 9, 33)], dtype=object)

    
    fitresults = wrap_fitNL(raw_data,exptimes,col_labels,times=dtobjs,TrackFlux=True,subBgd=True)
    

if __name__ == '__main__':
    test_wrap_fitNL()
    