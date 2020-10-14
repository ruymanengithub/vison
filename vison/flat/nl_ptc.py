#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""


Module with tools used in NL analysis, from PTC meta-analysis.

Created on Wed Oct 14 16:40:00 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
from scipy import stats

# END IMPORT

def get_mse_var(samples):
    """ 
    https://en.wikipedia.org/wiki/Mean_squared_error#Variance

    'In statistics, the mean squared error (MSE)[1][2] or mean squared deviation (MSD) of 
    an estimator (of a procedure for estimating an unobserved quantity) measures the 
    average of the squares of the errors—that is, the average squared difference between 
    the estimated values and the actual value. 

    The MSE is a measure of the quality of an estimator —it is always non-negative, and 
    values closer to zero are better.'
    
    MSE(S2_n-1) = 1/n (gamma_2 + 2n/(n-1)) sigma_4
     
    """
    fsamples = samples.flatten()
    n = fsamples.size
    sigma = np.nanstd(fsamples,ddof=1)
    mu4 = stats.moment(fsamples,moment=4)
    gamma2 = mu4 / sigma**4.-3. # excess kurtosis
    
    MSEvar = 1./n * (gamma2 + 2*n/(n-1) ) * sigma**4.
    
    return MSEvar


def f_extract_PTC(self, ccdobjcol, medcol, varcol, evarcol, binfactor=1):
        """ """

    # HARDWIRED VALUES
    wpx = self.window['wpx']
    hpx = self.window['hpx']

    indices = copy.deepcopy(self.dd.indices)

    nObs, nCCD, nQuad = indices.shape[0:3]

    Quads = indices.get_vals('Quad')
    CCDs = indices.get_vals('CCD')

    tile_coos = dict()
    for Quad in Quads:
        tile_coos[Quad] = self.ccdcalc.get_tile_coos(Quad, wpx, hpx)
    Nsectors = tile_coos[Quads[0]]['Nsamps']
    sectornames = np.arange(Nsectors)

    Sindices = copy.deepcopy(self.dd.indices)
    if 'Sector' not in Sindices.names:
        Sindices.append(core.vIndex('Sector', vals=sectornames))

    # Initializing new columns

    valini = 0.
    self.dd.initColumn(medcol, Sindices, dtype='float32', valini=valini)
    self.dd.initColumn(varcol, Sindices, dtype='float32', valini=valini)
    self.dd.initColumn(evarcol, Sindices, dtype='float32', valini=valini)

    # labels should be the same accross CCDs. PATCH.
    label = self.dd.mx['label'][:, 0].copy()
    ulabels = np.unique(label)
    ObsIDs = self.dd.mx['ObsID'][:].copy()

    # Pairing ObsIDs

    self.dd.initColumn(
        'ObsID_pair', self.dd.mx['ObsID'].indices, dtype='int64', valini=0)


    for ulabel in ulabels:
        six = np.where(label == ulabel)
        nsix = len(six[0])
        if nsix % 2 ==0:
            ixeven = np.arange(0, nsix, 2)
            ixodd = np.arange(1, nsix, 2)

            self.dd.mx['ObsID_pair'][six[0][ixeven]] = ObsIDs[six[0][ixodd]]

    if not self.drill:

        # if self.proc_histo['Masked']:
        #    estimators = dict(median=np.ma.median,std=np.ma.std)
        # else:
        #    estimators = dict(median=np.median,std=np.std)

        dpath = self.inputs['subpaths']['ccdpickles']

        misspairs = []

        for iObs in range(nObs):

            _ObsID_pair = self.dd.mx['ObsID_pair'][iObs]
            if _ObsID_pair == 0:
                continue
            iObs_pair = np.where(ObsIDs == _ObsID_pair)[0][0]
            

            flu_i = self.dd.mx['flu_med_img'][iObs, ...].mean()
            flu_p = self.dd.mx['flu_med_img'][iObs_pair, ...].mean()

            flus = np.array([flu_i, flu_p])

            if flus.std() / flus.mean() > 0.1:

                self.dd.mx[medcol][iObs, ...] = np.nan
                self.dd.mx[varcol][iObs, ...] = np.nan

                misspairs.append((self.dd.mx['ObsID'][iObs], self.dd.mx['ObsID'][iObs_pair]))

                continue

            for jCCD, CCDk in enumerate(CCDs):

                ccdobj_odd_f = os.path.join(
                    dpath, '%s.pick' % self.dd.mx[ccdobjcol][iObs, jCCD])
                ccdobj_eve_f = os.path.join(
                    dpath, '%s.pick' % self.dd.mx[ccdobjcol][iObs_pair, jCCD])

                if (not os.path.exists(ccdobj_odd_f)) or \
                    (not os.path.exists(ccdobj_eve_f)):
                    continue

                ccdobj_odd = copy.deepcopy(
                    cPickleRead(ccdobj_odd_f))
                ccdobj_eve = copy.deepcopy(
                    cPickleRead(ccdobj_eve_f))

                evedata = ccdobj_eve.extensions[-1].data.copy()

                # easy way to subtract one image from the other
                ccdobj_sub = copy.deepcopy(ccdobj_odd)
                ccdobj_sub.sub_bias(evedata, extension=-1)

                for kQ in range(nQuad):

                    Quad = Quads[kQ]

                    _tile_coos = tile_coos[Quad]

                    if binfactor >1:

                        ccdobj_odd = PTC0Xaux.CCDclone(ccdobj_odd)

                        _meds = ccdobj_odd.get_tiles_stats(
                            Quad, _tile_coos, 'median', extension=-1, binfactor=binfactor)

                        # IT'S A SUBTRACTION, SO WE HAVE TO DIVIDE BY 2 THE VARIANCE!

                        ccdobj_sub = PTC0Xaux.CCDclone(ccdobj_sub)

                        _vars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, 'std', extension=-1, binfactor=binfactor)**2. / 2.

                        _evars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, statkey=None, 
                            extractor=f_extract_mse_var,
                            extension=-1, binfactor=binfactor) / 2.**2.

                    else:
                        _meds = ccdobj_odd.get_tiles_stats(
                            Quad, _tile_coos, 'median', extension=-1)

                        # IT'S A SUBTRACTION, SO WE HAVE TO DIVIDE BY 2 THE VARIANCE!
                        _vars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, 'std', extension=-1)**2. / 2.

                        _evars = ccdobj_sub.get_tiles_stats(
                            Quad, _tile_coos, statkey=None, 
                            extractor=f_extract_mse_var,
                            extension=-1) / 2.**2.

                    self.dd.mx[medcol][iObs, jCCD, kQ, :] = _meds.copy()
                    self.dd.mx[varcol][iObs, jCCD, kQ, :] = _vars.copy()
                    self.dd.mx[evarcol][iObs, jCCD, kQ, :] = _evars.copy()
        
        if len(misspairs) > 0:
            if self.report is not None:
                self.report.add_Text(
                    'Pairs with unequal fluence skipped: %s' %
                    misspairs.__repr__())
            if self.log is not None:
                self.log.info('Pairs with unequal fluence skipped: %s' % \
                    misspairs.__repr__())    