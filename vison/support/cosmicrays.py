"""
Cosmic Rays
===========

This simple class can be used to include cosmic ray events to an image.
By default the cosmic ray events are drawn from distributions describing
the length and energy of the events. Such distributions can be generated
for example using Stardust code (http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=04636917).
The energy of the cosmic ray events can also be set to constant for
testing purposes. The class can be used to draw a single cosmic ray
event or up to a covering fraction.

:requires: NumPy
:requires: SciPy

:version: 0.25

:author: Sami-Matias Niemi, R. Azzollini
:contact: r.azzollini@ucl.ac.uk


:History:
0.21 - corrected a typo in _cosmicRayIntercepts (nested loop used index of nesting loop)
0.25 - added addToFluxTime
0.3 (3/Jul/2017) - corrected bug, as lum. of CRs was doubled by mistake

"""

from pdb import set_trace as stop

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import interp1d


class cosmicrays():
    """
    Cosmic ray generation class. Can either draw events from distributions or
    set the energy of the events to a constant.

    :param log: logger instance
    :param image: image to which cosmic rays are added to (a copy is made not to change the original numpy array)
    :param crInfo: column information (cosmic ray file)
    :param information: cosmic ray track information (file containing track length and energy information) and
                        exposure time.
    """
    def __init__(self, log, image, crInfo=None, information=None):
        """
        Cosmic ray generation class. Can either draw events from distributions or
        set the energy of the events to a constant.

        :param log: logger instance
        :param image: image to which cosmic rays are added to (a copy is made not to change the original numpy array)
        :param crInfo: column information (cosmic ray file)
        :param information: cosmic ray track information (file containing track length and energy information) and
                            exposure time.
        """
        #setup logger
        self.log = log
        
        #image and size
        self.image = image.copy()
        self.ysize, self.xsize = self.image.shape

        #set up the information dictionary, first with defaults and then overwrite with inputs if given
        self.information = (dict(#cosmicraylengths='data/cdf_cr_length.dat',
                                 cosmicraylength='data/TableauHistogLgTraj_CDF.txt',
                                 cosmicraydistance='data/cdf_cr_total.dat',
                                 exptime=565.))
        if information is not None:
            self.information.update(information)

        if crInfo is not None:
            self.cr = crInfo
        else:
            self._readCosmicrayInformation()


    def _readCosmicrayInformation(self):
        self.log.info('Reading in cosmic ray information from %s and %s' % (self.information['cosmicraylengths'],
                                                                            self.information['cosmicraydistance']))
        #read in the information from the files
        crLengths = np.loadtxt(self.information['cosmicraylengths'])
        crDists = np.loadtxt(self.information['cosmicraydistance'])

        #set up the cosmic ray information dictionary
        self.cr = dict(cr_u=crLengths[:, 0], cr_cdf=crLengths[:, 1], cr_cdfn=np.shape(crLengths)[0],
                       cr_v=crDists[:, 0], cr_cde=crDists[:, 1], cr_cden=np.shape(crDists)[0])

        return self.cr


    def _cosmicRayIntercepts(self, lum, x0, y0, l, phi):
        """
        Derive cosmic ray streak intercept points.

        :param lum: luminosities of the cosmic ray tracks
        :param x0: central positions of the cosmic ray tracks in x-direction
        :param y0: central positions of the cosmic ray tracks in y-direction
        :param l: lengths of the cosmic ray tracks
        :param phi: orientation angles of the cosmic ray tracks

        :return: cosmic ray map (image)
        :rtype: nd-array
        """
        #create empty array
        crImage = np.zeros((self.ysize, self.xsize), dtype=np.float64)

        #x and y shifts
        dx = l * np.cos(phi) / 2. # beware! 0<phi< pi, dx < 0
        dy = l * np.sin(phi) / 2.
        mskdx = np.abs(dx) < 1e-8
        mskdy = np.abs(dy) < 1e-8
        dx[mskdx] = 0.
        dy[mskdy] = 0.

        #pixels in x-direction
        ilo = np.round(x0.copy() - dx)
        msk = ilo < 0.
        ilo[msk] = 0
        ilo = ilo.astype(np.int)

        #ihi = 1 + np.round(x0.copy() + dx) # <3/Jul/2017
        ihi = np.round(x0.copy() + dx)     # >3/Jul/2017   
        msk = ihi > self.xsize
        ihi[msk] = self.xsize
        ihi = ihi.astype(np.int)

        #pixels in y-directions
        jlo = np.round(y0.copy() - dy)
        msk = jlo < 0.
        jlo[msk] = 0
        jlo = jlo.astype(np.int)

        #jhi = 1 + np.round(y0.copy() + dy) < 3/Jul/2017
        jhi = np.round(y0.copy() + dy)      # > 3/Jul/2017
        msk = jhi > self.ysize
        jhi[msk] = self.ysize
        jhi = jhi.astype(np.int)
        
        
        #stop()
        #v_n = []
        #v_e = []
        #v_w = []
        
        #offending_delta = 0.5
        offending_delta = 1.
        
        #loop over the individual events
        for i, luminosity in enumerate(lum):
            n = 0  # count the intercepts

            u = []
            x = []
            y = []

            #Compute X intercepts on the pixel grid
            if ilo[i] < ihi[i]:
                for xcoord in xrange(ilo[i], ihi[i]):
                    ok = (xcoord - x0[i]) / dx[i]
                    if np.abs(ok) <= offending_delta:
                        n += 1
                        u.append(ok)
                        x.append(xcoord)
                        y.append(y0[i] + ok * dy[i])
            else:
                for xcoord in xrange(ihi[i], ilo[i]):
                    ok = (xcoord - x0[i]) / dx[i] # parametric representation of line
                    if np.abs(ok) <= offending_delta:
                        n += 1
                        u.append(ok)
                        x.append(xcoord)
                        y.append(y0[i] + ok * dy[i])

            #Compute Y intercepts on the pixel grid
            if jlo[i] < jhi[i]:
                for ycoord in xrange(jlo[i], jhi[i]):
                    ok = (ycoord - y0[i]) / dy[i]
                    if np.abs(ok) <= offending_delta:
                        n += 1
                        u.append(ok)
                        x.append(x0[i] + ok * dx[i])
                        y.append(ycoord)
            else:
                for ycoord in xrange(jhi[i], jlo[i]):
                    ok = (ycoord - y0[i]) / dy[i]
                    if np.abs(ok) <= offending_delta:
                        n += 1
                        u.append(ok)
                        x.append(x0[i] + ok * dx[i])
                        y.append(ycoord)

            #check if no intercepts were found
            if n < 1:
                xc = int(np.floor(x0[i]))
                yc = int(np.floor(y0[i]))
                crImage[yc, xc] += luminosity
                
                #v_n += [1]
                #v_e += [luminosity]
                #v_w += [1.]

            #Find the arguments that sort the intersections along the track
            u = np.asarray(u)
            x = np.asarray(x)
            y = np.asarray(y)

            args = np.argsort(u)

            u = u[args]
            x = x[args]
            y = y[args]
            
            #s_n = 0.
            #s_e = 0.
            #s_w = 0.
            
            
            #Decide which cell each interval traverses, and the path length
            # AZZO: beware! index i was j before! Guessed it was a typo.
            for j in xrange(1, n - 1):
                #w = u[j + 1] - u[j] # <3/Jul/2017
                w = (u[j + 1] - u[j])/2. # >3/Jul/2017
                cx = int(1 + np.floor((x[j + 1] + x[j]) / 2.))
                cy = int(1 + np.floor((y[j + 1] + y[j]) / 2.))

                if 0 <= cx < self.xsize and 0 <= cy < self.ysize:
                    crImage[cy, cx] += (w * luminosity)
                    
                    #s_n += 1
                    #s_e += w*luminosity
                    #s_w += w
            
            # if s_n>0:
            #    v_n += [s_n]
            #    v_e += [s_e]
            #    v_w += [s_w]

        
        return crImage


    def _drawCosmicRays(self, limit=None):
        """
        Add cosmic rays to the arrays based on a power-law intensity distribution for tracks.

        Cosmic ray properties (such as location and angle) are chosen from random Uniform distribution.
        """
        #estimate the number of cosmics
        cr_n = self.xsize * self.ysize * 0.014 / 43.263316 * 2. # AZZO: MEANING of constants??
        #scale with exposure time, the above numbers are for the nominal 565s exposure
        cr_n *= (self.information['exptime'] / 565.0)
        
        #assume a power-law intensity distribution for tracks
        fit = dict(cr_lo=1.0e3, cr_hi=1.0e5, cr_q=2.0e0)
        fit['q1'] = 1.0e0 - fit['cr_q']
        fit['en1'] = fit['cr_lo'] ** fit['q1']
        fit['en2'] = fit['cr_hi'] ** fit['q1']

        #pseudo-random numbers taken from a uniform distribution between 0 and 1
        luck = np.random.rand(int(np.floor(cr_n)))

        #draw the length of the tracks
        if self.cr['cr_cdfn'] > 1:
            ius = InterpolatedUnivariateSpline(self.cr['cr_cdf'], self.cr['cr_u'])
            self.cr['cr_l'] = ius(luck)
        else:
            self.cr['cr_l'] = np.sqrt(1.0 - luck ** 2) / luck
        

        #draw the energy of the tracks
        if self.cr['cr_cden'] > 1:
            ius = InterpolatedUnivariateSpline(self.cr['cr_cde'], self.cr['cr_v'])
            self.cr['cr_e'] = ius(luck)
        else:
            self.cr['cr_e'] = (fit['en1'] + (fit['en2'] - fit['en1']) *
                               np.random.rand(int(np.floor(cr_n)))) ** (1.0 / fit['q1'])

        #Choose the properties such as positions and an angle from a random Uniform dist
        cr_x = self.xsize * np.random.rand(int(np.floor(cr_n)))
        cr_y = self.ysize * np.random.rand(int(np.floor(cr_n)))
        cr_phi = np.pi * np.random.rand(int(np.floor(cr_n)))

        #find the intercepts
        if limit is None:
            self.cosmicrayMap = self._cosmicRayIntercepts(self.cr['cr_e'], cr_x, cr_y, self.cr['cr_l'], cr_phi)
            print('Number of cosmic ray events:', len(self.cr['cr_e']))
        else:
            #limit to electron levels < limit
            msk = self.cr['cr_e'] < limit
            print('Number of cosmic ray events: %i / %i' % (len(self.cr['cr_e'][msk]), int(np.floor(cr_n))))
            self.cosmicrayMap = self._cosmicRayIntercepts(self.cr['cr_e'][msk], cr_x[msk], cr_y[msk],
                                                          self.cr['cr_l'][msk], cr_phi[msk])

        #count the covering factor
        area_cr = np.count_nonzero(self.cosmicrayMap)
        text = 'The cosmic ray covering factor is %i pixels i.e. %.3f per cent' \
               % (area_cr, 100.*area_cr / (self.xsize*self.ysize))
        self.log.info(text)
        
        print(text)


    def _drawSingleEvent(self, limit=1000, cr_n=1):
        """
        Generate a single cosmic ray event and include it to a cosmic ray map (self.cosmicrayMap).

        :param limit: limiting energy for the cosmic ray event
        :type limit: float
        :param cr_n: number of cosmic ray events to include
        :type cr_n: int

        :return: None
        """
        #pseudo-random numbers taken from a uniform distribution between 0 and 1
        luck = np.random.rand(cr_n)

        #draw the length of the tracks
        ius = InterpolatedUnivariateSpline(self.cr['cr_cdf'], self.cr['cr_u'])
        self.cr['cr_l'] = ius(luck)

        #set the energy directly to the limit
        self.cr['cr_e'] = np.asarray([limit, ]*cr_n)

        #Choose the properties such as positions and an angle from a random Uniform dist
        cr_x = self.xsize * np.random.rand(int(np.floor(cr_n)))
        cr_y = self.ysize * np.random.rand(int(np.floor(cr_n)))
        cr_phi = np.pi * np.random.rand(int(np.floor(cr_n)))

        #find the intercepts
        self.cosmicrayMap = self._cosmicRayIntercepts(self.cr['cr_e'], cr_x, cr_y, self.cr['cr_l'], cr_phi)

        #count the covering factor
        area_cr = np.count_nonzero(self.cosmicrayMap)
        text = 'The cosmic ray covering factor is %i pixels i.e. %.3f per cent' \
               % (area_cr, 100.*area_cr / (self.xsize*self.ysize))
        self.log.info(text)
        print(text)


    def _drawEventsToCoveringFactor(self, coveringFraction=3.0, limit=1000, verbose=False):
        """
        Generate cosmic ray events up to a covering fraction and include it to a cosmic ray map (self.cosmicrayMap).

        :param coveringFraction: covering fraction of cosmic rya events in per cent of total number of pixels
        :type coveringFraction: float
        :param limit: limiting energy for the cosmic ray event [None = draw from distribution]
        :type limit: None or float
        :param verbose: print out information to stdout
        :type verbose: bool


        :return: None
        """
        self.cosmicrayMap = np.zeros((self.ysize, self.xsize))

        #how many events to draw at once, too large number leads to exceeding the covering fraction
        
        # AZZO
        
        cdf = self.cr['cr_cdf']
        ucr = self.cr['cr_u']
        aproxpdf = (cdf[1:]-cdf[0:-1])/(ucr[1:]-ucr[0:-1])
        averlength = (aproxpdf*ucr[1:]).sum()/aproxpdf.sum()
        
        #averlength = ((1.-(cdf[1:]+cdf[0:-1])/2.)*(ucr[1:]-ucr[0:-1])).sum() # AZZO: average length of events
        cr_tot_guess = coveringFraction/100.*self.xsize*self.ysize/averlength
        
        cr_n = max(int(cr_tot_guess * 0.05),1) # allocating for a max. 5% error in cover. fraction, aprox.
                                               # Notice that the minimum number of events will be one...
        
        # END AZZO
        
        #cr_n = max(int(295 * self.information['exptime'] / 565. * coveringFraction / 1.4),1)

        covering = 0.0
        
        lengths = []
        energies = []
        cr_tot = 0
        
        #cr_n = 7689
        
        coveringFraction = 0       
        while covering < coveringFraction:
        #while True:
            #pseudo-random numbers taken from a uniform distribution between 0 and 1
            luck = np.random.rand(cr_n)

            #draw the length of the tracks
            #ius = InterpolatedUnivariateSpline(self.cr['cr_cdf'], self.cr['cr_u'])
            ius = interp1d(self.cr['cr_cdf'], self.cr['cr_u'])
            self.cr['cr_l'] = ius(luck)
            
            if limit is None:
                ius = InterpolatedUnivariateSpline(self.cr['cr_cde'], self.cr['cr_v'])
                self.cr['cr_e'] = ius(luck)
            else:
                #set the energy directly to the limit
                self.cr['cr_e'] = np.asarray([limit,])

            lengths += self.cr['cr_l'].tolist()
            energies += self.cr['cr_e'].tolist()
            
            #Choose the properties such as positions and an angle from a random Uniform dist
            cr_x = self.xsize * np.random.rand(int(np.floor(cr_n)))
            cr_y = self.ysize * np.random.rand(int(np.floor(cr_n)))
            cr_phi = np.pi * np.random.rand(int(np.floor(cr_n)))

            #find the intercepts
            self.cosmicrayMap += self._cosmicRayIntercepts(self.cr['cr_e'], cr_x, cr_y, self.cr['cr_l'], cr_phi)

            #count the covering factor
            area_cr = np.count_nonzero(self.cosmicrayMap)
            covering = 100.*area_cr / (self.xsize*self.ysize)

            text = 'The cosmic ray covering factor is %i pixels i.e. %.3f per cent' % (area_cr, covering)
            self.log.info(text)
            #stop()
            cr_tot += cr_n
            
            if verbose:
                print(text)
    
    
    def _drawEventsToFluxTime(self, flux=2.6, limit=1000, verbose=False):
        """
        Generate cosmic ray events given flux [input] and exposure time, and include it to a cosmic ray map (self.cosmicrayMap).

        :param flux: flux of cosmic ray events per unit area of active detector and unit time [cm-2 s-1].
        :type coveringFraction: float
        :param limit: limiting energy for the cosmic ray event [None = draw from distribution]
        :type limit: None or float
        :param verbose: print out information to stdout
        :type verbose: bool


        :return: None
        """
        self.cosmicrayMap = np.zeros((self.ysize, self.xsize))
        
        pix_cm2 = (self.cr['pixel_pitch'] * 1.E-4)**2.
        
        area_active = self.ysize * self.xsize * pix_cm2
        cr_n = flux * area_active * self.cr['exptime']
        
        cr_n = int(np.ceil(cr_n))
        
        
        lengths = []
        energies = []
        
        luck = np.random.rand(cr_n)

        #draw the length of the tracks
        #ius = InterpolatedUnivariateSpline(self.cr['cr_cdf'], self.cr['cr_u'])
        ius = interp1d(self.cr['cr_cdf'], self.cr['cr_u'])
        self.cr['cr_l'] = ius(luck)
            
        if limit is None:
            ius = InterpolatedUnivariateSpline(self.cr['cr_cde'], self.cr['cr_v'])
            self.cr['cr_e'] = ius(luck)
        else:
            #set the energy directly to the limit
            self.cr['cr_e'] = np.asarray([limit,])

        lengths += self.cr['cr_l'].tolist()
        energies += self.cr['cr_e'].tolist()
            
        #Choose the properties such as positions and an angle from a random Uniform dist
        cr_x = self.xsize * np.random.rand(int(np.floor(cr_n)))
        cr_y = self.ysize * np.random.rand(int(np.floor(cr_n)))
        cr_phi = np.pi * np.random.rand(int(np.floor(cr_n)))

        #find the intercepts
        self.cosmicrayMap += self._cosmicRayIntercepts(self.cr['cr_e'], cr_x, cr_y, self.cr['cr_l'], cr_phi)

        #count the covering factor
        area_cr = np.count_nonzero(self.cosmicrayMap)
        covering = 100.*area_cr / (self.xsize*self.ysize)

        text = 'The cosmic ray covering factor is %i pixels i.e. %.3f per cent' % (area_cr, covering)
        self.log.info(text)
        #stop()
        
            
        if verbose:
            print(text)
    

    def addCosmicRays(self, limit=None):
        """
        Include cosmic rays to the image given.

        :return: image with cosmic rays
        :rtype: ndarray
        """
        self._drawCosmicRays(limit=limit)

        #paste cosmic rays
        self.image += self.cosmicrayMap

        return self.image


    def addSingleEvent(self, limit=None):
        """
        Include a single cosmic ray event to the image given.

        :return: image with cosmic rays
        :rtype: ndarray
        """
        self._drawSingleEvent(limit=limit)

        #paste cosmic rays
        self.image += self.cosmicrayMap

        return self.image


    def addUpToFraction(self, coveringFraction, limit=None, verbose=False):
        """
        Add cosmic ray events up to the covering Fraction.

        :param coveringFraction: covering fraction of cosmic rya events in per cent of total number of pixels
        :type coveringFraction: float
        :param limit: limiting energy for the cosmic ray event [None = draw from distribution]
        :type limit: None or float
        :param verbose: print out information to stdout
        :type verbose: bool

        :return: image with cosmic rays
        :rtype: ndarray
        """
        self._drawEventsToCoveringFactor(coveringFraction, limit=limit, verbose=verbose)

        #paste cosmic rays
        self.image += self.cosmicrayMap

        return self.image


    def addToFluxTime(self,flux,limit=None,verbose=False):
        """
        Add cosmic ray events given detected events flux [cm-2 s-1] and exposure time.

        :param flux: nr. of cosmic ray events per cm^2 [of active detector area] and second
        :type flux: float
        :param limit: limiting energy for the cosmic ray event [None = draw from distribution]
        :type limit: None or float
        :param verbose: print out information to stdout
        :type verbose: bool

        :return: image with cosmic rays
        :rtype: ndarray
        """
        self._drawEventsToFluxTime(flux,limit=limit,verbose=verbose)

        #paste cosmic rays
        self.image += self.cosmicrayMap

        return self.image



if __name__ == "__main__":
    from vison.support import logger as lg
    from scipy import ndimage
    import astropy.io.fits as fts
    
    #set up logger
    log = lg.setUpLogger('VISsim.log')
    
    #test section
    crImage = np.zeros((2066, 2048), dtype=np.float64)
    
    #cosmic ray instance
    cosmics = cosmicrays(log, crImage)
    
    #add cosmic rays up to the covering fraction
    CCD_cr = cosmics.addUpToFraction(1.4, limit=None, verbose=True)
    
    effected = np.count_nonzero(CCD_cr)
    print(effected, effected*100./(CCD_cr.shape[0]*CCD_cr.shape[1]))
    
    #save to FITS
    fts.writeto(CCD_cr, 'cosmicrayTest.fits')
    
    #smooth with a charge diffusion kernel
    smooth = ndimage.filters.gaussian_filter(CCD_cr, (0.32, 0.32))
    
    fts.writeto(smooth, 'cosmicrayTestSmoothed.fits')
    
    