# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: NL02

Similar to NL01, but using 2 wavelengths:
        - ND4 for low fluences
        - 880 nm for high fluences
        - Also possible to use with different values of RD (stability tests)
End-To-End Non-Linearity Curve

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check exposure time pattern matches test design.
    - Check quality of data (rough scaling of fluences with Exposure times).
    - Subtract offset level.
    - Divide by Flat-field.
    - Synoptic analysis:
        fluence ratios vs. extime ratios >> non-linearity curve
    - extract: Non-Linearity curve for each CCD and quadrant
    - produce synoptic figures
    - Save results.


Created on Tue Oct 23 15:22:00 2018

:author: raf

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import pandas as pd
import string as st

from vison.pipe.task import HKKeys
from vison.support import context
from vison.datamodel import cdp
from vison.datamodel import scriptic as sc
from vison.support import files
#from vison.pipe.task import Task
from .FlatTask import FlatTask
from vison.datamodel import inputs, core
from vison.support import utils
from vison.flat import NL01
from vison.flat import NL01aux
from . import nl as nllib
from . import nl_ptc

from pylab import plot, show
# END IMPORT

isthere = os.path.exists

# HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


NL02_commvalues = dict(program='CALCAMP', test='NL02',
                       flushes=7, siflsh=1, siflsh_p=500,
                       inisweep=1,
                       vstart=0, vend=2086,
                       exptime=0., shuttr=1, e_shuttr=0,
                       mirr_on=0,
                       motr_on=0,
                       source='flat',
                       comments='')


class NL02_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptimesA', ([dict, list], 'Exposure times for each fluence, PART A.')),
        ('framesA', ([list], 'Number of Frames for each fluence, PART A.')),
        ('wavelengthA', ([int], 'Wavelength, PART A')),
        ('exptimesB', ([dict, list], 'Exposure times for each fluence, PART B.')),
        ('framesB', ([list], 'Number of Frames for each fluence, PART B.')),
        ('wavelengthB', ([int], 'Wavelength, PART B')),
        ('exptinter', ([
            float], 'Exposure time for interleaved fluence stability control frames.'))
    ])))


class NL02(NL01.NL01):
    """ """

    inputsclass = NL02_inputs
    FLUDIVIDE = 20.  # pc

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        super(
            NL02,
            self).__init__(
            inputs=inputs,
            log=log,
            drill=drill,
            debug=debug,
            cleanafter=cleanafter)
        self.name = 'NL02'
        self.CDP_lib = NL01aux.get_CDP_lib(self.name)
        self.figdict = NL01aux.get_NL0Xfigs(self.name)

        self.subtasks = [('check', self.check_data), 
                         ('prep', self.prep_data),
                         ('extract', self.extract_stats),
                         ('simulNL', self.simulNL),
                         ('NL', self.produce_NLCs),
                         ('satCTE', self.do_satCTE),
                         ('extract_PTC', self.extract_PTC),
                         ('debugtask', self.debugtask)]
        
        self.pin = None # TESTS

    def set_inpdefaults(self, **kwargs):

        waveA = 0
        tFWCwA = self.ogse.profile['tFWC_flat']['nm%i' % waveA]

        ixLOWFLUlast = np.where(NL01.NL01_relfluences > self.FLUDIVIDE)[0][0]
        ixLOWFLU = (np.arange(ixLOWFLUlast),)
        exptsA = (NL01.NL01_relfluences[ixLOWFLU] / 100. *
                  tFWCwA).tolist()  # ms
        framesA = (np.ones(len(exptsA), dtype='int32') * 4).tolist()

        waveB = 880
        tFWCwB = self.ogse.profile['tFWC_flat']['nm%i' % waveB]
        ixHIFLUfirst = np.where(NL01.NL01_relfluences < self.FLUDIVIDE)[0][-1]
        #ixHIFLU = np.where(NL01.NL01_relfluences >= self.FLUDIVIDE)
        ixHIFLU = (np.arange(ixHIFLUfirst, len(NL01.NL01_relfluences)),)
        exptsB = (NL01.NL01_relfluences[ixHIFLU] / 100. *
                  tFWCwB).tolist()  # ms
        framesB = (np.ones(len(exptsB), dtype='int32') * 4).tolist()

        self.inpdefaults = dict(
            wavelengthA=waveA,
            exptimesA=exptsA,
            framesA=framesA,
            wavelengthB=waveB,
            exptimesB=exptsB,
            framesB=framesB,
            exptinter=0.5 * tFWCwB,
        )

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds NL02 script structure dictionary.

        """

        wavelengthA = self.inputs['wavelengthA']
        exptsA = self.inputs['exptimesA']
        framesA = self.inputs['framesA']

        wavelengthB = self.inputs['wavelengthB']
        exptsB = self.inputs['exptimesB']
        framesB = self.inputs['framesB']

        exptinter = self.inputs['exptinter']

        assert (len(exptsA) == len(framesA)) and \
               (len(exptsB) == len(framesB))

        FW_IDA = self.ogse.get_FW_ID(wavelengthA)
        FW_IDXA = int(FW_IDA[-1])

        FW_IDB = self.ogse.get_FW_ID(wavelengthB)
        FW_IDXB = int(FW_IDB[-1])

        #NL02_commvalues['wave'] = FW_IDX

        NL02_sdict = dict()

        NL02_sdict['col001'] = dict(frames=self.Nbgd, exptime=0, comments='BGD',
                                    wave=FW_IDXB)
        NL02_sdict['col002'] = dict(frames=self.Nstab0, exptime=exptinter, comments='STAB',
                                    wave=FW_IDXB)

        colcountbase = 3

        for ix, ifraA in enumerate(framesA):

            iexpA = exptsA[ix]

            colkeyFlu = 'col%03i' % (ix * 2 + colcountbase,)

            NL02_sdict[colkeyFlu] = dict(
                frames=ifraA, exptime=iexpA,
                wave=FW_IDXA,
                comments='Fluence%i' % (ix + 1,))

            colkeySta = 'col%03i' % (ix * 2 + colcountbase + 1,)

            NL02_sdict[colkeySta] = dict(
                frames=1, exptime=exptinter,
                wave=FW_IDXB,
                comments='STAB')

        colcountbase = colcountbase + 2 * len(framesA)

        for jx, jfraB in enumerate(framesB):

            jexpB = exptsB[jx]

            colkeyFlu = 'col%03i' % (jx * 2 + colcountbase,)

            NL02_sdict[colkeyFlu] = dict(
                frames=jfraB, exptime=jexpB,
                wave=FW_IDXB,
                comments='Fluence%i' % (jx + len(framesA) + 1,))

            colkeySta = 'col%03i' % (jx * 2 + colcountbase + 1,)

            NL02_sdict[colkeySta] = dict(
                frames=1, exptime=exptinter,
                wave=FW_IDXB,
                comments='STAB')

        Ncols = len(list(NL02_sdict.keys()))
        NL02_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(NL02_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except BaseException:
                diffvalues = diffvalues = dict()

        NL02_sdict = sc.update_structdict(NL02_sdict, commvalues, diffvalues)

        return NL02_sdict

    def prep_data(self):
        """

        Takes Raw Data and prepares it for further analysis.

        **METACODE**

        ::

            f.e. ObsID:
                f.e.CCD:
                    f.e.Q:
                        mask-out bad pixels
                        mask-out detector cosmetics
                        subtract offset
                        opt: [sub bias frame]

        """

        super(NL02, self).prepare_images(doExtract=True,
                                         doBadPixels=True,
                                         doMask=True,
                                         doOffset=False,
                                         doBias=False,
                                         doFF=False)


    def debugtask(self):

        from pylab import plot,show
        
        indices = copy.deepcopy(self.dd.indices)

        nObs, nCCD, nQuad = indices.shape[0:3]
        CCDs = indices.get_vals('CCD')
        Quads = indices.get_vals('Quad')

        binfactor = 10

        for jCCD, CCD in enumerate(CCDs):

            for kQ, Q in enumerate(Quads):

                jkoffset = self.dd.mx['offset_pre'][:, jCCD, kQ]

                medsbin = self.dd.mx['sec_med_bin{:d}'.format(binfactor)][:,jCCD,kQ,:].copy()
                medsbinsuboff = medsbin-np.expand_dims(jkoffset,-1)
                varsbin = self.dd.mx['sec_var_bin{:d}'.format(binfactor)][:,jCCD,kQ,:].copy()
                evarsbin = self.dd.mx['sec_evar_bin{:d}'.format(binfactor)][:,jCCD,kQ,:].copy()

                ixsel = np.where((medsbin>0) & (varsbin>0) & (evarsbin>0) & ~np.isnan(evarsbin))

                medsbinsuboff = medsbinsuboff[ixsel]
                varsbin = varsbin[ixsel]
                evarsbin = evarsbin[ixsel]

                indata = dict(mu_nle=medsbinsuboff,
                    var_nle = varsbin,
                    evar_nle=evarsbin,
                    ron=1.2,
                    gain=3.5,
                    binfactor=binfactor)

                NLres = nl_ptc.forward_PTC_LM(indata, npol=6)

                stop()

    def simulNL(self):
        """ """
        from matplotlib import pyplot as plt

        NLdeg = 4

        wpx = self.window['wpx']
        hpx = self.window['hpx']
        indices = copy.deepcopy(self.dd.indices)

        ishape = indices.shape

        nObs = ishape[0]

        Quads = indices.get_vals('Quad')
        CCDs = indices.get_vals('CCD')
        
        nC = len(CCDs)
        nQ = len(Quads)
        NP = nC * nQ


        tile_coos = dict()
        for Q in Quads:
            tile_coos[Q] = self.ccdcalc.get_tile_coos(Q, wpx, hpx, noedges=False)

        Nsectors = tile_coos[Quads[0]]['Nsamps']
        sectornames = np.arange(Nsectors)

        Sindices = copy.deepcopy(self.dd.indices)
        if 'Sector' in Sindices.names:
            Sindices.pop(Sindices.find('Sector'))
            self.dd.indices.pop(self.dd.indices.find('Sector'))
        Sindices.append(core.vIndex('Sector', vals=sectornames))

        valini = 0.
        
        self.dd.initColumn('sec_med_sim', Sindices, dtype='float32', valini=valini)
        self.dd.initColumn('sec_var_sim', Sindices, dtype='float32', valini=valini)

        iCCD = 0
        mexptimes = self.dd.mx['exptime'][:, iCCD].copy()
        wave = self.dd.mx['wave'][:, iCCD].copy()

        texptimes = self.recalibrate_exptimes(mexptimes)

        nonzeroexptimes = mexptimes[mexptimes > 0.]
        unonzeroexptimes = np.unique(nonzeroexptimes) # unique non-zero exposure times

        _ixstab = np.array([(mexptimes == iexptime).sum() for iexptime in unonzeroexptimes]).argmax()
        exptimestab = unonzeroexptimes[_ixstab]

        # boolean indices of the different types of exposures

        uwaves = np.unique(wave) # unique wavelengths used (2)

        fluences1E = self.dd.mx['sec_med'][:, 0, 0, :].copy()

        ixboo_bgd = mexptimes == 0.
        ixboo_stab = mexptimes == exptimestab
        ixboo_fluA = (mexptimes != 0.) & (mexptimes != exptimestab) & (wave == uwaves[0])
        ixboo_fluB = (mexptimes != 0.) & (mexptimes != exptimestab) & (wave == uwaves[1])

        fluxes = [np.nanmean((np.nanmean(fluences1E[ixboo_fluA,:],axis=1)/texptimes[ixboo_fluA])),
            np.nanmean((np.nanmean(fluences1E[ixboo_fluB,:],axis=1)/texptimes[ixboo_fluB]))]

        # identifying which is the high-flux wavelength / set, and which one is the low-flux

        if fluxes[0]>fluxes[1]:
            ixboo_fluHI = ixboo_fluA
            ixboo_fluLO = ixboo_fluB
            waveHI = uwaves[0]
            waveLO = uwaves[1]
            fluxHI = fluxes[0]
            fluxLO = fluxes[1]
        else:
            ixboo_fluLO = ixboo_fluA
            ixboo_fluHI = ixboo_fluB
            waveHI = uwaves[1]
            waveLO = uwaves[0]
            fluxHI = fluxes[1]
            fluxLO = fluxes[0]

        gain = 3.5
        pin = np.zeros(3+NLdeg+1)

        pin[0] = 2.
        pin[1] = 0.05
        pin[2] = .1
        #p[-5] = 0.
        pin[-4] = 0.6
        pin[-3] = -0.3
        pin[-2] = 0.5
        pin[-1] = 0.


        def f_non_Lin(x, p):
            """
            fNL_wExp(x, *p)
            return p[0] * np.exp(-(x - p[1]) / p[2]) + np.poly1d(p[3:])(x)"""
            xn = x / 2.**16

            return p[0] * np.exp(-(xn - p[1]) / p[2]) + np.poly1d(p[3:])(xn)

        fkx = np.linspace(100.,2.**16,100)
        fky = f_non_Lin(fkx, pin)

        doShowNL = False

        if doShowNL:

            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(fkx, fky, 'b-')
            ax.set_xlabel('fluence DN')
            ax.set_ylabel('NL pc')
            ax.set_title('Input Simulated NL curve')
            plt.show()


        for iCCD, CCDkey in enumerate(CCDs):

            for jQ, Q in enumerate(Quads):

                kk = iCCD * nQ + jQ
                
                raw_X = self.dd.mx['sec_X'][0, iCCD, jQ, :].copy()
                raw_Y = self.dd.mx['sec_Y'][0, iCCD, jQ, :].copy()

                spatDist = ((raw_X-raw_X.mean())**2.+(raw_Y-raw_Y.mean())**2.)**0.5

                #ijoffset = np.median(self.dd.mx['offset_pre'][:, iCCD, jQ])
                ijoffset = self.dd.mx['offset_pre'][:, iCCD, jQ]

                spatScale = (1.-spatDist/spatDist.max()*0.) # x% spatial drop in flux

                ijflu = np.zeros((nObs,Nsectors),dtype='float32')

                # ixboo_fluLO
                

                ijflu[ixboo_fluLO,:] = fluxLO * np.expand_dims(texptimes[ixboo_fluLO],1)
                ijflu[ixboo_fluLO,:] *= np.expand_dims(spatScale,0)

                ijflu[ixboo_fluHI,:] = fluxHI * np.expand_dims(texptimes[ixboo_fluHI],1)
                ijflu[ixboo_fluHI,:] *= np.expand_dims(spatScale,0)

                ijflu[ixboo_stab,:] = fluxHI * np.expand_dims(texptimes[ixboo_stab], 1)
                ijflu[ixboo_stab,:] *= np.expand_dims(spatScale,0)

                ijnlpc = f_non_Lin(ijflu, pin)

                ijvar = (np.sqrt(ijflu * gain)/gain)**2.
                ijvar += 1. # adding RON

                ijflu *= (1.+ijnlpc/100.) # adding non-linearity

                # adding noise
                noisescale = np.sqrt(ijvar) / np.sqrt(wpx*hpx)
                ijflu += np.random.normal(loc=0.0, scale=noisescale, size=ijflu.shape)
                # adding offset
                ijflu += np.expand_dims(ijoffset,1)


                self.dd.mx['sec_med_sim'][:, iCCD, jQ, :] = ijflu.copy()
                self.dd.mx['sec_var_sim'][:, iCCD, jQ, :] = ijvar.copy()

        self.pin = pin
        

    def produce_NLCs(self):
        """

        **METACODE**

        ::

            Obtains Best-Fit Non-Linearity Curve

            f.e. CCD:
                f.e. Q:

                    [opt] apply correction for source variability (interspersed exposure
                      with constant exptime)
                    Build NL Curve (NLC) - use stats and exptimes
                    fit poly. shape to NL curve

            plot NL curves for each CCD, Q
            report max. values of NL (table)

        """

        doExptimeCalib = True
        NLdeg = 4
        debug = False  # TESTS
        useSims = True # TESTS

        if self.report is not None:
            self.report.add_Section(
                keyword='NL', Title='Non-Linearity Analysis', level=0)
            if useSims:
                self.report.add_Text('Using SIMULATED DATA!')

        if doExptimeCalib:
            niceshutterprofname = self.ogse.profile['SHUTTER_CALIB'].replace('_', '\_')

            if self.report is not None:
                self.report.add_Text('Exposure times corrected using: %s' %
                                     niceshutterprofname)
            if self.log is not None:
                self.log.info('Exposure times corrected using: %s' %
                              niceshutterprofname)

        if self.report is not None:
            self.report.add_Text(['\nModel:',
                                  '\\begin{equation}',
                                  'NL=A \cdot e^{-(x-x_0)/s}} + P_%s(x)' % NLdeg,
                                  '\end{equation}'])


        dIndices = copy.deepcopy(self.dd.indices)

        CCDs = dIndices.get_vals('CCD')
        Quads = dIndices.get_vals('Quad')

        nC = len(CCDs)
        nQ = len(Quads)
        NP = nC * nQ

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        prodspath = self.inputs['subpaths']['products']

        # INITIALISATIONS

        # NON-LINEARITY (SUMMARY) TABLE

        NL_TB = OrderedDict()

        NL_TB['CCD'] = np.zeros(NP, dtype='int32')
        NL_TB['Q'] = np.zeros(NP, dtype='int32')
        NL_TB['MAXNLPC'] = np.zeros(NP, dtype='float32')
        NL_TB['FLU_MAXNLPC'] = np.zeros(NP, dtype='float32')

        # NON-LINEARITY (FULL) TABLE
        NL_FULL_TB = OrderedDict()

        NL_FULL_TB['FLUENCE'] = np.zeros(1000, dtype='float32') + np.nan
        for CCDk in CCDs:
            for Q in Quads:
                NL_FULL_TB['NLPC_%s%s' % (CCDk, Q)] = np.zeros(1000, dtype='float32') + np.nan

        # NON LINEARITY RESULTS

        NLall_mx = OrderedDict()

        for CCDkey in CCDs:
            NLall_mx[CCDkey] = OrderedDict()
            for Quad in Quads:
                NLall_mx[CCDkey][Quad] = OrderedDict()

        # NL CURVES

        curves_cdp = self.CDP_lib['NL_CURVES']
        curves_cdp.header = CDP_header.copy()
        curves_cdp.path = prodspath
        curves_cdp.data = OrderedDict()

        for CCDk in CCDs:
            curves_cdp.data[CCDk] = OrderedDict()
            for Q in Quads:
                curves_cdp.data[CCDk][Q] = OrderedDict()
                curves_cdp.data[CCDk][Q]['x'] = OrderedDict()
                curves_cdp.data[CCDk][Q]['y'] = OrderedDict()

        curves_cdp.data['labelkeys'] = ['data', 'fit']

        # Fitting the NL curves

        curves_res_plot = OrderedDict(labelkeys=[],
            x=dict(),
            y=dict())

        for iCCD, CCDkey in enumerate(CCDs):

            for jQ, Q in enumerate(Quads):

                ckey = '%s%s' % (CCDk,Q)

                kk = iCCD * nQ + jQ

                curves_res_plot['labelkeys'].append(ckey)

                if useSims:
                    raw_med = self.dd.mx['sec_med_sim'][:, iCCD, jQ, :].copy()
                    raw_var = self.dd.mx['sec_var_sim'][:, iCCD, jQ, :].copy()
                else:
                    raw_med = self.dd.mx['sec_med'][:, iCCD, jQ, :].copy()
                    raw_var = self.dd.mx['sec_var'][:, iCCD, jQ, :].copy()
                raw_X = self.dd.mx['sec_X'][:, iCCD, jQ, :].copy()
                raw_Y = self.dd.mx['sec_Y'][:, iCCD, jQ, :].copy()
                #col_labels = self.dd.mx['label'][:, iCCD].copy()
                exptimes = self.dd.mx['exptime'][:, iCCD].copy()
                wave = self.dd.mx['wave'][:, iCCD].copy()
                dtobjs = self.dd.mx['time'][:, iCCD].copy()
                ObsIDs = self.dd.mx['ObsID'][:].copy()

                #ijoffset = np.median(self.dd.mx['offset_pre'][:, iCCD, jQ])
                ijoffset = self.dd.mx['offset_pre'][:, iCCD, jQ]

                if doExptimeCalib:
                    nexptimes = self.recalibrate_exptimes(exptimes)
                else:
                    nexptimes = copy.deepcopy(exptimes)
                

                #nexptimes[nexptimes>0] -= .2 # ADHOC TESTS

                # fitresults = OrderedDict(coeffs, NLdeg, maxNLpc,flu_maxNLpc, bgd)
                if debug:
                    print(('\n%s%s\n' % (CCDkey, Q)))
                #print('WITH shutter nl correction...')
                
                _fitresults = nllib.wrap_fitNL_TwoFilters_Tests(raw_med, 
                        raw_var, 
                        nexptimes, 
                        wave,
                        dtobjs,
                        TrackFlux=True,
                        debug=debug and (iCCD==0) and (jQ==0),
                        ObsIDs=ObsIDs,
                        NLdeg=NLdeg,
                        # offset=0.)
                        offset=ijoffset,
                        XX=raw_X, YY=raw_Y,
                        pin=self.pin)
#                print('WITHOUT shutter nl correction...')
#                __fitresults = nllib.wrap_fitNL_TwoFilters_Alt(raw_med, raw_var, exptimes, wave,
#                                            dtobjs,
#                                            TrackFlux=True,
#                                            debug=debug,
#                                            ObsIDs=ObsIDs)

                NLall_mx[CCDkey][Q].update(_fitresults)

                NL_TB['CCD'][kk] = iCCD + 1
                NL_TB['Q'][kk] = jQ + 1
                NL_TB['MAXNLPC'][kk] = _fitresults['maxNLpc']
                NL_TB['FLU_MAXNLPC'][kk] = _fitresults['flu_maxNLpc']

                if kk == 0:
                    NL_FULL_TB['FLUENCE'] = _fitresults['outputcurve']['X'].copy()
                NL_FULL_TB['NLPC_%s%s' % (CCDkey, Q)] = _fitresults['outputcurve']['Y'].copy()

                curves_cdp.data[CCDkey][Q]['x']['data'] = _fitresults['inputcurve']['X'].copy() / \
                    1.E3
                curves_cdp.data[CCDkey][Q]['y']['data'] = _fitresults['inputcurve']['Y'].copy()
                curves_cdp.data[CCDkey][Q]['x']['fit'] = _fitresults['outputcurve']['X'].copy() / \
                    1.E3
                curves_cdp.data[CCDkey][Q]['y']['fit'] = _fitresults['outputcurve']['Y'].copy()


                # Residuals

                curves_res_plot['x'][ckey] = _fitresults['xres'].copy()
                curves_res_plot['y'][ckey] = _fitresults['yres'].copy()

        self.dd.products['NL'] = copy.deepcopy(NLall_mx)

        # Build Tables

        NL_TB_dddf = OrderedDict(NL_TB=pd.DataFrame.from_dict(NL_TB),
                                 NL_FULL_FIT=pd.DataFrame.from_dict(NL_FULL_TB))

        nl_tb_cdp = self.CDP_lib['NL_TB']
        nl_tb_cdp.path = prodspath
        nl_tb_cdp.ingest_inputs(
            data=NL_TB_dddf.copy(),
            meta=dict(),
            header=CDP_header.copy()
        )

        nl_tb_cdp.init_wb_and_fillAll(header_title='NL02: RESULTS TABLE')
        self.save_CDP(nl_tb_cdp)
        self.pack_CDP_to_dd(nl_tb_cdp, 'NL_TB_CDP')

        if self.report is not None:

            def fccd(x): return CCDs[x - 1]

            def fq(x): return Quads[x - 1]

            def ff(x): return '%.2f' % x

            formatters = [fccd, fq, ff, ff]

            caption = 'NL02 results TABLE'
            Ntex = nl_tb_cdp.get_textable(sheet='NL_TB',
                                          caption=caption,
                                          fitwidth=True,
                                          tiny=True,
                                          formatters=formatters)

            self.report.add_Text(Ntex)

        # Do plots

        self.save_CDP(curves_cdp)
        self.pack_CDP_to_dd(curves_cdp, 'NL_CURVES_CDP')

        fdict_NL = self.figdict['NL0X_fit_curves'][1]
        fdict_NL['data'] = curves_cdp.data.copy()
        #fdict_NL['caption'] = fdict_NL['caption'].replace('NL01', 'NL02')
        #fdict_NL['meta']['suptitle'] = fdict_NL['meta']['suptitle'].replace('NL01', 'NL02')
        #fdict_NL['figname'] = fdict_NL['figname'].replace( 'NL01', 'NL02')

        if self.report is not None:
            self.addFigures_ST(figkeys=['NL0X_fit_curves'],
                               dobuilddata=False)

        # All NL curves in a single plot, to see them better

        curves_single_plot = OrderedDict(labelkeys=[],
            x=dict(),
            y=dict())

        for CCDk in CCDs:
            for Q in Quads:
                ckey = '%s%s' % (CCDk,Q)
                curves_single_plot['labelkeys'].append(ckey)

                curves_single_plot['x'][ckey] = curves_cdp.data[CCDk][Q]['x']['fit'].copy()
                curves_single_plot['y'][ckey] = curves_cdp.data[CCDk][Q]['y']['fit'].copy()


        fdict_singNL = self.figdict['NL0X_fit_curves_single'][1]
        fdict_singNL['data'] = curves_single_plot.copy()

        if self.report is not None:
            self.addFigures_ST(figkeys=['NL0X_fit_curves_single'],
                               dobuilddata=False)        

        # Linearity residuals after applying NL curves


        fdict_resNL = self.figdict['NL0X_fit_curves_res'][1]
        fdict_resNL['data'] = curves_res_plot.copy()

        if self.report is not None:
            self.addFigures_ST(figkeys=['NL0X_fit_curves_res'],
                               dobuilddata=False)        

        self.canbecleaned = True

    def do_satCTE(self):
        """

        **METACODE**

        ::

            select ObsIDs with fluence(exptime) >~ 0.5 FWC

            f.e. ObsID:
                CCD:
                    Q:
                        measure CTE from amount of charge in over-scan relative to fluence

            f.e. CCD:
                Q:
                    get curve of CTE vs. fluence
                    measure FWC from curve in ADU

            report FWCs in electrons [via gain in inputs] f.e. CCD, Q (table)

        """
        raise NotImplementedError


    def extract_PTC(self):
        """ """

        from vison.flat.PTC0X import PTC0X

        # HARDWIRED VALUES
        wpx = self.window['wpx']
        hpx = self.window['hpx']

        if self.report is not None:
            self.report.add_Section(
                keyword='extract', Title='(binned) PTC Extraction', level=0)
            self.report.add_Text('Segmenting on %i x %i windows...' % (wpx, hpx))


        binfactor=10
        medcol = 'sec_med_bin{:d}'.format(binfactor)
        varcol = 'sec_var_bin{:d}'.format(binfactor)
        evarcol = 'sec_evar_bin{:d}'.format(binfactor)
        ccdobjcol = 'ccdobj_name'

        nl_ptc.f_extract_PTC(self, ccdobjcol, medcol, varcol, evarcol,
            binfactor=binfactor)

