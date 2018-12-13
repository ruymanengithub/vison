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
from FlatTask import FlatTask
from vison.datamodel import inputs, core
from vison.support import utils
from vison.flat import NL01
from vison.flat import NL01aux
import nl as nllib
# END IMPORT

isthere = os.path.exists

#HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
#          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
#          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
#          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
#          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']


NL02_commvalues = dict(program='CALCAMP',test='NL02',
                       flushes=7, siflsh=1, siflsh_p=500,
                       inisweep=1,
                       vstart=0, vend=2086,
                       exptime=0., shuttr=1,e_shuttr=0,
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
    FLUDIVIDE = 20. # pc

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        super(NL02, self).__init__(inputs, log, drill, debug)
        self.name = 'NL02'

    def set_inpdefaults(self, **kwargs):
        
        
        waveA = 0
        tFWCwA = self.ogse.profile['tFWC_flat']['nm%i' % waveA]
        
        ixLOWFLU = np.where(NL01.NL01_relfluences<self.FLUDIVIDE)
        exptsA = (NL01.NL01_relfluences[ixLOWFLU]/100. *
                 tFWCwA).tolist()  # ms
        framesA = (np.ones(len(exptsA), dtype='int32')*4).tolist()
        
        waveB = 880
        tFWCwB = self.ogse.profile['tFWC_flat']['nm%i' % waveB]
        ixHIFLU = np.where(NL01.NL01_relfluences>=self.FLUDIVIDE)
        exptsB = (NL01.NL01_relfluences[ixHIFLU]/100. *
                 tFWCwB).tolist()  # ms
        framesB = (np.ones(len(exptsB), dtype='int32')*4).tolist()
        
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

            colkeyFlu = 'col%03i' % (ix*2+colcountbase,)

            NL02_sdict[colkeyFlu] = dict(
                frames=ifraA, exptime=iexpA, 
                wave=FW_IDXA,
                comments='Fluence%i' % (ix+1,))

            colkeySta = 'col%03i' % (ix*2+colcountbase+1,)

            NL02_sdict[colkeySta] = dict(
                frames=1, exptime=exptinter, 
                wave=FW_IDXB,
                comments='STAB')
        
        colcountbase = colcountbase + 2*len(framesA)
        
        for jx, jfraB in enumerate(framesB):

            jexpB = exptsB[jx]

            colkeyFlu = 'col%03i' % (jx*2+colcountbase,)

            NL02_sdict[colkeyFlu] = dict(
                frames=jfraB, exptime=jexpB, 
                wave=FW_IDXB,
                comments='Fluence%i' % (jx+len(framesA)+1,))

            colkeySta = 'col%03i' % (jx*2+colcountbase+1,)

            NL02_sdict[colkeySta] = dict(
                frames=1, exptime=exptinter, 
                wave=FW_IDXB,
                comments='STAB')

        Ncols = len(NL02_sdict.keys())
        NL02_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(NL02_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        NL02_sdict = sc.update_structdict(NL02_sdict, commvalues, diffvalues)

        return NL02_sdict



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

        if self.report is not None:
            self.report.add_Section(
                keyword='NL', Title='Non-Linearity Analysis', level=0)
        
        debug=False # TESTS

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
        
        
        # NON-LINEARITY TABLE
        
        NL_TB = OrderedDict()
        
        NL_TB['CCD'] = np.zeros(NP,dtype='int32')
        NL_TB['Q'] = np.zeros(NP,dtype='int32')
        NL_TB['MAXNLPC'] = np.zeros(NP,dtype='float32')
        NL_TB['FLU_MAXNLPC'] = np.zeros(NP,dtype='float32')
        
        
        # NON LINEARITY RESULTS

        NLall_mx = OrderedDict()

        for CCDkey in CCDs:
            NLall_mx[CCDkey] = OrderedDict()
            for Quad in Quads:
                NLall_mx[CCDkey][Quad] = OrderedDict()
                
        # NL CURVES
        
        curves_cdp = cdp.CDP()
        curves_cdp.header = CDP_header.copy()
        curves_cdp.path = prodspath
        curves_cdp.data = OrderedDict()
        
        for CCDk in CCDs:
            curves_cdp.data[CCDk] = OrderedDict()            
            for Q in Quads:
                curves_cdp.data[CCDk][Q] = OrderedDict()
                curves_cdp.data[CCDk][Q]['x'] = OrderedDict()
                curves_cdp.data[CCDk][Q]['y'] = OrderedDict()
        
        curves_cdp.data['labelkeys'] = ['data','fit']
        
        # Fitting the NL curves
        
        for iCCD, CCDkey in enumerate(CCDs):
            
            for jQ, Q in enumerate(Quads):
                
                kk = iCCD * nQ + jQ

                raw_med = self.dd.mx['sec_med'][:, iCCD, jQ, :].copy()
                raw_var = self.dd.mx['sec_var'][:,iCCD,jQ, :].copy()
                #col_labels = self.dd.mx['label'][:, iCCD].copy()
                exptimes = self.dd.mx['exptime'][:, iCCD].copy()
                wave = self.dd.mx['wave'][:, iCCD].copy()
                dtobjs = self.dd.mx['time'][:, iCCD].copy()
                
                
                # fitresults = OrderedDict(coeffs, NLdeg, maxNLpc,flu_maxNLpc, bgd)
                if debug:
                    print('\n%s%s\n' % (CCDkey,Q))
                _fitresults = nllib.wrap_fitNL_TwoFilters(raw_med, raw_var, exptimes, wave, 
                                            dtobjs, 
                                            TrackFlux=True,
                                            subBgd=True, 
                                            debug=debug) 
                
                NLall_mx[CCDkey][Q].update(_fitresults)
                
                NL_TB['CCD'][kk] = iCCD+1
                NL_TB['Q'][kk] = jQ+1
                NL_TB['MAXNLPC'][kk] = _fitresults['maxNLpc']
                NL_TB['FLU_MAXNLPC'][kk] = _fitresults['flu_maxNLpc']
                
                curves_cdp.data[CCDkey][Q]['x']['data'] = _fitresults['inputcurve']['X'].copy()
                curves_cdp.data[CCDkey][Q]['y']['data'] = _fitresults['inputcurve']['Y'].copy()
                curves_cdp.data[CCDkey][Q]['x']['fit'] = _fitresults['outputcurve']['X'].copy()
                curves_cdp.data[CCDkey][Q]['y']['fit'] = _fitresults['outputcurve']['Y'].copy()
        
        self.dd.products['NL'] = copy.deepcopy(NLall_mx)

        # Build Tables
        
        NL_TB_dddf = OrderedDict(NL_TB = pd.DataFrame.from_dict(NL_TB))
        
        nl_tb_cdp = self.CDP_lib['NL_TB']
        nl_tb_cdp.path = prodspath
        nl_tb_cdp.ingest_inputs(
                data = NL_TB_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )

        nl_tb_cdp.init_wb_and_fillAll(header_title='NL02: RESULTS TABLE')
        self.save_CDP(nl_tb_cdp)
        self.pack_CDP_to_dd(nl_tb_cdp, 'NL_TB_CDP')
        
        if self.report is not None:
            
            fccd = lambda x: CCDs[x-1]
            fq = lambda x: Quads[x-1]
            ff = lambda x: '%.2f' % x
            
            formatters=[fccd,fq,ff,ff]
            
            caption = 'NL02 results TABLE' 
            Ntex = nl_tb_cdp.get_textable(sheet='NL_TB', caption=caption,
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=formatters)
            
            
            self.report.add_Text(Ntex)        

        # Do plots

        fdict_NL = self.figdict['NL01_fit_curves'][1]
        fdict_NL['data'] = curves_cdp.data.copy()
        fdict_NL['caption'] = st.replace(fdict_NL['caption'],'NL01','NL02')
        fdict_NL['meta']['suptitle'] = st.replace(fdict_NL['meta']['suptitle'],'NL01','NL02')
        fdict_NL['figname'] = st.replace(fdict_NL['figname'],'NL01','NL02')
        
        if self.report is not None:
            self.addFigures_ST(figkeys=['NL01_fit_curves'], 
                               dobuilddata=False)
        
        

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
