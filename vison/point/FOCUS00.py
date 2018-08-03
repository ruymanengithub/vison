# -*- coding: utf-8 -*-
"""

VIS Ground Calibration
TEST: FOCUS00

Focus analysis script

Tasks:

    - Select exposures, get file names, get metadata (commandig, HK).
    - Check quality of data (integrated fluxes are roughly constant, matching expected level).
    - Subtract offset level.
    - Divide by Flat-field.
    - Crop stamps of the sources on each CCD/Quadrant.
       - save snapshot figures of sources.
    - for each source (5 x Nquadrants):
       - measure shape using Gaussian Fit
    - Find position of mirror that minimizes PSF sizes
    - Produce synoptic figures:
        source size and ellipticity across combined FOV (of 3 CCDs)
    - Save results.

Created on Mon Apr 03 16:21:00 2017

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import numpy as np
from pdb import set_trace as stop
import os
import copy
from collections import OrderedDict
import pandas as pd
from scipy import stats
from scipy import interpolate

#from vison.pipe import lib as pilib
from vison.support import context
from vison.flat import FlatFielding as FFing
from vison.support.report import Report
from vison.support import files
from vison.point import lib as polib
from vison.datamodel import ccd
from vison.datamodel import EXPLOGtools as ELtools
from vison.datamodel import HKtools
from vison.datamodel import ccd
from vison.datamodel import generator
from vison.datamodel import scriptic as sc
from vison.point.spot import Spot as SpotObj
from vison.point import display as pdspl
from vison.support import vistime
from vison.support import utils
#from vison.pipe.task import Task
import PointTask as PT
#from PointTask import PointTask, BGD_lims
import FOCUS00_lib as F00lib
from vison.image import performance
from vison.datamodel import inputs
import F00aux
# END IMPORT

isthere = os.path.exists

HKKeys = ['CCD1_OD_T', 'CCD2_OD_T', 'CCD3_OD_T', 'COMM_RD_T',
          'CCD2_IG1_T', 'CCD3_IG1_T', 'CCD1_TEMP_T', 'CCD2_TEMP_T', 'CCD3_TEMP_T',
          'CCD1_IG1_T', 'COMM_IG2_T', 'FPGA_PCB_TEMP_T', 'CCD1_OD_B',
          'CCD2_OD_B', 'CCD3_OD_B', 'COMM_RD_B', 'CCD2_IG1_B', 'CCD3_IG1_B', 'CCD1_TEMP_B',
          'CCD2_TEMP_B', 'CCD3_TEMP_B', 'CCD1_IG1_B', 'COMM_IG2_B']

dtobj_default = vistime.dtobj_default

stampw = 25

FOCUS00_commvalues = dict(program='CALCAMP', test='FOCUS_%i',
                          IPHI1=1, IPHI2=1, IPHI3=1, IPHI4=0,
                          rdmode='fwd_bas',
                          flushes=7,
                          shuttr=1,
                          vstart=0, vend=2086,
                          siflsh=0, siflsh_p=500,
                          mirr_on=1,
                          motr_on=0,
                          source='point')

FWHM_lims = OrderedDict(CCD1=OrderedDict(
    E=OrderedDict(ALPHA=[0.5, 10.])))  # CCD-Q-Spot, pixels
for Spotname in polib.Point_CooNom['names'][1:]:
    FWHM_lims['CCD1']['E'][Spotname] = copy.deepcopy(
        FWHM_lims['CCD1']['E']['ALPHA'])
for Q in ['F', 'G', 'H']:
    FWHM_lims['CCD1'][Q] = copy.deepcopy(FWHM_lims['CCD1']['E'])
for iCCD in [2, 3]:
    FWHM_lims['CCD%i' % iCCD] = copy.deepcopy(FWHM_lims['CCD1'])

Flu_lims = OrderedDict(CCD1=OrderedDict(E=OrderedDict(
    ALPHA=[2.*40./100.*2.**16, 2.*75./100.*2.**16])))  # CCD-Q-Spot
for Spotname in polib.Point_CooNom['names'][1:]:
    Flu_lims['CCD1']['E'][Spotname] = copy.deepcopy(
        Flu_lims['CCD1']['E']['ALPHA'])
for Q in ['F', 'G', 'H']:
    Flu_lims['CCD1'][Q] = copy.deepcopy(Flu_lims['CCD1']['E'])
for iCCD in [2, 3]:
    Flu_lims['CCD%i' % iCCD] = copy.deepcopy(Flu_lims['CCD1'])


class FOCUS00_inputs(inputs.Inputs):
    manifesto = inputs.CommonTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptime', ([float], 'Exposure time.')),
        ('wavelength', ([int], 'Wavelength')),
        ('deltafocus', ([float], 'delta-focus step, mm'))
    ])))


class FOCUS00(PT.PointTask):
    """ """

    inputsclass = FOCUS00_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False):
        """ """
        self.subtasks = [('check', self.check_data), 
                         #('prep', self.prep_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis)]
        super(FOCUS00, self).__init__(inputs, log, drill, debug)
        self.name = 'FOCUS00'
        self.type = 'Simple'
        
        self.HKKeys = HKKeys
        self.figdict = F00aux.gt_F00figs(self.inputs['wavelength'])
        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')
        

    def set_inpdefaults(self, **kwargs):

        tFWC800 = self.ogse.profile['tFWC_point']['nm%i' % 800]
        self.inpdefaults = dict(wavelength=800,
                                exptime=60./100.*tFWC800,
                                deltafocus=0.1)

    def set_perfdefaults(self, **kwargs):
        super(FOCUS00, self).set_perfdefaults(**kwargs)
        self.perfdefaults['BGD_lims'] = copy.deepcopy(PT.BGD_lims)
        self.perfdefaults['FWHM_lims'] = copy.deepcopy(FWHM_lims)
        self.perfdefaults['Flu_lims'] = copy.deepcopy(Flu_lims)

    def build_scriptdict(self, diffvalues=dict(), elvis=context.elvis):
        """Builds FOCUS00 script structure dictionary.

        #:param wavelength: int, [nm], wavelength.
        #:param exptime: int, [ms], exposure time.
        :param diffvalues: dict, opt, differential values.


        """

        wavelength = self.inputs['wavelength']
        exptime = self.inputs['exptime']
        delta_focus = self.inputs['deltafocus']

        FW_ID = self.ogse.get_FW_ID(wavelength)
        FW_IDX = int(FW_ID[-1])
        mirror_nom = self.ogse.profile['mirror_nom'][FW_ID]

        # FOCUS00_sdict = dict(col1=dict(frames=5,wave=FW_IDX,exptime=0,
        #                               mirr_pos=mirror_nom-5,
        #                               comments='BGD'))

        FOCUS00_sdict = dict()

        for i, j in enumerate(range(-3, 4, 1)):
            FOCUS00_sdict['col%i' % (i+1,)] = dict(frames=2,
                                                   test='FOCUS00_%i' % wavelength,
                                                   exptime=exptime,
                                                   mirr_pos=mirror_nom +
                                                   float(j)*delta_focus,
                                                   wave=FW_IDX,
                                                   comments='F%.1f' % float(j))

        Ncols = len(FOCUS00_sdict.keys())
        FOCUS00_sdict['Ncols'] = Ncols

        commvalues = copy.deepcopy(sc.script_dictionary[elvis]['defaults'])
        commvalues.update(FOCUS00_commvalues)

        if len(diffvalues) == 0:
            try:
                diffvalues = self.inputs['diffvalues']
            except:
                diffvalues = diffvalues = dict()

        FOCUS00_sdict = sc.update_structdict(
            FOCUS00_sdict, commvalues, diffvalues)

        return FOCUS00_sdict

    def filterexposures(self, structure, explog, OBSID_lims):
        """ """
        wavedkeys = ['motr_siz']
        return super(FOCUS00, self).filterexposures(structure, explog, OBSID_lims, colorblind=False,
                                                    wavedkeys=wavedkeys)
        # return pilib.filterexposures(structure,explogf,datapath,OBSID_lims,colorblind=False,
        #                      wavedkeys=wavedkeys,elvis=elvis)

    def prep_data(self):
        """ """
        raise NotImplementedError

    def basic_analysis(self):
        """ 
        This is just an assignation of values measured in check_data.        
        """

        CHKindices = copy.deepcopy(self.dd.indices)
        assert 'CCD' in CHKindices.names
        assert 'Quad' in CHKindices.names
        assert 'Spot' in CHKindices.names
        
        newkeys = ['x','y','x_ccd','y_ccd','fluence','fwhmx','fwhmy']
        for key in newkeys:
            chkkey = 'chk_%s' % key
            self.dd.addColumn(self.dd.mx[chkkey][:].copy(), key, CHKindices)

    def meta_analysis(self):
        """ """
        # fit fwhmZ vs mirror_pos across focal plane (CCDxQxSpot)> 
        #                          best_mirr_pos(x,y)
        #   repeat for Z=x,y
        # find & report mean(medians) of best_mirror_posX/Y across focal plane> 
        #                          Best_mirr_pos(beam)
        # plot fwhmZ vs. mirror_pos with fits and Best_mirr_pos(beam)
        # plot predicted (surface) delta-fwhmZ across focal plane for 
        #                                         Best_mirr_pos(beam)
        # save Best Mirr pos as a "product"
        
        wave = self.inputs['wavelength']
        
        poldegree = 2
        
        CCDs = self.dd.indices.get_vals('CCD')
        Quads = self.dd.indices.get_vals('Quad')
        Spots = self.dd.indices.get_vals('Spot')
        
        xmirr = self.dd.mx['mirr_pos'][:].copy()
        
        nC = len(CCDs)
        nQ = len(Quads)
        nS = len(Spots)
        NP = nC * nQ * nS
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        
        Focus_dd = OrderedDict()
        Focus_dd['CCD'] = np.zeros(NP,dtype='int32')
        Focus_dd['Q'] = np.zeros(NP,dtype='int32')
        Focus_dd['Spot'] = np.zeros(NP,dtype='int32')
        Focus_dd['x_ccd'] =  np.zeros(NP,dtype='float32')
        Focus_dd['y_ccd'] =  np.zeros(NP,dtype='float32')
        Focus_dd['focus'] =  np.zeros(NP,dtype='float32')
        Focus_dd['belocal_fwhm'] =  np.zeros(NP,dtype='float32')
        
        zres = OrderedDict()
        
        for iCCD, CCDk in enumerate(CCDs):
            
            zres[CCDk] = OrderedDict()
            
            for jQ, Q in enumerate(Quads):
                
                zres[CCDk][Q] = OrderedDict()
                
                for lSpot, SpotName in enumerate(Spots):

                    ix = iCCD  * (nQ*nS) + jQ * nS + lSpot
                    
                    fwhmx = self.dd.mx['fwhmx'][:,iCCD,jQ,lSpot].copy()
                    fwhmy = self.dd.mx['fwhmy'][:,iCCD,jQ,lSpot].copy()
                    
                    fwhm = np.sqrt(fwhmx**2.+fwhmy**2.) / np.sqrt(2.)
                    
                    try:
                        p_zres = F00lib.fit_focus_single(xmirr[:,iCCD], 
                                fwhm, yerror=None, 
                                 degree=poldegree, doplot=False)
                    except ValueError:
                        if self.log is not None:
                            self.log.info('Focus fit did not converge for %s-%s-%s' %
                                  (CCDk,Q,SpotName))
                        p_zres = dict(coeffs=np.zeros(poldegree+1)+np.nan,
                         ecoeffs=np.zeros(poldegree+1)+np.nan,
                         focus=np.nan)
                    
                    Focus_dd['CCD'][ix] = iCCD 
                    Focus_dd['Q'][ix] = jQ
                    Focus_dd['Spot'][ix] = lSpot
                    
                    av_x_ccd = np.nanmean(self.dd.mx['x_ccd'][:,iCCD,jQ,lSpot])
                    av_y_ccd = np.nanmean(self.dd.mx['y_ccd'][:,iCCD,jQ,lSpot])
                                
                    Focus_dd['x_ccd'][ix] = av_x_ccd
                    Focus_dd['y_ccd'][ix] = av_y_ccd
                    
                    Focus_dd['focus'][ix] = p_zres['focus']
                    
                    Focus_dd['belocal_fwhm'][ix] = np.polyval(p_zres['coeffs'],p_zres['focus'])
                            
                    zres[CCDk][Q][SpotName] = p_zres

        # Find centroid of focus values, as Best Focus
        
        cbe_focus = np.nanmean(stats.sigmaclip(Focus_dd['focus'],4,4).clipped)
        
        # Obtaining delta-focus map
        
        beglobal_fwhm = np.zeros_like(Focus_dd['belocal_fwhm'],dtype='float32')
        
        for iCCD, CCDk in enumerate(CCDs):
            
            for jQ, Q in enumerate(Quads):
                                
                for lSpot, SpotName in enumerate(Spots):

                    ix = iCCD  * (nQ*nS) + jQ * nS + lSpot
                    
                    pcoeffs = zres[CCDk][Q][SpotName]['coeffs']
                    beglobal_fwhm[ix] = np.polyval(pcoeffs,cbe_focus)
        
        cbe_fwhm = np.nanmean(stats.sigmaclip(beglobal_fwhm,4,4).clipped)
        
        delta_fwhm = beglobal_fwhm - Focus_dd['belocal_fwhm']
        x_ccd = Focus_dd['x_ccd'].copy()
        y_ccd = Focus_dd['y_ccd'].copy()
        
        
        delta_fwhm_dict = OrderedDict()
        
        ccalc = self.ccdcalc # ALIAS
        xlims = (ccalc.prescan, 
                         ccalc.prescan+ccalc.NcolsCCD-1.) 
        ylims = (0., ccalc.NrowsCCD-1.) 
        
        for iCCD, CCDk in enumerate(CCDs):
            
            delta_fwhm_dict[CCDk] = OrderedDict()
            
            for jQ, Q in enumerate(Quads):
                
                ixsel = (np.arange(nS)+iCCD*(nQ*nS)+jQ*nS,)
                
                Qsel = np.zeros_like(x_ccd[ixsel],dtype='S1')
                Qsel[:] = Q
                xQ, yQ = self.ccdcalc.cooconv_CCD_2_Qrel(x_ccd[ixsel],
                                                           y_ccd[ixsel],
                                                        Qsel)
                
                # TESTS        
                print 'FOCUS00.meta_lib running in TEST-HACKED MODE'
                xQ = np.array([100.,2000.,1000.,100.,2000.])
                yQ = np.array([2000.,2000.,1000.,100.,1000.])
                delta_fwhm_HK = xQ/xQ.max()+1.
                # END TESTS
                
                delta_fwhm_dict[CCDk][Q] = dict(
                        img = F00lib.build_fwhm_map_CQ(delta_fwhm_HK, # [ixsel],
                                                xQ, yQ,
                                                xlims, ylims))
        
        # Figure
        
        fdict = self.figdict['F00meta_deltafwhm'][1]
        fdict['data'] = delta_fwhm_dict.copy()
        if self.report is not None:
            self.addFigures_ST(figkeys=['F00meta_deltafwhm'], 
                               dobuilddata=False)
        
        
        # CDP with CBE FOCUS and related tables
        
        meta_contents = ['cbe_focus', 'cbe_fwhm', 'poldegree', 'CCDs', 'Quads', 'Spots']
        focus_meta = OrderedDict([(name, eval(name)) for name in meta_contents])
        
        FOCUS_dddf = OrderedDict(FOCUS=pd.DataFrame.from_dict(Focus_dd))
        
        focus_cdp = F00aux.CDP_lib['FOCUS']
        focus_cdp.rootname += '_%snm' % wave
        focus_cdp.path = self.inputs['subpaths']['products']
        focus_cdp.ingest_inputs(data = FOCUS_dddf.copy(),                             
                              meta=focus_meta.copy(),
                              header=CDP_header.copy())
        
        focus_cdp = self.init_and_fill_CDP(focus_cdp, 
            header_title='FOCUS00_%i: FOCUS' % wave)
        focus_cdp.init_wb_and_fillAll(header_title='FOCUS00_%i: FOCUS' % wave)
        self.save_CDP(focus_cdp)
        self.pack_CDP_to_dd(focus_cdp, 'FOCUS_CDP')

        if self.report is not None:
            FOCUStex = focus_cdp.get_textable(sheet='FOCUS', 
                caption='FOCUS00 @ %i nm' % wave)
            self.report.add_Text(FOCUStex)

