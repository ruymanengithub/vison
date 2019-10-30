#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

TEST: FPA-FWD (WARM)

Created on Thu Sep 26 16:11:14 2019

@author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import os
import numpy as np
import sys
from collections import OrderedDict
import pandas as pd
import string as st

from vison.support import utils
from vison.image import bits
from vison.other import MOT_FFaux
from vison.fpatests import fpatask
from vison.datamodel import inputs
from vison.fpatests import FW_aux

# END IMPORT

#def _basic_onLE1(mode,kwargs_retrieve,kwargs_apply):
#    print(mode)
#    if mode == 'retrieve':
#        
#        queue = kwargs_retrieve['queue']
#        CCDID = kwargs_retrieve['CCDID']
#        LE1 = kwargs_retrieve['LE1']
#        print('processing %s' % CCDID)
#        kccdobj = LE1.get_ccdobj(CCDID)
#        
#        reply = dict(CCDID=CCDID,
#                     Quads=kccdobj.Quads)
#        for Q in kccdobj.Quads:
#            reply[Q] = kccdobj.get_stats(Q,sector='img',statkeys=['mean'])
#        
#        queue.put(reply)
#
#    elif mode == 'apply':
#        
#        Quads = kwargs_apply['Quads']
#        CCDID = kwargs_apply['CCDID']
#        holder = kwargs_apply['holder']
#        print('applying %s' % CCDID)
#        for Q in Quads:
#            holder[CCDID][Q] = reply[Q]

class FWD_WARM_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()

class FWD_WARM(fpatask.FpaTask):
    
    inputsclass = FWD_WARM_inputs
    
    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """
        
        self.subtasks = [('check', self.check_data), 
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis),
                         ('appendix', self.appendix)]
        
        super(FWD_WARM, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                        cleanafter=cleanafter)
        
        self.name = 'FWD_WARM'
        self.type = 'Simple'
        
        self.figdict = FW_aux.get_FWfigs()
        self.CDP_lib = FW_aux.get_CDP_lib()
        
        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')
    

    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(preprocessing=dict(
                                            offsetkwargs=dict(
                                            ignore_pover= True, 
                                            trimscan = [25, 5], 
                                            method = 'median',
                                            extension= -1, 
                                            scan = 'pre' 
                                                )))    
    
    def _get_RAMPslope(self,vQ):
        """ """
        
        rawy = vQ.data['y'].copy()
        rownum = np.arange(len(rawy))
        rowsat = np.where(rawy==2.**16-1)[0][0]
        ixsel = np.where((rawy<5.E4) & (rownum<rowsat))
        
        pol = np.polyfit(rownum[ixsel],rawy[ixsel],1)
        
        slope, intercept = pol[0], pol[1]
        
        return slope,intercept
        
    
    def _basic_onLE1(self,**kwargs):
        """ """
        
        CCDID = kwargs['CCDID']
        LE1 = kwargs['LE1']
        vstart = kwargs['vstart']
        vend = kwargs['vend']
        kccdobj = LE1.get_ccdobj(CCDID)
        debug = kwargs['debug']
        
        
        HERprofs= MOT_FFaux.extract_overscan_profiles(kccdobj, 
                                            [1.E3, 4.E4], 
                                            direction='serial')
        
        for Q in self.Quads:
            
            # HERprofile
            
            # save HERprof
            
            self.dd.products['HER'][CCDID][Q] = HERprofs[Q].copy()
            
            
            # save HER value
            
            ixjump = HERprofs['ixjump']            
            self.dd.products['HERval'][CCDID][Q] = HERprofs[Q]['y'][ixjump]            
            
            
            if not debug:
                
                # Stats on pre and over scan
                
                # PENDING
                
                
                # vertical profile: RAMP
                
                vQ = kccdobj.get_1Dprofile(Q=Q, orient='ver', area='img', stacker='median',
                                           vstart=vstart, vend=vend)
                        
                self.dd.products['profiles1D'][CCDID][Q] = vQ.data.copy()
            
                # extract and save slope+intercept of profile for comparison
            
                rampslopeQ, rampinterQ = self._get_RAMPslope(vQ)
            
                self.dd.products['RAMPfits'][CCDID][Q] = (rampslopeQ,rampinterQ)
            
            
                # RAMP; bit-histograms extraction
                        
                bitsmean = bits.get_histo_bits(kccdobj, Q, vstart=vstart, vend=vend)
                bitsbin = np.arange(0,16)+0.5
            
                # save bit histograms
            
                self.dd.products['BITS'][CCDID][Q] = dict(bins=bitsbin, H=bitsmean)
            
    
    def basic_analysis(self):
        """        
        To-Do:
            Extract average profiles along columns (image area)
            Extract HER profiles and metrics
            Extract bits histograms
        """
        
        if self.report is not None:
            self.report.add_Section(keyword='extract', Title='Extraction', level=0)
                
        iObs = 0
        vstart = 0
        vend = 2086
                
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        LE1fits = self.dd.mx['File_name'][iObs]
        

        fullLE1fits = os.path.join(self.inputs['datapath'], LE1fits)
        
        LE1 = self.load_LE1(fullLE1fits)
        
        
        # DISPLAY IMAGE
        
        FWImgDict = self.get_ImgDictfromLE1(LE1,doequalise=True)
        
        self.figdict['FW_img'][1]['data'] = FWImgDict
        self.figdict['FW_img'][1]['plotter'] = self.metacal.plot_ImgFPA
                
        if self.report is not None:            
            self.addFigures_ST(figkeys=['FW_img'],dobuilddata=False)
        
        prodkeys = ['profiles1D','RAMPfits',
                    'HER','HERval','BITS']
        
        for prodkey in prodkeys:
        
            self.dd.products[prodkey] = dict()
            
            for jY in range(1,LE1.fpamodel.NSLICES+1):
                for iX in range(1,LE1.fpamodel.NSLICES+1):                
                    CCDID = 'C_%i%i' % (jY,iX)
                    self.dd.products[prodkey][CCDID] = dict()
        
        
                
        Bkwargs = dict(vstart=vstart,
                      vend=vend,
                      debug=False)
                
        
        self.iterate_over_CCDs(LE1, FWD_WARM._basic_onLE1, **Bkwargs)
        
        if self.report is not None:
        
            for prodkey in prodkeys:
                self.report.add_Text('product: %s, all extracted!' % prodkey)
        
            
        # Matrices: HERvalues, RAMPslopes
        
        # HERvalues matrix
        
        
        her_tb_cdp = self.CDP_lib['HER_TB']
        hercdpdict = dict(
                TBkey = 'HER_TB',
                meta=dict(),
                CDP_header=CDP_header,
                header_title='FWD WARM: HER TABLE',
                CDP_KEY='HER_TB_CDP',
                caption='HER value at first pixel.',
                valformat='%.2e'
                )
        
        def _getHERval(self,Ckey,Q):
            return self.dd.products['HERval'][Ckey][Q]
        
        self.add_StandardQuadsTable(extractor=_getHERval, 
                                    cdp=her_tb_cdp, 
                                    cdpdict=hercdpdict)
        
        # RAMP Slopes
        
        
        rslope_tb_cdp = self.CDP_lib['RSLOPE_TB']
        rslopecdpdict = dict(
                TBkey = 'RSLOPE_TB',
                meta=dict(),
                CDP_header=CDP_header,
                header_title='FWD WARM: RAMP SLOPES TABLE',
                CDP_KEY='RSLOPE_TB_CDP',
                caption='RAMP Slope in ADU/pixel row',
                valformat='%.1f'
                )
        
        
        def _getRSlope(self,Ckey,Q):
            return self.dd.products['RAMPfits'][CCDID][Q][0]
        
        self.add_StandardQuadsTable(extractor=_getRSlope, 
                                    cdp=rslope_tb_cdp, 
                                    cdpdict=rslopecdpdict)
        
        
        
    def _get_Rslopes_MAP(self, inData, Ckey, Q):
        return inData[Ckey][Q][0]
        
        
    def meta_analysis(self):
        """ """
        
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)
        
        
        # Display FPA-Map of ramp profiles
        
        def _assignProf1D(profdict, Ckey):
            """ """
            Cprofs = profdict[Ckey]
            
            res = dict(x=dict(),y=dict())
            
            for Q in self.Quads:
                res['x'][Q] = Cprofs[Q]['x']
                res['y'][Q] = Cprofs[Q]['y']
                
            return res
        
        
        FWRampsDict = self.iter_overCCDs(self.dd.products['profiles1D'], _assignProf1D)
        FWRampsDict['labelkeys'] = self.Quads
        
        
        self.figdict['FW_RAMPS'][1]['data'] = FWRampsDict
        self.figdict['FW_RAMPS'][1]['meta']['plotter'] = self.metacal.plot_XYMAP
        
        if self.report is not None:            
            self.addFigures_ST(figkeys=['FW_RAMPS'],dobuilddata=False)
        
        
        # Display FPA-Map of ramp slopes
        
        #RslopesMap = self.get_FPAMAP(self.dd.products['RAMPfits'],
        #                             extractor=self._get_Rslopes_MAP)
        
        #self.plot_SimpleMAP(RslopesMap, kwargs=dict(
        #                suptitle='FPA\_WARM: RAMP SLOPES [ADU/ROW]',
        #                figname = ''
        #                ))
        
        # Display FPA-Map of ramp-slope differences with CALCAMP
        
        # Dispay HER profiles
        
        # Display FPA-Map of HER-values
        
        # Display FPA-Map of differences with CALCAMP HER-values
        
        # Display all bit-histogram maps together
        
        
    def appendix(self):
        """Adds Appendices to Report."""
        
        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level=0)
   
        # TABLE: reference values of OFFSETS
        
        # TABLE: reference values of SLOPES
        
        # TABLE: reference values of HER
        
        
                
        