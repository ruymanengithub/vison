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
                         ('meta', self.meta_analysis)]
        
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
        
        
        HERprofs= MOT_FFaux.extract_overscan_profiles(kccdobj, 
                                            [1.E3, 4.E4], 
                                            direction='serial')
        
        for Q in self.Quads:
            
            # vertical profile: RAMP
            
            #vQ = kccdobj.get_1Dprofile(Q=Q, orient='ver', area='img', stacker='median',
            #                         vstart=vstart, vend=vend)
            
            
            #self.dd.products['profiles1D'][CCDID][Q] = vQ.data.copy()
            
            # extract and save slope+intercept of profile for comparison
            
            #rampslopeQ, rampinterQ = self._get_RAMPslope(vQ)
            
            #self.dd.products['RAMPfits'][CCDID][Q] = (rampslopeQ,rampinterQ)
            
            # HERprofile
            
            # save HERprof
            
            self.dd.products['HER'][CCDID][Q] = HERprofs[Q].copy()
            
            
            # save HER value: PENDING
            
            ixjump = HERprofs['ixjump']            
            self.dd.products['HERval'][CCDID][Q] = HERprofs[Q]['y'][ixjump]
            
            
            # RAMP; bit-histograms extraction
            
            #bitsmean = bits.get_histo_bits(kccdobj, Q, vstart=vstart, vend=vend)
            #bitsbin = np.arange(0,16)+0.5
            
            # save bit histograms
            
            #self.dd.products['BITS'][CCDID][Q] = dict(bins=bitsbin, H=bitsmean)
            
    
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
        
        figspath = self.inputs['subpaths']['figs']
        
        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()
        
        LE1fits = self.dd.mx['File_name'][iObs]
        

        fullLE1fits = os.path.join(self.inputs['datapath'], LE1fits)
        
        LE1 = self.load_LE1(fullLE1fits)
        
        
        # DISPLAY IMAGE
        
        FWImgDict = self.get_ImgDictfromLE1(LE1,doequalise=True)
        FWImg_figname = os.path.join(figspath,self.figdict['FW_img'][1]['figname'])
        imgkwargs = dict(figname=FWImg_figname)
        self.plot_ImgFPA(FWImgDict, imgkwargs)
        
        
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
                      vend=vend)
        
        self.iterate_over_CCDs(LE1, FWD_WARM._basic_onLE1, **Bkwargs)
        
        if self.report is not None:
        
            for prodkey in prodkeys:
                self.report.add_Text('product: %s, all extracted!' % prodkey)
        
            
        # Matrices: HERvalues, RAMPslopes
        
        NQuads = len(self.Quads)
        NCCDs = len(self.CCDs)
        NBlocks = len(self.fpa.flight_blocks)
        NP = NBlocks * NCCDs * NQuads
        
        HER_TB = OrderedDict()
        HER_TB['CCDID'] = np.zeros(NP,dtype='S4')
        HER_TB['BLOCK'] = np.zeros(NP,dtype='S4')
        HER_TB['CCD'] = np.zeros(NP,dtype='S4')
        HER_TB['Q'] = np.zeros(NP,dtype='S1')
        for Q in self.Quads:
            HER_TB[Q] = np.zeros(NP,dtype='float32')
        
        
        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                
                Ckey  = 'C_%i%i' % (jY+1,iX+1)
                
                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]
                
                for iQ,Q in enumerate(self.Quads):
                    indx = (jY * self.NCOLS_FPA + iX)*NQuads + iQ
                
                    HER_TB['CCDID'][indx] = Ckey
                    HER_TB['BLOCK'][indx] = block[0:4]
                    HER_TB['CCD'][indx] = CCDk
                    HER_TB['Q'][indx] = Q
                    HER_TB[Q][indx] = self.dd.products['HERval'][CCDID][Q]
        
        HER_TB_dddf = OrderedDict(HER_TB = pd.DataFrame.from_dict(HER_TB))
        
        her_tb_cdp = self.CDP_lib['HER_TB']
        her_tb_cdp.path = self.inputs['subpaths']['products']
        her_tb_cdp.ingest_inputs(
                data = HER_TB_dddf.copy(),
                meta=dict(),
                header=CDP_header.copy()
                )

        her_tb_cdp.init_wb_and_fillAll(header_title='FWD WARM: HER TABLE')
        self.save_CDP(her_tb_cdp)
        self.pack_CDP_to_dd(her_tb_cdp, 'HER_TB_CDP')

        if self.report is not None:
            
            fstr = lambda x: '%s' % x            
            fE = lambda x: '%.2E' % x
            
            her_formatters=[fstr,fstr,fstr,fstr]
            for iQ,Q in enumerate(self.Quads):
                her_formatters.append(fE)
            
            caption = 'Hard Edge Response, first pixel relative drift.'
            nicecaption = st.replace(caption,'_','\_')
            Htex = her_tb_cdp.get_textable(sheet='HER_TB', caption=nicecaption,
                                               fitwidth=True,
                                               tiny=True,
                                               formatters=her_formatters)
            
            
            self.report.add_Text(Htex)        
        
        
        
    def _get_Rslopes_MAP(self, inData, Ckey, Q):
        return inData[Ckey][Q][0]
        
        
    def meta_analysis(self):
        """ """
        
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)
        
        
        # Display FPA-Map of ramp profiles
        
        
        
        # Display FPA-Map of ramp slopes
        
        RslopesMap = self.get_FPAMAP(self.dd.products['RAMPfits'],
                                     extractor=self._get_Rslopes_MAP)
        
        self.plot_SimpleMAP(RslopesMap, kwargs=dict(
                        suptitle='FPA\_WARM: RAMP SLOPES [ADU/ROW]',
                        figname = ''
                        ))
        
        # Display FPA-Map of CALCAMP ramp slopes
        
        # Dispay HER profiles
        
        # Display FPA-Map of HER-values
        
        # Display FPA-Map of CALCAMP HER-values
        
        # Display all bit-histogram maps together
        
        
        
        
    
