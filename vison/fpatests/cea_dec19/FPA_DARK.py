# -*- coding: utf-8 -*-
"""

TEST: FPA-DARK (COLD)

Created on Thu Nov 14 16:27:00 2019

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

from vison.support import utils, files
from vison.image import bits
from vison.other import MOT_FFaux
from vison.fpatests import fpatask
from vison.datamodel import inputs
from vison.datamodel import cdp as cdpmod
from vison.fpatests.cea_dec19 import DK_aux
from vison.inject import lib as ilib
# END IMPORT


class DARK_inputs(inputs.Inputs):
    manifesto = inputs.CommonFpaTaskInputs.copy()
    manifesto.update(OrderedDict(sorted([
        ('exptime', ([float], 'Exposure Time')),
    ])))

class DARK(fpatask.FpaTask):

    inputsclass = DARK_inputs

    def __init__(self, inputs, log=None, drill=False, debug=False, cleanafter=False):
        """ """

        self.subtasks = [('check', self.check_data),
                         ('basic', self.basic_analysis),
                         ('meta', self.meta_analysis),
                         ('appendix', self.appendix)]

        super(DARK, self).__init__(inputs=inputs, log=log, drill=drill, debug=debug,
                                       cleanafter=cleanafter)

        self.name = 'DARK'
        self.type = 'Simple'


        self.figdict = DK_aux.get_DKfigs()
        self.CDP_lib = DK_aux.get_CDP_lib()

        self.inputs['subpaths'] = dict(figs='figs',
                                       products='products')

        if self.inputs['inCDPs'] is not None:
            if 'gain' in self.inputs['inCDPs']:
                allgains = files.cPickleRead(self.inputs['inCDPs']['gain'])

                self.incdps['GAIN'] = OrderedDict()
                
                for block in self.fpa.flight_blocks:
                    self.incdps['GAIN'][block] = allgains[block]['PTC01'].copy()



    def set_inpdefaults(self, **kwargs):
        self.inpdefaults = self.inputsclass(preprocessing=dict(
                                            offsetkwargs=dict(
                                                ignore_pover=True,
                                                trimscan=[25, 5],
                                                method='median',
                                                extension=-1,
                                                scan='pre'
                                            )))

    def load_references(self):
        """ """

        # Load Reference Offsets
        
        refoffkey, offreg = ('offsets_fwd_cold','ove')

        refOFF_incdp = cdpmod.Json_CDP()
        refOFF_incdp.loadhardcopy(filef=self.inputs['inCDPs']['references'][refoffkey])
        

        def _get_ref_OFFs_MAP(inData, Ckey, Q):
            return inData[Ckey][Q][offreg]
        
        
        RefOFFsMap = self.get_FPAMAP(refOFF_incdp.data.copy(),
                                        extractor=_get_ref_OFFs_MAP)

        self.dd.products['REF_OFFs'] = RefOFFsMap.copy()


    def _basic_onLE1(self, **kwargs):
        """ """

        CCDID = kwargs['CCDID']
        LE1 = kwargs['LE1']
        vstart = kwargs['vstart']
        vend = kwargs['vend']

        kccdobj = LE1.get_ccdobj(CCDID)
        # print(CCDID)
        debug = kwargs['debug']        

        trimscans = dict(pre=[25, 5],
                         img=[5, 5],
                         ove=[5, 5])

        if not debug:

            for Q in self.Quads:

                # produce average profile along rows
                hor1Dprof = kccdobj.get_1Dprofile(
                                Q=Q, orient='hor', area='all', stacker='mean', vstart=vstart, vend=vend)

                self.dd.products['profiles_H_1D'][CCDID][Q] = hor1Dprof


                for reg in ['pre','img','ove']:

                    # produce average profile along cols
                    _ver1Dprof = kccdobj.get_1Dprofile(
                                Q=Q, orient='ver', area=reg, stacker='mean', 
                                vstart=vstart, vend=vend)
            
                    self.dd.products['profiles_V_1D_%s' % reg][CCDID][Q] = _ver1Dprof


                for reg in ['pre','img','ove']:

                    stats = kccdobj.get_stats(
                                    Q,sector=reg,
                                    statkeys=[
                                        'median','std'],
                                    trimscan=trimscans[reg],
                                    ignore_pover=True,
                                    extension=-1,
                                    VSTART=vstart,
                                    VEND=vend)

                    self.dd.products['MED_%s' % reg.upper()][CCDID][Q] = stats[0]
                    self.dd.products['STD_%s' % reg.upper()][CCDID][Q] = stats[1]

                # DARK Measurement: PENDING
                
                dark_meas = []

                for reg in ['img','ove']:

                    dark_meas += kccdobj.get_stats(
                                        Q,sector=reg,
                                        statkeys=['mean'],
                                        trimscan=trimscans[reg],
                                        ignore_pover=True,
                                        extension=-1,
                                        VSTART=vstart,VEND=vend,
                                        clip=[3,3])

                dark_val = dark_meas[0]-dark_meas[1]

                edark_val = (self.dd.products['STD_IMG'][CCDID][Q]**2.+\
                            self.dd.products['STD_OVE'][CCDID][Q]**2.)**0.5


                ncols_ove = kccdobj.overscan-trimscans['ove'][1]-trimscans['ove'][0]
                edark_val /= np.sqrt((vend-vstart)*ncols_ove) # limited by size of overscan region

                # MISSING: convert to e-/hr/pix

                BLOCK, CCDk, _, _ = self.fpa.FPA_MAP[CCDID]

                gain = self.incdps['GAIN'][BLOCK][CCDk][Q][0]

                convfactor = gain * 3600. / self.inputs['exptime']

                dark_val_epixhr = dark_val * convfactor
                edark_val_epixhr = edark_val * convfactor

                #if CCDID == 'C_66' and Q =='E':
                #    stop()

                self.dd.products['DARK'][CCDID][Q] = (dark_val_epixhr,edark_val_epixhr)




    def basic_analysis(self):
        """
        To-Do:
            Extract STDs and MEDs in pre/img/ove

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

        prodkeys = []
        for reg in ['pre','img','ove']:
            prodkeys.append('profiles_V_1D_%s' % reg)

        prodkeys += ['profiles_H_1D', 
            'STD_PRE', 'STD_IMG','STD_OVE',
            'MED_PRE', 'MED_IMG', 'MED_OVE',
            'DARK', 'EDARK']

        for prodkey in prodkeys:

            self.dd.products[prodkey] = dict()

            for jY in range(1, self.NSLICES_FPA + 1):
                for iX in range(1, self.NCOLS_FPA + 1):
                    CCDID = 'C_%i%i' % (jY, iX)
                    self.dd.products[prodkey][CCDID] = dict()

        Bkwargs = dict(vstart=vstart,
                       vend=vend,
                       debug=False)

        self.iterate_over_CCDs_inLE1(LE1, DARK._basic_onLE1, **Bkwargs)

        if self.report is not None:
            for prodkey in prodkeys:
                nprodkey = prodkey.replace('_','\_')
                self.report.add_Text('product: %s, all extracted!' % nprodkey)



        # Subtract OFFSET for display

        def offset_subtractor(ccdobj, **kwargs):
            """ """
            for Q in ccdobj.Quads:
                ccdobj.sub_offset(Q, **kwargs)
            return ccdobj
        
        suboffkwargs = dict(method='row', scan='ove', trimscan=[5, 5],
                   ignore_pover=True, extension=-1)

        LE1.apply_function_to_ccds(offset_subtractor, **suboffkwargs)

        # DISPLAY IMAGE

        SubOffImgDict = self.get_ImgDictfromLE1(LE1, doequalise=True)

        self.figdict['DK_img'][1]['data'] = SubOffImgDict
        self.figdict['DK_img'][1]['meta']['plotter'] = self.metacal.plot_ImgFPA


        # Matrix of offsets (over-scan)

        off_tb_cdp = cdpmod.Tables_CDP()
        offcdpdict = dict(
            caption='Offsets in the over-scan region.',
            valformat='%.2f')

        def _getOFF_OVEval(self, Ckey, Q):
            return self.dd.products['MED_OVE'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getOFF_OVEval,
                                    cdp=off_tb_cdp,
                                    cdpdict=offcdpdict)


        # Matrix of RONs (over-scan)

        ron_tb_cdp = cdpmod.Tables_CDP()
        roncdpdict = dict(
            caption='RON in the over-scan region.',
            valformat='%.2f')

        def _getRON_OVEval(self, Ckey, Q):
            return self.dd.products['STD_OVE'][Ckey][Q]

        self.add_StandardQuadsTable(extractor=_getRON_OVEval,
                                    cdp=ron_tb_cdp,
                                    cdpdict=roncdpdict)


        if self.report is not None:
            self.addFigures_ST(figkeys=['DK_img'], dobuilddata=False)

    def _get_XYdict_VPROFS(self):
        """ """
        
        data = dict(pre=self.dd.products['profiles_V_1D_pre'],
                    img=self.dd.products['profiles_V_1D_img'],
                    ove=self.dd.products['profiles_V_1D_ove'])

        regions = ['pre','img','ove']

        XYdict_VPROFS = dict(x=dict(),
                            y=dict(),
                            labelkeys=regions)
        for reg in regions:
            XYdict_VPROFS['x'][reg] = []
            XYdict_VPROFS['y'][reg] = []


        def assigner(VPROFS, data, Ckey):

            for Q in self.Quads:

                cenY = []

                for reg in regions:
                
                    _x = data[reg][Ckey][Q].data['x'].copy()
                    ixorder = np.argsort(_x)
                    _x = _x[ixorder]-_x.min()
                    _y = data[reg][Ckey][Q].data['y'][ixorder].copy()
                    VPROFS['x'][reg].append(_x)
                    VPROFS['y'][reg].append(_y)
                    cenY.append(np.median(_y))
                    #VPROFS['x'][reg] = [_x]
                    #VPROFS['y'][reg] =[_y]

                cenY = np.mean(cenY)

                for reg in regions: # recentering profiles
                    VPROFS['y'][reg][-1] -= cenY

            return VPROFS

        XYdict_VPROFS = self.iter_overCCDs(data, assigner, RetDict=XYdict_VPROFS)

        for reg in regions:
            XYdict_VPROFS['x'][reg] = np.array(XYdict_VPROFS['x'][reg]).T.tolist()
            XYdict_VPROFS['y'][reg] = np.array(XYdict_VPROFS['y'][reg]).T.tolist()


        return XYdict_VPROFS


    def _get_Histos_DARKS(self):
        """ """


        HistosDict = dict(y=[])


        def assigner(Histos, data, Ckey):
            for Q in self.Quads:
                Histos['y'].append(data[Ckey][Q][0])
            return Histos

        HistosDict = self.iter_overCCDs(self.dd.products['DARK'], assigner, RetDict=HistosDict)

        minval = -1
        maxval = +1

        bin_edges = np.linspace(minval,maxval,20)
        
        HistosDict['x'] = bin_edges

        return HistosDict


    def meta_analysis(self):
        """ """
        
        if self.report is not None:
            self.report.add_Section(keyword='meta', Title='Meta-Analysis', level=0)


        # Matrix of DARK CURRENT LEVELS

        DK_UNCs = [] 
        DKs = []

        for JY in range(self.NCOLS_FPA):
            for iX in range(self.NSLICES_FPA):
                Ckey = 'C_%i%i' % (JY + 1, iX + 1)
                for Q in self.Quads:
                    DKs.append(self.dd.products['DARK'][Ckey][Q][0])
                    DK_UNCs.append(self.dd.products['DARK'][Ckey][Q][1])


        avgDK = np.median(DKs)
        avgDK_UNK = np.median(DK_UNCs)

        dk_tb_cdp = cdpmod.Tables_CDP()
        dkcdpdict = dict(
            caption='Dark Current estimates in e-/pix/hr. '+\
                'Median (Value, Unc) = %.1e $\pm$ %.1e e-/pix/hr.' % (avgDK, avgDK_UNK),
            valformat='%s')


        def _getDK_val(self, Ckey, Q):
            return '%.2e' % self.dd.products['DARK'][Ckey][Q][0]

        self.add_StandardQuadsTable(extractor=_getDK_val,
                                    cdp=dk_tb_cdp,
                                    cdpdict=dkcdpdict)


        # Histogram of dark current values

        DARKS_HISTOS = self._get_Histos_DARKS()

        self.figdict['DK_HISTO'][1]['data'] = DARKS_HISTOS
        self.figdict['DK_HISTO'][1]['meta']['plotter'] = self.metacal.plot_HISTO

        # show average profiles along columns (single plot): pre-/img-/ove-
        #       Is there stray-light?

        XYdict_VPROFS = self._get_XYdict_VPROFS()

        self.figdict['VPROFS'][1]['data'] = XYdict_VPROFS
        self.figdict['VPROFS'][1]['meta']['plotter'] = self.metacal.plot_XY

        # PLOTTING

        if self.report is not None:
            self.addFigures_ST(figkeys=['DK_HISTO', 'VPROFS'], dobuilddata=False)



    def appendix(self):
        """Adds Appendices to Report."""

        if self.report is not None:
            self.report.add_Section(keyword='appendix', Title='Appendix', level=0)

        # TABLE: reference values of OFFSETS
        
        def _getRefOffs(self, Ckey, Q):
            return self.dd.products['REF_OFFs'][Ckey][Q]
        
        cdpdictoff = dict(
            caption = 'Reference OFFSETS [ADU] (GRCALCAMP).',
            valformat = '%.1f')
        
        self.add_StandardQuadsTable(extractor=_getRefOffs,
                                    cdp=None,
                                    cdpdict=cdpdictoff)
        