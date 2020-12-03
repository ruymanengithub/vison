#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 14:19:00 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
import copy
import numpy as np
from collections import OrderedDict
import string as st
import os
import matplotlib.cm as cm
from sklearn import linear_model

from vison.datamodel import cdp
from vison.fpa import fpa as fpamod
from vison.support import utils
from vison.metatests.metacal import MetaCal
from vison.plot import plots_fpa as plfpa
from vison.plot import baseplotclasses as bpc

from vison.support import vcal
from vison.datamodel import core as vcore
from vison.ogse import ogse
from vison.support import files
from vison.xtalk import xtalk as xtalkmod

from matplotlib import pyplot as plt
plt.switch_backend('TkAgg')
from matplotlib.colors import Normalize
# END IMPORT

cols2keep = [
    'test',
    'sn_ccd1',
    'sn_ccd2',
    'sn_ccd3',
    'sn_roe',
    'sn_rpsu',
    'exptime',
    'vstart',
    'vend',
    'rdmode',
    'flushes',
    'siflsh',
    'siflsh_p',
    'swellw',
    'swelldly',
    'inisweep',
    'cdpu_clk',
    'chinj',
    'chinj_on',
    'chinj_of',
    'id_wid',
    'id_dly',
    'chin_dly',
    'v_tpump',
    's_tpump',
    'v_tp_mod',
    's_tp_mod',
    'v_tp_cnt',
    's_tp_cnt',
    'dwell_v',
    'dwell_s',
    'toi_fl',
    'toi_tp',
    'toi_ro',
    'toi_ch',
    'motr',
    'motr_cnt',
    'motr_siz',
    'source',
    'wave',
    'mirr_on',
    'mirr_pos',
    'R1C1_TT',
    'R1C1_TB',
    'R1C2_TT',
    'R1C2_TB',
    'R1C3_TT',
    'R1C3_TB',
    'IDL',
    'IDH',
    'IG1_1_T',
    'IG1_2_T',
    'IG1_3_T',
    'IG1_1_B',
    'IG1_2_B',
    'IG1_3_B',
    'IG2_T',
    'IG2_B',
    'OD_1_T',
    'OD_2_T',
    'OD_3_T',
    'OD_1_B',
    'OD_2_B',
    'OD_3_B',
    'RD_T',
    'RD_B',
    'time',
    'HK_CCD1_TEMP_T',
    'HK_CCD2_TEMP_T',
    'HK_CCD3_TEMP_T',
    'HK_CCD1_TEMP_B',
    'HK_CCD2_TEMP_B',
    'HK_CCD3_TEMP_B',
    'HK_CCD1_OD_T',
    'HK_CCD2_OD_T',
    'HK_CCD3_OD_T',
    'HK_CCD1_OD_B',
    'HK_CCD2_OD_B',
    'HK_CCD3_OD_B',
    'HK_COMM_RD_T',
    'HK_COMM_RD_B',
    'HK_CCD1_IG1_T',
    'HK_CCD2_IG1_T',
    'HK_CCD3_IG1_T',
    'HK_CCD1_IG1_B',
    'HK_CCD2_IG1_B',
    'HK_CCD3_IG1_B',
    'HK_COMM_IG2_T',
    'HK_COMM_IG2_B',
    'HK_FPGA_BIAS_ID2',
    'HK_VID_PCB_TEMP_T',
    'HK_VID_PCB_TEMP_B',
    'HK_RPSU_TEMP1',
    'HK_FPGA_PCB_TEMP_T',
    'HK_FPGA_PCB_TEMP_B',
    'HK_RPSU_TEMP_2',
    'HK_RPSU_28V_PRI_I',
    'chk_NPIXOFF',
    'chk_NPIXSAT',
    'offset_pre',
    'offset_ove',
    'bgd_img',
    'std_pre',
    'std_ove',
    'chk_x',
    'chk_x_ccd',
    'chk_y',
    'chk_y_ccd',
    'chk_peak',
    'chk_fluence',
    'chk_fwhmx',
    'chk_fwhmy']


class PsfPlot(bpc.BasicPlot):
    def __init__(self, data, **kwargs):

        super(BeamPlot, self).__init__(**kwargs)

        meta = dict(suptitle='',
                    #ccdtitles=dict(CCD1='CCD1', CCD2='CCD2', CCD3='CCD3'),
                    doLegend=False)
        meta.update(kwargs)

        self.figsize = (12, 12)
        self.measkeys = ['MEAN', 'SLOPE']
        self.bfkeys = ['BFE', 'NOBFE']
        self.data = copy.deepcopy(data)
        self.meta = dict()
        self.meta.update(meta)
        self.handles = []
        self.labels = []
        self.fig = None
        self.axs = dict()
        self.axsarr = []

        self.corekwargs = dict()
        if 'corekwargs' in kwargs:
            self.corekwargs.update(kwargs['corekwargs'])

    def init_fig(self):
        self._init_fig_and_axes()

    def _init_fig_and_axes(self):
        """ """
        plt.close('all')
        fig, axsarr = plt.subplots(
            2, 2, sharex=True, sharey=True, figsize=self.figsize)
        self.fig = fig

        self.axsarr = axsarr

        # initialisation of self.axs

        k = 0

        for bfkey in self.bfkeys:
            self.axs[bfkey] = dict()
            for measkey in self.measkeys:
                self.axs[bfkey][measkey] = self.axsarr.flatten()[k]
                k+=1

    def _ax_core_funct(self, ax, MBdict, key=''):
        """ """
        """ """

        ckwargs = self.corekwargs.copy()

        if key != '':
            bins = MBdict['x'][key]
            h = MBdict['y'][key]

            label = key.replace('_', '\_')
            kwargs = dict(label=label, weights=None, cumulative=False,
                    histtype='step', align='mid', orientation='vertical', log=False)
            if key in ckwargs:
                kwargs.update(ckwargs[key])
            else:
                kwargs.update(ckwargs)
            #kwargs = dict()
            _, _, handle = ax.hist(h, bins=bins, **kwargs)
            
            ax.axvline(x=MBdict['mean'][key],color=kwargs['color'],ls='--')

        else:
            bins = MBdict['x']
            h = MBdict['y']
            kwargs = dict(weights=None, cumulative=False,
                    histtype='step', align='mid', orientation='vertical', log=False)
            kwargs.update(ckwargs)

            ax.hist(h, bins=bins, **kwargs)

            ax.axvline(x=MBdict['mean'],color=kwargs['color'],ls='--')

            handle, label = None, None
            

        return handle, label

    def populate_axes(self):
        """ """

        try:
            labelkeys = self.data['labelkeys']
        except KeyError:
            labelkeys = []

        k=0
        for bfkey in self.bfkeys:
            for measkey in self.measkeys:

                ax = self.axs[bfkey][measkey]
                BMdict = self.data[bfkey][measkey]

                if len(labelkeys) > 0:
                    for labelkey in labelkeys:
                        handle, label = self._ax_core_funct(ax, BMdict, labelkey)
                        if k==0:
                            self.handles += handle
                            self.labels.append(label)
                else:
                    _, _ = self._ax_core_funct(ax, CQdict)

                k += 1

                #if Q in ['E', 'H']:
                #    ax.text(0.05, 0.9, Q, horizontalalignment='left',
                #            transform=self.axs[CCDkey][Q].transAxes)
                #elif Q in ['F', 'G']:
                #    ax.text(0.9, 0.9, Q, horizontalalignment='right',
                #            transform=self.axs[CCDkey][Q].transAxes)

                #if Q == 'E':
                #    ax.set_title(CCDkey, x=1)

                #if 'xlabel' in self.meta and Q in ['H', 'G']:
                #    ax.set_xlabel(self.meta['xlabel'])
                #if 'ylabel' in self.meta and Q in ['E', 'H'] and CCD == 1:
                #    ax.set_ylabel(self.meta['ylabel'])

                #if 'ylim' in self.meta:
                #    ax.set_ylim(self.meta['ylim'])
                #if 'xlim' in self.meta:
                #    ax.set_xlim(self.meta['xlim'])

        # self.axs[CCDkey][Q].locator_params(nticks=4,axis='x')

    def plt_trimmer(self):

        for measkey in self.measkeys:
            plt.setp(self.axs[bfkey]['BFE'].get_xticklabels(), visible=False)

            if measkey == 'SLOPE':
                for bfkey in self.bfkeys:
                    plt.setp(self.axs[measkey][bfkey].get_yticklabels(), visible=False)

        if self.meta['doLegend']:
            plt.figlegend(self.handles, self.labels, loc='center right')

        plt.locator_params(axis='y', nbins=5, prune='both')

        try:
            plt.locator_params(axis='x', nbins=4, prune='both')
        except BaseException:
            pass

        plt.subplots_adjust(hspace=0.0)
        plt.subplots_adjust(wspace=0.0)

        plt.margins(0.05)

        plt.suptitle(self.meta['suptitle'])
        # plt.tight_layout()
        plt.subplots_adjust(top=0.85)
        if self.meta['doLegend']:
            plt.subplots_adjust(right=0.85)

        


class MetaPsf(MetaCal):
    """ """

    def __init__(self, *args, **kwargs):
        """ """

        super(MetaPsf, self).__init__(*args, **kwargs)

        self.Spots = ['ALPHA', 'BRAVO', 'CHARLIE', 'DELTA', 'ECHO']
        self.testnames = ['PSF01_590', 'PSF01_730', 'PSF01_800', 'PSF01_880']
        #self.testnames = ['PSF01_730', 'PSF01_800', 'PSF01_880'] # TESTS
        self.incols = cols2keep
        self.ParsedTable = OrderedDict()
        # self.blocks = self.blocks[1:] # TESTS!
        #self.blocks = ['BORN', 'CURIE', 'DIRAC', 'ERWIN', 'FOWLER', 'GUYE', 'JULES2', 'KRAMERS'\
        #'LORENTZ', 'NIELS', 'OWEN', 'EINSTEIN', 'SKLODOWSKA', 'MAX', 'HEISENBERG', \
        #'JULES'] # TESTS
        allgains = files.cPickleRead(kwargs['cdps']['gain'])

        self.cdps['GAIN'] = OrderedDict()
        for block in self.blocks:
            self.cdps['GAIN'][block] = allgains[block]['PTC01'].copy()

        self.products['XTALK'] = OrderedDict()
        self.products['XTALK_RT'] = files.cPickleRead(kwargs['cdps']['xtalk_roetab'])

        self.init_fignames()
        self.init_outcdpnames()

    def parse_single_test(self, jrep, block, testname, inventoryitem):
        """ """

        NCCDs = len(self.CCDs)
        NQuads = len(self.Quads)
        session = inventoryitem['session']

        CCDkeys = ['CCD%i' % CCD for CCD in self.CCDs]

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        IndexCQ = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
                                     vcore.vIndex('CCD',vals=self.CCDs),
                                     vcore.vIndex('Quad',vals=self.Quads)])

        # IndexCQS = vcore.vMultiIndex([vcore.vIndex('ix',vals=[0]),
        #                             vcore.vIndex('CCD',vals=self.CCDs),
        #                             vcore.vIndex('Quad',vals=self.Quads),
        #                             vcore.vIndex('Spot', vals=self.Spots)])

        #idd = copy.deepcopy(inventoryitem['dd'])
        sidd = self.parse_single_test_gen(jrep, block, testname, inventoryitem)

        # TEST SCPECIFIC
        # TO BE ADDED:
        #   BLOCK, TEST, REPEAT
        #   wavenm, calibrated HK (voltages),

        CHAMBER = sidd.meta['inputs']['CHAMBER']

        CHAMBER_key = CHAMBER[0]
        chamber_v = np.array([CHAMBER_key])
        sidd.addColumn(chamber_v, 'CHAMBERKEY', IndexS, ix=0)

        ogseobj = ogse.Ogse(CHAMBER=CHAMBER)

        wave = sidd.mx['wave'][0, 0]

        wave_v = np.array([ogseobj.get_wavelength(wave)])
        sidd.addColumn(wave_v, 'WAVENM', IndexS, ix=0)

        block_v = np.array([block])
        sidd.addColumn(block_v, 'BLOCK', IndexS, ix=0)

        test_v = np.array([jrep + 1])
        sidd.addColumn(test_v, 'REP', IndexS, ix=0)

        test_v = np.array([session])
        sidd.addColumn(test_v, 'SESSION', IndexS, ix=0)

        test_v = np.array([testname])
        sidd.addColumn(test_v, 'TEST', IndexS, ix=0)

        xtalkpath = os.path.join(inventoryitem['resroot'], 'xtalk')

        ctalkcdp_pick = os.path.join(xtalkpath, os.path.split(sidd.products['CTALK'])[-1])
        ctalkcdp = files.cPickleRead(ctalkcdp_pick)

        ctalkkey = 'XTALK_%s_%s_%s_%i' % (testname, block, session, jrep + 1)
        self.products['XTALK'][ctalkkey] = ctalkcdp['data'].copy()
        ctalkkey_v = np.array([ctalkkey])
        sidd.addColumn(ctalkkey_v, 'XTALK', IndexS, ix=0)

        # XTALK-from ROE-TALK


        # PSF fits w/wo B-F correction
        tmp_v_CQ = np.zeros((1, NCCDs, NQuads), dtype='float32')


        #i00, fwhmx, fwhmy, didfit
        fwhmx_bfe_mean_v = tmp_v_CQ.copy()
        fwhmx_bfe_slope_v = tmp_v_CQ.copy()

        fwhmy_bfe_mean_v = tmp_v_CQ.copy()
        fwhmy_bfe_slope_v = tmp_v_CQ.copy()

        fwhmx_nobfe_mean_v = tmp_v_CQ.copy()
        fwhmx_nobfe_slope_v = tmp_v_CQ.copy()

        fwhmy_nobfe_mean_v = tmp_v_CQ.copy()
        fwhmy_nobfe_slope_v = tmp_v_CQ.copy()

        dd = inventoryitem['dd']

        for iCCD, CCDk in enumerate(CCDkeys):

            for kQ, Q in enumerate(self.Quads):

                cXbfe = self._get_fwhm_flu_bfit(dd, iCCD, kQ, 'fwhmx', bfecorr=False)
                fwhmx_bfe_slope_v[0, iCCD, kQ] = cXbfe[0]
                fwhmx_bfe_mean_v[0, iCCD, kQ] = cXbfe[1]

                cYbfe = self._get_fwhm_flu_bfit(dd, iCCD, kQ, 'fwhmy', bfecorr=False)
                fwhmy_bfe_slope_v[0, iCCD, kQ] = cYbfe[0]
                fwhmy_bfe_mean_v[0, iCCD, kQ] = cYbfe[1]

                cXnobfe = self._get_fwhm_flu_bfit(dd, iCCD, kQ, 'fwhmx', bfecorr=True)
                fwhmx_nobfe_slope_v[0, iCCD, kQ] = cXnobfe[0]
                fwhmx_nobfe_mean_v[0, iCCD, kQ] = cXnobfe[1]

                cYnobfe = self._get_fwhm_flu_bfit(dd, iCCD, kQ, 'fwhmy', bfecorr=True)
                fwhmy_nobfe_slope_v[0, iCCD, kQ] = cYnobfe[0]
                fwhmy_nobfe_mean_v[0, iCCD, kQ] = cYnobfe[1]


        sidd.addColumn(fwhmx_bfe_mean_v, 'FWHMX_BFE_MEAN', IndexCQ)
        sidd.addColumn(fwhmx_bfe_slope_v, 'FWHMX_BFE_SLOPE', IndexCQ)
        sidd.addColumn(fwhmy_bfe_mean_v, 'FWHMY_BFE_MEAN', IndexCQ)
        sidd.addColumn(fwhmy_bfe_slope_v, 'FWHMY_BFE_SLOPE', IndexCQ)

        sidd.addColumn(fwhmx_nobfe_mean_v, 'FWHMX_NOBFE_MEAN', IndexCQ)
        sidd.addColumn(fwhmx_nobfe_slope_v, 'FWHMX_NOBFE_SLOPE', IndexCQ)
        sidd.addColumn(fwhmy_nobfe_mean_v, 'FWHMY_NOBFE_MEAN', IndexCQ)
        sidd.addColumn(fwhmy_nobfe_slope_v, 'FWHMY_NOBFE_SLOPE', IndexCQ)



        # flatten sidd to table

        sit = sidd.flattentoTable()

        return sit


    def _get_fwhm_flu_bfit(self, dd, iCCD, kQ, fwhmkey, bfecorr):
        """ """

        xkey = 'gau_i00'
        ykey = 'gau_%s' % fwhmkey
        fitkey = 'gau_didfit'
        if bfecorr:
            xkey = 'nobfe_%s' % xkey
            ykey = 'nobfe_%s' % ykey
            fitkey = 'nobfe_%s' % fitkey

        xdata = dd.mx[xkey][:,iCCD,kQ,:].copy()
        ydata = dd.mx[ykey][:,iCCD,kQ,:].copy()
        didfit = dd.mx[fitkey][:,iCCD,kQ,:].copy()

        nSpots = xdata.shape[1]

        slopes = []
        intercepts = []

        for i in range(nSpots):
            _x = xdata[:,i]
            _y = ydata[:,i]
            _df = didfit[:,i]

            ixsel = np.where((_x>1.e3) & (_df==1))

            if len(ixsel[0])>5:
                ransac = linear_model.RANSACRegressor()

                x2fit = (_x[ixsel]-2.**15)/1.e4

                x2fit = np.expand_dims(x2fit, 1)

                ransac.fit(x2fit, np.expand_dims(_y[ixsel], 1))

                slopes.append(ransac.estimator_.coef_[0][0])
                intercepts.append(ransac.estimator_.intercept_[0])

        slope = np.nanmean(slopes)
        intercept = np.nanmean(intercepts)
        coeffs = [slope,intercept]

        #p = np.poly1d(coeffs)
        #xfit = (np.linspace(1.,2.**16,5)-2.**15)/1.e4
        #ybest = np.polyval(p,xfit)

        #xbest = xfit+2.**15/1.e4

        #return xbest, ybest, coeffs
        return coeffs

    def get_XTALKDICT_from_PT(self, testname):
        """ """

        XTALKdict = OrderedDict()
        XTALKdict['blocks'] = self.flight_blocks

        PT = self.ParsedTable[testname]

        for block in self.flight_blocks:
            ixblock = np.where(PT['BLOCK'] == block)[0][0]
            ctalkkey = PT['XTALK'][ixblock]
            XTALKdict[block] = self.products['XTALK'][ctalkkey].copy()

        return XTALKdict

    def _get_xtalk(self, xtalk_dict, mode='sign'):

        coefs = xtalk_dict['coefs']
        xsource = np.linspace(0, 2.**16, 200)
        yvictim = xtalkmod.f_fitXT(xsource, *coefs)

        ixmax = np.argmax(np.abs(yvictim))
        ixtalk = yvictim[ixmax]
        if mode == 'sign':
            return ixtalk
        elif mode == 'abs':
            return np.abs(ixtalk)

    def _get_XYdict_XT(self, TEST1, TEST2, mode='sign'):
        """ """

        x = dict()
        y = dict()

        PT1 = self.ParsedTable[TEST1]

        if TEST2 != 'RT':
            PT2 = self.ParsedTable[TEST2]

        for block in self.flight_blocks:

            _x = []
            _y = []

            ixblock1 = np.where(PT1['BLOCK'] == block)[0][0]
            ctalkkey1 = PT1['XTALK'][ixblock1]
            ctalk1 = self.products['XTALK'][ctalkkey1].copy()

            if TEST2 == 'RT':
                ctalk2 = self.products['XTALK_RT'][block].copy()
            else:
                ixblock2 = np.where(PT2['BLOCK'] == block)[0][0]
                ctalkkey2 = PT2['XTALK'][ixblock2]
                ctalk2 = self.products['XTALK'][ctalkkey2].copy()

            for iCr, CCDref in enumerate(self.CCDs):
                for iQr, Qref in enumerate(self.Quads):

                    for iC, CCD in enumerate(self.CCDs):
                        for iQ, Q in enumerate(self.Quads):
                            CCDreftag = 'CCD%i' % CCDref
                            CCDtag = 'CCD%i' % CCD
                            try:
                                _x.append(self._get_xtalk(ctalk1[CCDreftag][Qref][CCDtag][Q], mode))
                                _y.append(self._get_xtalk(ctalk2[CCDreftag][Qref][CCDtag][Q], mode))
                            except KeyError:
                                pass
            x[block] = np.array(_x)
            y[block] = np.array(_y)

        x['oneone'] = [-100, 100]
        y['oneone'] = [-100, 100]

        labelkeys = self.flight_blocks + ['oneone']

        XTdict = dict(x=x, y=y, labelkeys=labelkeys)

        return XTdict


    def _get_Hdict_FWHM(self, testname):
        """ """
        # secondary paramters: fwhmkey, bfkey, measkey

        measkeys = ['MEAN','SLOPE']
        bfkeys = ['BFE','NOBFE']

        Hdict = dict()

        for bfkey in bfkeys:
            Hdict[bfkey] = dict()
            for measkey in measkeys:
                Hdict[bfkey][measkey] = dict(x=dict(),y=dict())

        Hdict['labelkeys'] = ['fwhmx','fwhmy']
            

        PT = self.ParsedTable[testname]

        for measkey in measkeys:
            for bfkey in bfkeys:
        
                _valy = []
                _valx = []

                for block in self.flight_blocks:

                    ixblock = np.where(PT['BLOCK'] == block)[0][0]
                    
                    for iCCD, CCD in enumerate(self.CCDs):
                        for iQ, Q in enumerate(self.Quads):

                            datakeyX = 'FWHMX_%s_%s_CCD%s_Quad%s' % \
                                (bfkey,measkey,CCD,Q)

                            _valx.append(PT[datakeyX][ixblock])

                            datakeyY = 'FWHMY_%s_%s_CCD%s_Quad%s' % \
                                (bfkey,measkey,CCD,Q)

                            _valy.append(PT[datakeyY][ixblock])


                histoX = np.histogram(_valx, bins=20)
                Hdict[bfkey][measkey]['x']['fwhmx'] = histoX[0]
                Hdict[bfkey][measkey]['y']['fwhmx'] = histoX[1]
                Hdict[bfkey][measkey]['mean']['fwhmx'] = np.nanmean(_valx)

                histoY = np.histogram(_valy, bins=20)
                Hdict[bfkey][measkey]['x']['fwhmy'] = histoY[0]
                Hdict[bfkey][measkey]['y']['fwhmy'] = histoY[1]
                Hdict[bfkey][measkey]['mean']['fwhmy'] = np.nanmean(_valy)


        return Hdict

    def OLD_get_Hdict_FWHM(self, testname, fwhmkey, bfkey, measkey):
        """ """


        x = dict()
        y = dict()

        PT = self.ParsedTable[testname]

        _y = []

        for block in self.flight_blocks:

            ixblock = np.where(PT['BLOCK'] == block)[0][0]

            
            for iCCD, CCD in enumerate(self.CCDs):
                for iQ, Q in enumerate(self.Quads):

                    datakey = '%s_%s_%s_CCD%s_Quad%s' % \
                        (fwhmkey.upper(),bfkey.upper(),measkey.upper(),CCD,Q)

                    _y.append(PT[datakey][ixblock])

        histo = np.histogram(_y, bins=20)
        
        Hdict = dict(x=histo[1], y=histo[0])

        return Hdict



    def plot_XtalkMAP(self, XTALKs, **kwargs):
        """ """

        xtalkmap = plfpa.XtalkPlot(XTALKs, **kwargs)
        if 'figname' in kwargs:
            figname = kwargs['figname']
        else:
            figname = ''
        xtalkmap.render(figname=figname)

    def plot_PSF_HISTOS(self, fwhmxy_H, **psf_kwargs):
        """ """

        psfhistos = PsfPlot(fwhmxy_H, **psf_kwargs)
        if 'figname' in psf_kwargs:
            figname = psf_kwargs['figname']
        else:
            figname = ''
        psfhistos.render(figname=figname)

    def init_fignames(self):
        """ """

        if not os.path.exists(self.figspath):
            os.system('mkdir %s' % self.figspath)

        self.figs['XTALK_MAP_RT'] = os.path.join(self.figspath,
                                                 'XTALK_MAP_RT.png')

        self.figs['XTALK_RTvs800'] = os.path.join(self.figspath,
                                                  'XTALK_RT_vs_800nm.png')

        self.figs['XTALK_RTvs800_ABS'] = os.path.join(self.figspath,
                                                      'XTALK_RT_vs_800nm_abs.png')

        for wave in [590, 730, 880]:

            self.figs['XTALK_%ivs800' % wave] = os.path.join(self.figspath,
                                                             'XTALK_%inm_vs_800nm.png' % wave)

        for testname in self.testnames:
            self.figs['XTALK_MAP_%s' % testname] = os.path.join(self.figspath,
                                                                'XTALK_MAP_%s.png' % testname)

        for testname in self.testnames:

            self.figs['PSF_FWHM_MS_%s' % (testname,)] = \
                os.path.join(self.figspath,'PSF_FWHM_MS_%s.png' % (testname,))

            #for bfkey in ['bfe', 'nobfe']:
            #    self.figs['PSF_FWHMX_SLOPE_%s_%s' % (bfkey.upper(),testname)] = \
            #        os.path.join(self.figspath,'PSF_FWHMX_SLOPE_%s_%s.png' % (bfkey.upper(),testname))

            #    self.figs['PSF_FWHMX_MEAN_%s_%s' % (bfkey.upper(),testname)] = \
            #            os.path.join(self.figspath,'PSF_FWHMX_MEAN_%s_%s.png' % (bfkey.upper(),testname))

            #    self.figs['PSF_FWHMY_SLOPE_%s_%s' % (bfkey.upper(),testname)] = \
            #            os.path.join(self.figspath,'PSF_FWHMY_SLOPE_%s_%s.png' % (bfkey.upper(),testname))

            #    self.figs['PSF_FWHMY_MEAN_%s_%s' % (bfkey.upper(),testname)] = \
            #            os.path.join(self.figspath,'PSF_FWHMY_MEAN_%s_%s.png' % (bfkey.upper(),testname))

    def init_outcdpnames(self):
        """ """

        if not os.path.exists(self.cdpspath):
            os.system('mkdir %s' % self.cdpspath)

        for testname in self.testnames:
            self.outcdps['XTALK_MX_%s' % testname] = 'XTALK_MX_%s.json' % testname

    def _serialise_XTALKs(self,XTALKs):
        """ """
        blocks = XTALKs['blocks']
        dd = OrderedDict(blocks=blocks)

        for block in blocks:

            dd[block] = OrderedDict()

            for CCDso in self.CCDs:

                CCDsok = 'CCD%i' % CCDso

                dd[block][CCDsok] = OrderedDict()

                for Qso in self.Quads:

                    dd[block][CCDsok][Qso] = OrderedDict()

                    for CCDvi in self.CCDs:

                        CCDvik = 'CCD%i' % CCDvi
                        dd[block][CCDsok][Qso][CCDvik]=OrderedDict()

                        for Qvi in self.Quads:

                            inval = XTALKs[block][CCDsok][Qso][CCDvik][Qvi]

                            if len(inval)==0:
                                val=(None,None)
                            else:
                                val = ('%.3e' % inval['xtalk'],
                                    '%.3e' % inval['std_xtalk'])

                            dd[block][CCDsok][Qso][CCDvik][Qvi]=val

        return dd

    def verify_reqs(self, XTALKs):
        """ 
        for each victim channel, the maximum (dominant) coupling with any
            other channel shall be abs(c) <= 6x10-4.
        for each victim channel, only one other channel, the “dominant”, can
            have a coupling factor abs(c)> 7.6x10-5.
        """

        xtalks_VI = OrderedDict() # reformating the crosstalk matrix for easy
                                  # validadation
        xtalks_VI['blocks'] = XTALKs['blocks']

        for block in xtalks_VI['blocks']:  

            xtalks_VI[block] = OrderedDict()

            for CCDso in self.CCDs:
                CCDsok = 'CCD%i' % CCDso
                for Qso in self.Quads:
                    for CCDvi in self.CCDs:
                        CCDvik = 'CCD%i' % CCDvi
                        for Qvi in self.Quads:

                            _mx_elem = XTALKs[block][CCDsok][Qso][CCDvik][Qvi]

                            if len(_mx_elem) > 0:

                                ixtalk = _mx_elem['xtalk']

                                sok = '%s\\_%s' % (CCDsok, Qso)
                                vik = '%s\\_%s' % (CCDvik, Qvi)


                                if vik not in xtalks_VI[block]:
                                    xtalks_VI[block][vik] = OrderedDict()
                                    xtalks_VI[block][vik]['xtalk'] = []
                                    xtalks_VI[block][vik]['sok'] = []

                                xtalks_VI[block][vik]['xtalk'].append(ixtalk)
                                xtalks_VI[block][vik]['sok'].append(sok)


        validation = []
        req1 = 6.E-4
        req2 = 7.6E-5

        for block in xtalks_VI['blocks']:
            viks = xtalks_VI[block].keys()
            for vik in viks:
                ixorder = np.argsort(np.abs(xtalks_VI[block][vik]['xtalk']))[::-1]

                _xtalks = np.array(xtalks_VI[block][vik]['xtalk'])[ixorder]
                _soks = np.array(xtalks_VI[block][vik]['sok'])[ixorder]

                if np.abs(_xtalks[0]) > req1:
                    validation.append('REQ1 viol.: $|%.2e| >$ %.2e %s; source:%s; victim:%s\n' %\
                        (_xtalks[0], req1, block, _soks[0], vik))

                NREQ2 = len([val for val in _xtalks[1:] if np.abs(val)>req2])
                if NREQ2>1:
                    for i in range(1,len(_xtalks)):
                        if np.abs(_xtalks[i]) > req2:
                            validation.append('REQ2 viol.: $|%.2e| >$ %.2e %s; source:%s; victim:%s\n' %\
                                (_xtalks[i], req2, block, _soks[i], vik))
                        else:
                            break


        return validation



    def dump_aggregated_results(self):
        """ """

        if self.report is not None:
            self.report.add_Section(keyword='dump', Title='Aggregated Results', level=0)

            self.add_DataAlbaran2Report()

        function, module = utils.get_function_module()
        CDP_header = self.CDP_header.copy()
        CDP_header.update(dict(function=function, module=module))
        CDP_header['DATE'] = self.get_time_tag()

        doXtalk = False
        doPsf = True

        if doXtalk:

            # XTALK MAP (ROE-TAB)

            XTALKs_RT = self.products['XTALK_RT'].copy()

            figkey0 = 'XTALK_MAP_RT'
            figname0 = self.figs[figkey0]

            self.plot_XtalkMAP(XTALKs_RT, **dict(
                scale='ADU',
                showvalues=False,
                title='XTALK [ADU] - ROE-TAB',
                figname=figname0))

            captemp0 = 'ROE-TAB: Cross-Talk matrix [in ADU] obtained '+\
                'using the ROE-TAB [electrical injection].'

            if self.report is not None:
                self.addFigure2Report(figname0, 
                    figkey=figkey0, 
                    caption= captemp0, 
                    texfraction=0.7)

            # XTALK MAPS (OPTICAL)

            for testname in self.testnames:

                XTALKs = self.get_XTALKDICT_from_PT(testname)

                figkey1 = 'XTALK_MAP_%s' % testname
                figname1 = self.figs[figkey1]

                stestname = testname.replace('_', '\_')
                self.plot_XtalkMAP(XTALKs, **dict(
                    scale='ADU',
                    showvalues=False,
                    title='%s: XTALK [ADU]' % stestname,
                    figname=figname1))

                captemp1 = '%s: Cross-Talk as amplitude in ADU of ghost in victim channel,'+\
                ' for a saturating signal in the source channel.'

                if self.report is not None:
                    self.addFigure2Report(figname1, 
                            figkey=figkey1, 
                            caption= captemp1 % stestname, 
                            texfraction=0.7)


                if self.report is not None and testname =='PSF01_800':

                    req_ver_txt = [
                    '\nTEST: %s\n' % stestname,
                    'Req.:\n',
                    '    dominant coupling $|c1|<=$ 6E-4\n',
                    '    only 1 sub-dominant coupling 6E-4$>|c2|>$7.6E-5\n'
                    ]

                    req_ver_txt += self.verify_reqs(XTALKs)

                    self.report.add_Text(req_ver_txt)

                xt_header = OrderedDict()
                xt_header['title'] = 'XTALK_MATRIX:%s' % testname
                xt_header.update(CDP_header)

                xt_cdp = cdp.Json_CDP(rootname=self.outcdps['XTALK_MX_%s' % testname],
                    path=self.cdpspath)

                XTALKs_CDP_MX = self._serialise_XTALKs(XTALKs)

                xt_cdp.ingest_inputs(data=XTALKs_CDP_MX,
                        header = xt_header,
                        meta=dict(units='ADIM',
                            structure='block:ccd_source:Q_source:'+\
                            'ccd_victim:Q_victim:(xtalk,e_xtalk)'))

                xt_cdp.savehardcopy()



            # XTALK: 800-optical vs. RT  (with SIGN)

            XT_RTvs800 = self._get_XYdict_XT('PSF01_800', 'RT', mode='sign')

            figkey2 = 'XTALK_RTvs800'
            figname2 = self.figs[figkey2]


            XTkwargs = dict(
                title='Cross-Talk Direct Comparison',
                doLegend=False,
                xlabel='Xtalk - Opt. 800nm',
                ylabel='Xtalk - ROE-TAB',
                xlim=[-20, 40],
                ylim=[-20, 40],
                figname=figname2)

            BLOCKcolors = cm.rainbow(np.linspace(0, 1, len(self.flight_blocks)))

            xtcorekwargs = dict()
            for jblock, block in enumerate(self.flight_blocks):
                jcolor = BLOCKcolors[jblock]
                xtcorekwargs['%s' % (block,)] = dict(linestyle='',
                                                     marker='.', color=jcolor)

            xtcorekwargs['oneone'] = dict(linestyle='--', marker='', color='k')

            XTkwargs['corekwargs'] = xtcorekwargs

            self.plot_XY(XT_RTvs800, **XTkwargs)


            if self.report is not None:

                captemp2 = '%s: Cross-talk coupling factor measured using the ROE-TAB, vs. using '+\
                'optical stimulation (at 800 nm). Sign of coupling preserved.'

                self.addFigure2Report(figname2, 
                            figkey=figkey2, 
                            caption=captemp2 % stestname, 
                            texfraction=0.7)

            # XTALK: 800-optical vs. RT (ABS-VALUE)

            XT_RTvs800_abs = self._get_XYdict_XT('PSF01_800', 'RT', mode='abs')

            figkey3 = 'XTALK_RTvs800_ABS'
            figname3 = self.figs[figkey3]

            XTABSkwargs = dict(
                title='Cross-Talk Comparison - ABS. Value',
                doLegend=False,
                xlabel='Abs(Xtalk - Opt. 800 nm)',
                ylabel='Abs(Xtalk - ROE-TAB)',
                xlim=[-20, 50],
                ylim=[-20, 50],
                figname=figname3)

            XTABSkwargs['corekwargs'] = xtcorekwargs

            self.plot_XY(XT_RTvs800_abs, **XTABSkwargs)

            if self.report is not None:

                captemp3 = '%s: Cross-talk coupling factor (in absolute value) measured '+\
                'using the ROE-TAB, vs. using optical stimulation (at 800 nm).'

                self.addFigure2Report(figname3, 
                            figkey=figkey3, 
                            caption=captemp3 % stestname, 
                            texfraction=0.7)

            # XTALK: 800-optical vs. OTHER-opt (with SIGN)

            for wave in [590, 730, 880]:

                XT_NMvs800 = self._get_XYdict_XT('PSF01_800', 'PSF01_%i' % wave, mode='sign')

                figkey4 = 'XTALK_%ivs800' % wave
                figname4 = self.figs[figkey4]

                XTNMvs800kwargs = dict(
                    title='Cross-Talk Comparison - With Sign',
                    doLegend=False,
                    xlabel='Xtalk - Opt. 800 nm',
                    ylabel='Xtalk - Opt. %i nm' % wave,
                    xlim=[-20, 50],
                    ylim=[-20, 50],
                    figname=figname4)

                XTNMvs800kwargs['corekwargs'] = xtcorekwargs

                self.plot_XY(XT_NMvs800, **XTNMvs800kwargs)

                if self.report is not None:

                    captemp4 = 'Cross-talk coupling factor measured at %i nm '+\
                    'vs. 800 nm.'

                    self.addFigure2Report(figname4, 
                            figkey=figkey4, 
                            caption=captemp4 % wave,
                            texfraction=0.7)

                    std_comparisons = []

                    for block in self.flight_blocks:
                        std_comparisons.append(np.nanstd(XT_NMvs800['y'][block]-\
                                XT_NMvs800['x'][block]))

                    avg_std = np.nanmedian(std_comparisons)

                    self.report.add_Text('\nAvg. STD of comparison (%i vs. 800 nm): %.2e ADU' %\
                            (wave, avg_std))

        # PSF analysis

        if doPsf:

            # Histogram of slopes of FWHWMx/y vs. fluence w/o BF correction

            for testname in self.testnames:

                stestname = testname.replace('_','\_')

                figkey5 = 'PSF_FWHM_MS_%s' % testname
                figname5 = self.figs[figkey5]

                fwhmxy_H = self._get_Hdict_FWHM(testname)

                psf_kwargs = dict(title='%s: FWHMxy trends' % (stestname,),
                    doLegend=True,
                    xlabel='FWHMX [pixels]',
                    ylabel='N',
                    #xlim=[],
                    #ylim=[],
                    figname=figname5)

                self.plot_PSF_HISTOS(fwhmxy_H, **psf_kwargs)

                captemp5 = '%s'

                if self.report is not None:
                    self.addFigure2Report(figname5, 
                        figkey=figkey5, 
                        caption= captemp5 % (stestname, ),
                        texfraction=1.0)


                skip = True

                if not skip:

                    for bfkey in ['bfe', 'nobfe']:


                        # Mean FWHMz

                        figkey5 = 'PSF_FWHMX_MEAN_%s_%s' % (bfkey.upper(),testname)
                        figname5 = self.figs[figkey5]

                        fwhmx_mean_H = self._get_Hdict_FWHM(testname, 'fwhmx', bfkey, 'mean')

                        XMEAN_kwargs = dict(title='%s: FWHMX-mean (%s)' % (stestname,bfkey.upper()),
                            doLegend=False,
                            xlabel='FWHMX [pixels]',
                            ylabel='N',
                            #xlim=[],
                            #ylim=[],
                            figname=figname5)

                        self.plot_HISTO(fwhmx_mean_H, **XMEAN_kwargs)

                        captemp5 = '%s: (%s)'

                        if self.report is not None:
                            self.addFigure2Report(figname5, 
                                figkey=figkey5, 
                                caption= captemp5 % (stestname, bfkey),
                                texfraction=0.9)

                        figkey6 = 'PSF_FWHMY_MEAN_%s_%s' % (bfkey.upper(),testname)
                        figname6 = self.figs[figkey6]

                        fwhmy_mean_H = self._get_Hdict_FWHM(testname, 'fwhmy', bfkey, 'mean')

                        YMEAN_kwargs = dict(title='%s: FWHMY-mean (%s)' % (stestname,bfkey.upper()),
                            doLegend=False,
                            xlabel='FWHMY [pixels]',
                            ylabel='N',
                            #xlim=[],
                            #ylim=[],
                            figname=figname6)

                        self.plot_HISTO(fwhmy_mean_H, **YMEAN_kwargs)

                        captemp6 = '%s: (%s)'

                        if self.report is not None:
                            self.addFigure2Report(figname6, 
                                figkey=figkey6, 
                                caption= captemp6 % (stestname, bfkey),
                                texfraction=0.9)
                        

                        # SLOPES - FWHMz vs. fluence


                        figkey7 = 'PSF_FWHMX_SLOPE_%s_%s' % (bfkey.upper(),testname)
                        figname7 = self.figs[figkey7]

                        fwhmx_slope_H = self._get_Hdict_FWHM(testname, 'fwhmx', bfkey, 'slope')

                        XSLOPE_kwargs = dict(title='%s: FWHMX-Slope (%s)' % (stestname, bfkey.upper()),
                            doLegend=False,
                            xlabel='DeltaFWHMX/fluence [pixels/(10 kADU)]',
                            ylabel='N',
                            #xlim=[],
                            #ylim=[],
                            figname=figname7)

                        self.plot_HISTO(fwhmx_slope_H, **XSLOPE_kwargs)

                        captemp7 = '%s: (%s)'

                        if self.report is not None:
                            self.addFigure2Report(figname7, 
                                figkey=figkey7, 
                                caption= captemp7 % (stestname, bfkey),
                                texfraction=0.9)

                        figkey8 = 'PSF_FWHMY_SLOPE_%s_%s' % (bfkey.upper(),testname)
                        figname8 = self.figs[figkey8]

                        fwhmy_slope_H = self._get_Hdict_FWHM(testname, 'fwhmy', bfkey, 'slope')

                        YSLOPE_kwargs = dict(title='%s: FWHMY-Slope (%s)' % (stestname,bfkey.upper()),
                            doLegend=False,
                            xlabel='DeltaFWHMY/fluence [pixels/(10 kADU)]',
                            ylabel='N',
                            #xlim=[],
                            #ylim=[],
                            figname=figname8)

                        self.plot_HISTO(fwhmy_slope_H, **YSLOPE_kwargs)

                        captemp8 = '%s: (%s)'

                        if self.report is not None:
                            self.addFigure2Report(figname8, 
                                figkey=figkey8, 
                                caption= captemp8 % (stestname, bfkey),
                                texfraction=0.9)