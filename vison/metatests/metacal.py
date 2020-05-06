#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 16:06:40 2019

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import os
import copy
import numpy as np
import os
from astropy import table
import glob
import gc
import string as st
import pandas as pd

from vison import __version__
from vison.support import files
from vison.support import vjson
from vison.datamodel import core as vcore
from vison.datamodel import cdp as cdpmod
from vison.support import vcal
from vison.fpa import fpa as fpamod
from vison.plot import plots_fpa as plfpa
from vison.plot import baseplotclasses as basepl
from vison.support.report import Report

from matplotlib import pyplot as plt
plt.switch_backend('TkAgg')
from matplotlib.colors import Normalize

# END IMPORT


class MetaCal(object):

    def __init__(self, **kwargs):
        """ """

        if 'design' in kwargs:
            design = kwargs['design']
        else:
            design = 'final'
        self.fpa = fpamod.FPA(design)
        self.blocks = copy.deepcopy(self.fpa.all_blocks)
        self.flight_blocks = copy.deepcopy(self.fpa.flight_blocks)
        # self.blocks = ['BORN','CURIE','DIRAC','FOWLER','GUYE','KRAMERS'] # TESTS
        # self.flight_blocks = ['BORN','CURIE','DIRAC','FOWLER','GUYE','KRAMERS'] # TESTS
        self.CCDs = [1, 2, 3]
        self.Quads = ['E', 'F', 'G', 'H']

        self.NSLICES_FPA = self.fpa.NSLICES
        self.NCOLS_FPA = self.fpa.NCOLS

        if 'vcalfile' in kwargs:
            self.vcalfile = kwargs['vcalfile']
        else:
            self.vcalfile = ''
        if 'respathroot' in kwargs:
            self.respathroot = kwargs['respathroot']
        else:
            self.respathroot = ''
        if 'jsonf' in kwargs:
            self.jsonf = kwargs['jsonf']
        else:
            self.jsonf = ''

        if 'testkey' in kwargs:
            self.testkey = kwargs['testkey']
        else:
            self.testkey = ''

        if 'outparent' in kwargs:
            outparent = kwargs['outparent']
        else:
            outparent = ''
        self.outpathroot = os.path.join(outparent, '%s_FPA' % self.testkey.upper())

        self.inventory = OrderedDict()
        self.results = OrderedDict()
        self.products = OrderedDict()
        self.ParsedTable = None
        self.roeVCals = dict()
        self.init_roeVCals()
        self.cdps = OrderedDict()
        self.outcdps = OrderedDict()
        self.figs = dict()
        self.report = None
        self.figspath = os.path.join(self.outpathroot, 'figs')
        self.cdpspath = os.path.join(self.outpathroot, 'cdps')

        self.CDP_header = OrderedDict()
        self.CDP_header['vison']=__version__
        self.CDP_header['fpa_design'] = design
        self.CDP_header['FPA'] = self.fpa.FPA_MAP.copy()
        self.CDP_header['vcalfile'] = os.path.split(self.vcalfile)[-1]

    def get_time_tag(self):
        from vison.support.vistime import get_time_tag
        return get_time_tag()

    def init_fignames(self):
        raise NotImplementedError("Subclass must implement abstract method")
    
    def init_outcdpnames(self):
        raise NotImplementedError("Subclass must implement abstract method")

    def run(self, doLoad=True, doParse=True, doDump=True, doReport=True):
        """ """

        if not os.path.exists(self.outpathroot):
            os.system('mkdir -p %s' % self.outpathroot)

        dictiopick = os.path.join(self.outpathroot,
                                  '%s_dictionary.pick' % self.testkey.lower())


        if doLoad:
            print('Loading block(s) results...')
            self.load_block_results()
            files.cPickleDumpDictionary(self.inventory, dictiopick)
        else:
            print('Re-loading block(s) results')
            self.inventory = files.cPickleRead(dictiopick)

        parsedpick = os.path.join(self.outpathroot,
                                  '%s_parsed.pick' % self.testkey.lower())

        if doParse:
            print('Parsing results')
            for testname in self.testnames:
                self.parse_test_results(testname)
            parsedbundle = dict(PT=self.ParsedTable,
                                products=self.products)
            files.cPickleDumpDictionary(parsedbundle, parsedpick)

        else:
            print('Re-loading parsed results')
            parsedbundle = files.cPickleRead(parsedpick)
            self.ParsedTable = parsedbundle['PT'].copy()
            self.products = parsedbundle['products'].copy()

        if doReport:
            self.report = Report(TestName=self.testkey.upper(), 
                            Model='FM',
                            Reference='7-XXX',
                            doDraft = False)
        else:
            self.report = None

        if doDump:
            self.dump_aggregated_results()

        if doReport:
            cleantexafter = False # hard-wired by now
            reportroot = '%s_TR' % self.testkey.upper()
            #self.report()
            outfiles = self.report.doreport(
                reportroot,
                cleanafter=cleantexafter,
                silent=True)  # commented on TESTS

            for outfile in outfiles:
                os.system('mv %s %s/' % (outfile, self.outpathroot))


    def init_roeVCals(self):

        for block in self.blocks:

            ROE_SN = fpamod.ROE_SNs[block]

            roevcal = vcal.RoeVCal(self.vcalfile)
            try:
                roevcal.load_VCal(ROE_SN)
            except IOError:
                continue

            self.roeVCals[block] = copy.deepcopy(roevcal)

    def _update_inventory(self, block, testline):

        test = testline[0]
        session = testline[1]
        repeat = testline[2]

        sesspath = os.path.join(self.respathroot, block, 'ANALYSIS', 'results_atCALDATA',
                                session)

        #alltestpaths = glob.glob(os.path.join(sesspath,'%s*' % test))
        tester = glob.glob(os.path.join(sesspath, '%s.[0-9]' % test))

        if len(tester) > 0:
            testpath = '%s.%i' % (test, repeat)
        else:
            testpath = test

        respath = os.path.join(sesspath, testpath)

        DDfile = os.path.join(respath, '%s_DataDict.pick' % test)

        try:
            assert os.path.exists(respath)
        except BaseException:
            stop()
        try:
            assert os.path.exists(DDfile)
        except BaseException:
            stop()

        if block not in self.inventory:
            self.inventory[block] = OrderedDict()
        if test not in self.inventory[block]:
            self.inventory[block][test] = []

        self.inventory[block][test].append(OrderedDict())

        ii = self.inventory[block][test][-1]

        ii['DD'] = DDfile
        try:
            dd = files.cPickleRead(DDfile)
        except BaseException:
            print('Could not load %s\n\n' % DDfile)
            raise RuntimeError

        ii['dd'] = copy.deepcopy(dd)
        ii['resroot'] = respath
        ii['session'] = session
        ii['repeat'] = repeat

        return None

    def load_block_results(self, inventoryfile=None):

        if inventoryfile is None:
            inventoryfile = self.jsonf
        inventraw = vjson.load_jsonfile(inventoryfile, useyaml=True)

        for block in self.blocks:
            if block in inventraw['inventory'].keys():
                rawcargo = inventraw['inventory'][block]

                for testline in rawcargo:
                    self._update_inventory(block, testline)

        return None

    def dump_aggregated_results(self):
        raise NotImplementedError("Subclass must implement abstract method")

    def stack_dd(self, dd, cols, index2stack, indices2keep, stacker='median'):
        """ """

        #indices = copy.deepcopy(dd.indices)

        stkdd = vcore.DataDict()

        stkdd.compliances = dd.compliances.copy()
        stkdd.flags = copy.deepcopy(dd.flags)
        stkdd.meta = dd.meta.copy()
        stkdd.products = dd.products.copy()

        stacker_functions = dict(median=np.median,
                                 mean=np.mean)

        fstacker = stacker_functions[stacker]

        for ixc, icol in enumerate(cols):

            idtype = dd.mx[icol][:].dtype

            iindices = dd.mx[icol].indices

            niindices = []

            niindices = [index for index in iindices
                         if index.name in indices2keep]

            # STACKING

            _ix2stack = iindices.names.index(index2stack)

            if idtype.char not in ['S', 'O']:

                imx = fstacker(dd.mx[icol][:], axis=_ix2stack, keepdims=True)

            else:
                slicer = []
                for _ix, _i in enumerate(iindices):
                    if _i.name == index2stack:
                        slicer.append([0])
                    else:
                        slicer.append(slice(_i.len))
                imx = dd.mx[icol][tuple(slicer)]

            # TRIMMING UNNECESSARY AXES

            if idtype.char not in ['S', 'O']:

                _ix2ax = [_ix for _ix, ix in enumerate(iindices) if
                          ix.name not in indices2keep]

                imx = fstacker(imx, axis=_ix2ax, keepdims=False)
            else:

                slicer = []

                for _ix, _i in enumerate(iindices):
                    if _i.name not in indices2keep and _i.name != index2stack:
                        slicer.append(0)

                imx = imx[tuple(slicer)]

            # INDEXING

            niindices = []

            for ix in iindices:
                if ix.name not in indices2keep:
                    pass
                else:
                    if ix.name == index2stack:
                        _ix = copy.deepcopy(ix)
                        _ix.len = 1
                        _ix.vals = [ix.vals[0]]
                        niindices.append(_ix)
                    else:
                        niindices.append(copy.deepcopy(ix))

            stkdd.addColumn(imx, name=icol, indices=niindices)

        return stkdd

    def stackTables(self, t1, t2):
        """ """
        return table.vstack([t1, t2])

    def parse_single_test_gen(self, jrep, block, testname, inventoryitem):
        """ """

        IndexS = vcore.vMultiIndex([vcore.vIndex('ix', vals=[0])])

        idd = copy.deepcopy(inventoryitem['dd'])

        sidd = self.stack_dd(idd, self.incols,
                             indices2keep=['ix', 'CCD', 'Quad'],
                             index2stack='ix',
                             stacker='median')

        sidd.dropColumn('test')

        # rename the CCD index values in sidd, for convenience

        sidd.indices[sidd.indices.names.index('CCD')].vals = self.CCDs

        for col in sidd.colnames:
            if 'CCD' in sidd.mx[col].indices.get_names():
                _i = sidd.mx[col].indices.names.index('CCD')
                sidd.mx[col].indices[_i].vals = self.CCDs

        # CALIBRATED HK

        try:
            roeVCal = self.roeVCals[block]
        except KeyError:
            print('Voltage calibrations for block %s not found!' % block)
            roeVCal = None

        if roeVCal is not None:

            for commcal_key in ['IDL', 'IDH']:

                for iCCD, CCD in enumerate(self.CCDs):
                    for Q in self.Quads:

                        cEkey = '%s_CCD%s_Quad%s_CAL' % (commcal_key, CCD, Q)

                        EV = sidd.mx[commcal_key][0][iCCD]
                        EVcal = roeVCal.fcal_ELVIS_script(EV, commcal_key, CCD, Q)

                        sidd.addColumn(np.zeros(1, dtype=float) + EVcal,
                                       cEkey, IndexS)

            for hkcal_key in ['OD', 'RD', 'IG1', 'IG2']:

                for CCD in self.CCDs:
                    for Q in self.Quads:

                        rHKkey = roeVCal.get_HKkey(hkcal_key, CCD, Q)
                        HKkey = 'HK_%s' % rHKkey.upper()
                        cHKkey = '%s_CAL' % HKkey

                        if HKkey in sidd.mx:

                            HKV = sidd.mx[HKkey][0]

                            HKVcal = roeVCal.fcal_HK(HKV, hkcal_key, CCD, Q)

                            sidd.addColumn(np.zeros(1, dtype=float) + HKVcal,
                                           cHKkey, IndexS)

        else:

            dummy_roeVCal = vcal.RoeVCal()

            cHKkeys = []

            for cal_key in ['OD', 'RD', 'IG1', 'IG2']:

                for CCD in self.CCDs:
                    for Q in self.Quads:

                        rHKkey = dummy_roeVCal.get_HKkey(cal_key, CCD, Q)
                        HKkey = 'HK_%s' % rHKkey.upper()
                        cHKkey = '%s_CAL' % HKkey

                        cHKkeys.append(cHKkey)

            for cHKkey in cHKkeys:
                sidd.addColumn(np.zeros(1, dtype=float) + np.nan,
                               cHKkey, IndexS)

        return sidd

    def parse_test_results(self, testname):
        """ """

        for iblock, block in enumerate(self.blocks):

            try:
                Nreps = len(self.inventory[block][testname])
            except KeyError:
                print('block %s not found!' % block)
                continue

            for jrep in range(Nreps):

                print('Parsing %s:%s' % (block, jrep + 1))

                inventoryitem = self.inventory[block][testname][jrep]

                sit = self.parse_single_test(jrep, block, testname, inventoryitem)

                # MERGING WITH PREVIOUS DDs

                if (iblock == 0) and (jrep == 0):
                    pt = copy.deepcopy(sit)
                else:
                    pt = self.stackTables(pt, sit)

        self.ParsedTable[testname] = pt

    def add_DataAlbaran2Report(self):
        """Adds a data delivery note to the report."""

        albaran_dict = OrderedDict()
        cols = ['BLOCK','TEST','SESSION','DAY-FOLDER','OBSIDS']
        for col in cols:
            albaran_dict[col] = []

        for iblock, block in enumerate(self.blocks):

            for testname in self.testnames:

                try:
                    Nreps = len(self.inventory[block][testname])
                except KeyError:
                    print('block %s not found!' % block)
                    continue

                for jrep in range(Nreps):
                    
                    iitem = self.inventory[block][testname][jrep]

                    dayfolder = os.path.split(iitem['dd'].meta['inputs']['datapath'])[-1]
                    sdayfolder = st.replace(dayfolder,'_','\\_')
                    

                    OBSID_min = iitem['dd'].meta['data_inventory']['ObsID'].min()
                    OBSID_max = iitem['dd'].meta['data_inventory']['ObsID'].max()

                    OBSID_lims = '%i-%i' % (OBSID_min, OBSID_max)

                    albaran_dict['BLOCK'].append(block)
                    albaran_dict['TEST'].append(testname)
                    albaran_dict['SESSION'].append(iitem['session'])
                    albaran_dict['DAY-FOLDER'].append(sdayfolder)
                    albaran_dict['OBSIDS'].append(OBSID_lims)

        albaran_df = pd.DataFrame(albaran_dict)
        
        kwargs = dict(multicolumn=True, multirow=True, longtable=True, index=False)
        
        tex = albaran_df.to_latex(**kwargs)

        ncols = len(cols)

        caption = 'Data Inventory.'
        wtex = cdpmod.wraptextable(tex, ncols, caption, fitwidth=True, tiny=True, 
            longtable=True)

        self.report.add_Text(wtex)



    
    def get_ixblock(self, PT, block):
        return np.where(PT['BLOCK'].data==block)
    
    def get_stat_from_FPAMAP(self, M, numpystat):
        """ """

        vals = []

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                for Q in self.Quads:
                    vals.append(M[Ckey][Q])

        return numpystat(vals)

    def get_FPAMAP_from_PT(self, PT, extractor):
        """ """

        M = OrderedDict()

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                M[Ckey] = OrderedDict()

                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]

                for Q in self.Quads:

                    M[Ckey][Q] = extractor(PT, block, CCDk, Q)

        return M

    def plot_SimpleMAP(self, MAPdict, **kwargs):
        """ """

        VALs = []
        for ckey in MAPdict.keys():
            for Q in self.Quads:
                VALs.append(MAPdict[ckey][Q])

        normfunction = Normalize(vmin=np.nanmin(VALs), vmax=np.nanmax(VALs), clip=False)

        _kwargs = dict(doColorbar=True,
                       corekwargs=dict(
                           norm=normfunction
                       ))

        _kwargs.update(kwargs)

        if 'figname' in _kwargs:
            figname = _kwargs['figname']
        else:
            figname = ''


        with plfpa.FpaHeatMap(MAPdict, **_kwargs) as heatmap:
            heatmap.render(figname=figname)
        #heatmap = None

        gc.collect()


    def _get_plot_gen(self, plotclass):
        """ """

        def plot_gen(self, datadict, **kwargs):
            _kwargs = dict()
            _kwargs.update(kwargs)

            if 'figname' in _kwargs:
                figname = _kwargs['figname']
            else:
                figname = ''

            with plotclass(datadict, **_kwargs) as plot:
                plot.render(figname=figname)

        return plot_gen

        gc.collect()


    def plot_XY(self, XYdict, **kwargs):
        """ """
        plotobj = self._get_plot_gen(basepl.XYPlot)
        plotobj(self, XYdict, **kwargs)
        

    def plot_HISTO(self, Hdict, **kwargs):
        """ """
        plotobj = self._get_plot_gen(basepl.HistoPlot)
        plotobj(self, Hdict, **kwargs)


    def plot_XYMAP(self, XYMAP, **kwargs):

        plotobj = self._get_plot_gen(plfpa.FpaPlotYvsX)
        plotobj(self, XYMAP, **kwargs)


    def plot_XYCCD(self, XYCCD, **kwargs):

        plotobj = self._get_plot_gen(basepl.CCD2DPlotYvsX)
        plotobj(self, XYCCD, **kwargs)


    def plot_ImgFPA(self, BCdict, **kwargs):

        plotobj = self._get_plot_gen(plfpa.FpaImgShow)
        plotobj(self, BCdict, **kwargs)


    def add_StdQuadsTable2Report(self, extractor=None, Matrix=None, cdp=None, cdpdict=None):
        """ """

        assert (extractor is None) != (Matrix is None)
        
        if extractor is not None:
            def _get_val(Ckey,Q):
                return extractor(Ckey,Q)
        elif Matrix is not None:
            def _get_val(Ckey,Q):
                return Matrix[Ckey][Q]
        

        defcdpdict = dict(TBkey = 'TB',
            meta = dict(),
            CDP_header = dict(),
            header_title = 'generic title',
            CDP_KEY = 'UNK',
            caption = 'CAPTION PENDING',
            valformat = '%.2e')

        if cdpdict is not None:

            assert isinstance(cdpdict, dict)
            
            defcdpdict.update(cdpdict)
        
        TBkey = defcdpdict['TBkey']
        meta = defcdpdict['meta'].copy()
        CDP_header = defcdpdict['CDP_header'].copy()
        header_title = defcdpdict['header_title']
        CDP_KEY = defcdpdict['CDP_KEY']
        valformat = defcdpdict['valformat']
        caption = defcdpdict['caption']
        

        NCCDs = len(self.CCDs)
        NBlocks = len(self.fpa.flight_blocks)
        NP = NBlocks * NCCDs

        TB = OrderedDict()
        TB['CCDID'] = np.zeros(NP, dtype='S4')
        TB['BLOCK'] = np.zeros(NP, dtype='S4')
        TB['CCD'] = np.zeros(NP, dtype='S4')

        _dtype = np.array(_get_val('C_11','E')).dtype
        if 'S' in _dtype.str:
            Qdtype = 'S30'
        else:
            Qdtype = 'float32'

        for Q in self.Quads:
            TB[Q] = np.zeros(NP, dtype=Qdtype)

        for jY in range(self.NSLICES_FPA):
            for iX in range(self.NCOLS_FPA):

                Ckey = 'C_%i%i' % (jY + 1, iX + 1)

                locator = self.fpa.FPA_MAP[Ckey]
                block = locator[0]
                CCDk = locator[1]

                indx = jY * self.NCOLS_FPA + iX

                TB['CCDID'][indx] = Ckey
                TB['BLOCK'][indx] = block[0:4]
                TB['CCD'][indx] = CCDk

                for iQ, Q in enumerate(self.Quads):
                    TB[Q][indx] = _get_val(Ckey, Q)
        
        TB_ddf = OrderedDict()
        TB_ddf[TBkey] = pd.DataFrame.from_dict(TB)

        if cdp is not None:
            cdp.path = self.inputs['subpaths']['products']
            cdp.ingest_inputs(
                data=TB_ddf.copy(),
                meta=meta.copy(),
                header=CDP_header.copy()
            )

            cdp.init_wb_and_fillAll(header_title=header_title)
            self.save_CDP(cdp)
            self.pack_CDP_to_dd(cdp, CDP_KEY)
        else:
            cdp = cdpmod.Tables_CDP()
            cdp.ingest_inputs(
                data=TB_ddf.copy(),
                meta=meta.copy(),
                header=CDP_header.copy()
            )

        def fstr(x): return '%s' % x

        def ff(x): return valformat % x

        her_formatters = [fstr, fstr, fstr]
        for iQ, Q in enumerate(self.Quads):
            her_formatters.append(ff)

        nicecaption = st.replace(caption, '_', '\_')
        Ttex = cdp.get_textable(sheet=TBkey, caption=nicecaption,
                                fitwidth=True,
                                tiny=True,
                                formatters=her_formatters)

        self.report.add_Text(Ttex)

        return TB, cdp


    def addFigure2Report(self, png, figkey, caption='', texfraction=0.7):
        """ """

        if self.report is None:
            return

        assert os.path.exists(png)
        eps = '%s.eps' % os.path.splitext(png)[0]
        os.system('convert %s %s' % (png, eps))
        self.report.add_Figure(eps, texfraction=texfraction,
            caption=caption,
            label=figkey)

    def FITSify_CDP_header(self, CDP_header):
        """ """
        
        outCDP_header = CDP_header.copy()

        if 'fpa_design' in outCDP_header:
            fpa_design = outCDP_header.pop('fpa_design')
            outCDP_header['fpadesi'] = fpa_design

        if 'FPA' in outCDP_header:
            rawFPA = outCDP_header.pop('FPA')
            
            for jY in range(self.NSLICES_FPA):
                for iX in range(self.NCOLS_FPA):
                    Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                    outCDP_header[Ckey] = rawFPA[Ckey][-1]
                    
        return outCDP_header



