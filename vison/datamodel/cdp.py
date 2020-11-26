#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Classes to store Calibration Data Products.

:History:
Created on Tue Feb 27 10:58:42 2018

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import copy
import os
import string as st
import numpy as np
from astropy.io import fits as fts
import pandas as pd
from astropy import table as atable

from vison import __version__
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.datamodel import ccd as ccdmod
from vison.datamodel import fpa_dm
from vison.support.excel import ReportXL
from vison.support import vjson
# END IMPORT


def loadCDPfromPickle(pickf):
    """Function to load a CDP from a pickle file."""
    cdp = cPickleRead(pickf)
    return cdp

def wraptextable(tex, ncols=1, caption='', fitwidth=False, tiny=False, longtable=False):
    """Auxiliary function to Tables_CDP class"""

    tex = tex.split('\n')

    if fitwidth:

        beglongtabu = '\\begin{longtabu} to \\textwidth {|%s}' % (ncols * 'X|',)
        endlongtabu = '\end{longtabu}'
        tex[0] = beglongtabu
        tex[-2] = endlongtabu

    if not longtable:
        tex = ['\\begin{table}[!htb]'] + tex + ['\end{table}']

    if tiny:
        tex = ['\\tiny'] + tex + ['\\normalsize']

    #tex = ['\\begin{table}[!htb]','\center']+tex+['\end{table}']

    if caption != '':

        if fitwidth:
            ixendlongtabu = np.where(['\end{longtabu}' in item for item in tex])
            tex.insert(ixendlongtabu[0][-1], r'\caption{%s}' % caption)
        elif not fitwidth and longtable:
            ixendlongtable = np.where(['\end{longtable}' in item for item in tex])
            tex.insert(ixendlongtable[0][-1], r'\caption{%s}' % caption)
        elif not fitwidth and not longtable:
            ixendtable = np.where(['\end{table}' in item for item in tex])
            tex.insert(ixendtable[0][-1], r'\caption{%s}' % caption)

    return tex


class CDP(object):
    """Parent CDP Class."""

    rootname = 'Unknown'
    path = ''
    header = OrderedDict()
    meta = OrderedDict()
    data = None


    def __init__(self, *args, **kwargs):
        """ """

        defaults = dict(rootname='Unknown', path='', header=OrderedDict(),
                        meta=OrderedDict(), data=None)
        defaults.update(args)
        defaults.update(kwargs)

        for key in list(defaults.keys()):
            self.__setattr__(key, defaults[key])

    def savetopickle(self, pickf=''):
        """ """
        if pickf == '':
            pickf = os.path.join(self.path, '%s.pick' % self.rootname)

        outdict = copy.copy(self.__dict__)

        if 'report' in outdict:
            report = outdict['report']

            if isinstance(outdict['report'], ReportXL):
                report = outdict['report']
                outdict['report'] = None  # gives trouble... HACK
                self.__dict__['report'] = report

        cPickleDumpDictionary(outdict, pickf)

    def loadfrompickle(self, pickf=''):
        """ """
        if pickf == '':
            pickf = os.path.join(self.path, '%s.pick' % self.rootname)
        saved = cPickleRead(pickf)
        self.__dict__.update(saved)

    def savehardcopy(self, filef=''):
        """ """
        pass  # HACK
        # raise NotImplementedError(
        #    'Subclass implements abstract method (if needed).')


class Tables_CDP(CDP):
    """Table CDP. Can export to excel."""

    def __init__(self, *args, **kwargs):
        """ """
        super(Tables_CDP, self).__init__(*args, **kwargs)
        self.figs = OrderedDict()
        self.data = OrderedDict()
        self.report = None

    def ingest_inputs(self, data, meta=None, header=None, figs=None):
        """ """
        if meta is None:
            meta = OrderedDict()
        if header is None:
            header = OrderedDict()
        if figs is None:
            figs = OrderedDict()

        self.header.update(header)
        self.meta.update(meta)
        self.data.update(data)
        self.figs.update(figs)

    def init_workbook(self):
        """ """
        self.report = ReportXL(OrderedDict())
        self.report.wb.create_sheet('Header', 0)
        self.report.wb.create_sheet('Meta', 1)

        sheetcounter = 2

        for i, k in enumerate(self.data.keys()):
            sheetcounter += i
            self.report.wb.create_sheet(k, sheetcounter)

        if len(self.figs) > 0:
            self.report.wb.create_sheet('figs', sheetcounter + 1)

    def fill_Header(self, title=''):
        """ """
        headdict = self.header.copy()
        self.report.dict_to_sheet(headdict, 'Header', title=title)

    def fill_Meta(self):
        """ """
        metadict = self.meta.copy()
        self.report.dict_to_sheet(metadict, 'Meta', title='MetaInfo')

    def fill_Sheet(self, sheet):
        """ """
        df = self.data[sheet]
        self.report.df_to_sheet(df, sheet, index=True, header=True)

    def fill_allDataSheets(self):
        """ """
        for sheet in list(self.data.keys()):
            self.fill_Sheet(sheet)

    def fill_Figures(self):
        figsdict = self.figs
        if len(figsdict) == 0:
            return

        figskeys = figsdict['keys']
        jump = figsdict['jump']

        for ik, key in enumerate(figskeys):
            figpath = figsdict[key]
            cellnum = 1 + ik * jump
            self.report.add_image(figpath, 'figs', cell='A%i' % cellnum)

    def init_wb_and_fillAll(self, header_title=''):
        self.init_workbook()
        self.fill_Header(title=header_title)
        self.fill_Meta()
        self.fill_allDataSheets()
        self.fill_Figures()

    def get_textable(self, sheet, caption='', fitwidth=False, tiny=False,
                     **kwargs):
        """ """
        _kwargs = dict(multicolumn=True, multirow=True, longtable=True, index=False)
        _kwargs.update(kwargs)
        if 'columns' in _kwargs:
            data = self.data[sheet][:,_kwargs['columns']].copy()
        else:
            data = self.data[sheet].copy()

        tex = data.to_latex(**_kwargs)
        
        #if 'columns' not in _kwargs:
        #    ncols = len(self.data[sheet].columns)
        #else:
        ncols = len(_kwargs['columns'])

        if _kwargs['index']:
            ncols += 1

        wrapped = wraptextable(tex, ncols, caption, fitwidth, tiny, 
            longtable=_kwargs['longtable'])

        return wrapped


    def savehardcopy(self, filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path, '%s.xlsx' % self.rootname)
        if os.path.exists(filef):
            os.system('rm %s' % filef)
        self.report.save(filef)


class FitsTables_CDP(CDP):
    """Fits Table CDP."""

    formatsdict = dict(
            int16='I',
            int32='J',
            int64='K',
            char='A',
            float32='E',
            float64='D')

    def __init__(self, *args, **kwargs):
        """ """
        super(FitsTables_CDP, self).__init__(*args, **kwargs)
        self.data = OrderedDict()
        self.hdulist = None

    def ingest_inputs(self, data, meta=None, header=None, figs=None):
        """ """
        if meta is None:
            meta = OrderedDict()
        if header is None:
            header = OrderedDict()

        self.header.update(header)
        self.meta.update(meta)
        self.data.update(data)

    def init_HDUList(self):
        """ """
        self.hdulist = fts.HDUList()
        self.hdulist.append(fts.PrimaryHDU())
        self.hdulist.append(fts.TableHDU(name='META'))


    def fill_Header(self):
        """ """
        self.hdulist[0].header.update(self.header)
        

    def fill_Meta(self):
        """ """
        self.hdulist[1].header.update(self.meta)

    def fill_Table(self, sheet):
        """ """
        dd = self.data[sheet]

        columns = []
        for k in dd.keys():
            dtype = dd[k].dtype.name

            if 'string' in dtype:
                kformat = '%s%i' % (self.formatsdict['char'],
                    dd[k].dtype.itemsize)
            else:
                kformat = self.formatsdict[dtype]

            columns.append(fts.Column(name=k,
                format=kformat,
                array=dd[k].copy()))
            if kformat=='A':stop()

        self.hdulist.append(fts.BinTableHDU.from_columns(columns,
            name=sheet))

    def fill_allTables(self):
        """ """
        for sheet in self.data.keys():
            self.fill_Table(sheet)

    def init_HL_and_fillAll(self, *args):

        self.init_HDUList()
        self.fill_Header()
        self.fill_Meta()
        self.fill_allTables()

    def savehardcopy(self, filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path, '%s.fits' % self.rootname)
        self.hdulist.writeto(filef,overwrite=True)

    def checkInputHardcopy(self, inhdulist):
        """To be filled up by a child-class, if so desired."""
        pass

    def loadhardcopy(self, filef):
        """ """
        inhdulist = fts.open(filef)
        self.checkInputHardcopy(inhdulist)
        self.hdulist = inhdulist
    



class CCD_CDP(CDP):
    """CCD Calibration Data Product"""

    def __init__(self, *args, **kwargs):
        """ """
        super(CCD_CDP, self).__init__(*args, **kwargs)
        self.data = OrderedDict()
        self.ccdobj = ccdmod.CCD()

    def ingest_inputs(self, data, meta=None, header=None):
        """ """
        dmeta = OrderedDict()
        dmeta['ID'] = self.ID
        dmeta['BLOCKID'] = self.BLOCKID
        dmeta['CHAMBER'] = self.CHAMBER

        if meta is None:
            dmeta = dmeta.copy()
        else:
            dmeta.update(meta)
        if header is None:
            header = OrderedDict()

        self.meta.update(dmeta)
        self.header.update(header)
        self.data.update(data)

        if len(self.meta) > 0:
            self.ccdobj.add_extension(data=None,
                                      headerdict=self.meta,
                                      label='META')

        labels = self.data['labels']

        for i, label in enumerate(labels):

            if i == 0:
                self.ccdobj.add_extension(data=self.data[label],
                                          label=label,
                                          headerdict=self.header)
            else:
                self.ccdobj.add_extension(data=self.data[label],
                                          label=label)

    def savehardcopy(self, filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path, '%s.fits' % self.rootname)
        if os.path.exists(filef):
            os.system('rm %s' % filef)
        self.ccdobj.writeto(filef)

class LE1_CDP(CDP):
    """LE1 FPA Image CDP. One extension per Quadrant."""

    def __init__(self, *args, **kwargs):
        super(LE1_CDP, self).__init__(*args, **kwargs)
        self.fpaobj = fpa_dm.FPA_LE1()


    def ingest_inputs(self, data, header=None, inextension=-1, fillval=0):
        """ """

        if header is None:
            header = OrderedDict()

        self.header.update(header)

        self.fpaobj.fillval=fillval
        self.fpaobj.initialise_as_blank()

        self.fpaobj.set_extension(iext=0, data=None, header=None, headerdict=self.header, label=None)

        # Get the CCDs from data and fill up self.fpaobj

        NCOLS_FPA = self.fpaobj.fpamodel.NCOLS
        NROWS_FPA = self.fpaobj.fpamodel.NSLICES

        for jY in range(NCOLS_FPA):
            for iX in range(NROWS_FPA):
                Ckey = 'C_%i%i' % (jY + 1, iX + 1)
                ccdobj = copy.deepcopy(data[Ckey])

                self.fpaobj.set_ccdobj(ccdobj, Ckey, inextension=inextension)


    def savehardcopy(self, filef='', clobber=True, uint16=False):
        """ """
        if filef == '':
            filef = os.path.join(self.path, '%s.fits' % self.rootname)
        if os.path.exists(filef):
            os.system('rm %s' % filef)
        self.fpaobj.savetoFITS(filef,clobber=clobber,unsigned16bit=uint16)





class Json_CDP(CDP):
    """Generic Json Object CDP."""

    def __init__(self,*args,**kwargs):
        """ """
        super(Json_CDP, self).__init__(*args, **kwargs)

        self.data = OrderedDict()

    def ingest_inputs(self, data, meta=None, header=None):
        """ """
        if meta is None:
            meta = OrderedDict()
        if header is None:
            header = OrderedDict()

        self.header.update(header)
        self.meta.update(meta)
        self.data.update(data)


    def savehardcopy(self, filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path, self.rootname)
        if os.path.exists(filef):
            os.system('rm %s' % filef)

        outdict = OrderedDict()
        outdict['header'] = self.header.copy()
        outdict['meta'] = self.meta.copy()
        outdict['data'] = self.data.copy()

        vjson.save_jsonfile(outdict, filef)

    def loadhardcopy(self, filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path, self.rootname)

        outdict = vjson.load_jsonfile(filef, useyaml=True)

        self.header = outdict['header'].copy()
        self.meta = outdict['meta'].copy()
        self.data = outdict['data'].copy()




