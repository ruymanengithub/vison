#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Euclid Ground Calibration Campaign


Classes to store Calibration Data Products.

Created on Tue Feb 27 10:58:42 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import copy
import os
import string as st

from vison import __version__
from vison.support.files import cPickleDumpDictionary, cPickleRead
from vison.datamodel import ccd as ccdmod
from vison.support.excel import ReportXL
# END IMPORT


def loadCDPfromPickle(pickf):
    """ """
    cdp = cPickleRead(pickf)
    return cdp


class CDP(object):
    """ """

    rootname = 'Unknown'
    path = ''
    header = OrderedDict()
    meta = OrderedDict()
    data = None

    def __init__(self, *args, **kwargs):
        """ """

        defaults = dict(rootname='Uknown', path='', header=OrderedDict(),
                        meta=OrderedDict(), data=None)
        defaults.update(args)
        defaults.update(kwargs)

        for key in defaults.keys():
            self.__setattr__(key, defaults[key])

    def savetopickle(self, pickf=''):
        """ """
        if pickf == '':
            pickf = os.path.join(self.path, '%s.pick' % self.rootname)

        outdict = copy.deepcopy(self.__dict__)
        cPickleDumpDictionary(outdict, pickf)

    def loadfrompickle(self, pickf=''):
        """ """
        if pickf == '':
            pickf = os.path.join(self.path, '%s.pick' % self.rootname)
        saved = cPickleRead(pickf)
        self.__dict__.update(saved)

    def savehardcopy(self, filef=''):
        """ """
        raise NotImplementedError(
            'Subclass implements abstract method (if needed).')


class Tables_CDP(CDP):
    """ """

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

        for i, k in enumerate(self.data.iterkeys()):
            self.report.wb.create_sheet(k, i+2)

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
        for sheet in self.data.keys():
            self.fill_Sheet(sheet)
            
    def init_wb_and_fillAll(self, header_title=''):
        self.init_workbook()
        self.fill_Header(title=header_title)
        self.fill_Meta()
        self.fill_allDataSheets()

    def get_textable(self, sheet, caption=''):
        """ """
        tex = self.data[sheet].to_latex(
            multicolumn=True, multirow=True, longtable=True, index=True)
        tex = st.split(tex, '\n')
        if caption != '':
            tex.insert(-2, r'\caption{%s}' % caption)
        return tex


    def savehardcopy(self, filef=''):
        """ """
        if filef == '':
            filef = os.path.join(self.path, '%s.xlsx' % self.rootname)
        if os.path.exists(filef):
            os.system('rm %s' % filef)
        self.report.save(filef)


class CCD_CDP(CDP):

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
        
        if len(self.meta)>0:
            self.ccdobj.add_extension(data=None,
                                      headerdict=self.meta,
                                      label='META')
        
        labels = self.data['labels']
        
        for i, label in enumerate(labels):
            
            if i==0:
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
    