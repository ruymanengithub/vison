#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:22:36 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import matplotlib
matplotlib.use("TkAgg")
import glob
from astropy.table import Table

import os
import time
import string as st
import numpy as np

from vison.datamodel import EXPLOGtools as ELtools
from vison.support import context
from vison.eyegore.eyelib import get_bitsize

import Tkinter as tk
import ttk
import tkFont as tkFont

try:
    import pyds9
except ImportError:
    print 'pyds9 not installed in the system. Will not be possible to liaise with SAO-ds9.'


# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)
medium_font = ("Verdana", 10)


def isNumeric(s):
    """
    test if a string s is numeric
    """
    for c in s:
        if c in "1234567890-.+":
            numeric = True
        else:
            return False
    return numeric

def changeNumeric(data):
    """
    if the data to be sorted is numeric change to float
    """
    new_data = []
    if isNumeric(data[0][0]):
        # change child to a float
        for child, col in data:
            new_data.append((float(child), col))
        return new_data
    return data

class ExpLogDisplay(tk.Toplevel):
    """ """

    def __init__(self, parent, path, interval, elvis=context.elvis):

        self.path = path
        self.elvis = elvis
        self.date = '21-02-80'
        self.interval = interval
        self.explogfs = None
        self.EXPLOG = Table()
        self.nEL = 0  # Nr. lines in EXPLOG
        self.sEL = 0  # size of EXPLOG, bytes
        self.log = parent.log

        self.labels = {}
        self.elementHeader = []
        self.elementList = []
        self.info = ""
        self.tree = None
        #self.NObsID = 0
        #self.Nentries = 0
        self.updating = True

        tk.Toplevel.__init__(self, parent)

        # self.minsize(width=850,height=500)

        self.info = """\
        Click on header to sort by that column. To change width of column drag boundary.
        """

        self.wm_title('EXP-LOG')

        self.fr0 = tk.Frame(self)
        self.fr0.pack(fill='both', expand=False)
        # self.fr0.grid()

        # msg = ttk.Label(self,wraplength="4i", justify="left", anchor="n",
        #    padding=(10, 2, 10, 6), text=self.info)
        # msg.grid(row=0,column=0,columnspan=2,in_=self.frame)

        self.fr1 = tk.Frame(self)
        self.fr1.grid(row=0, in_=self.fr0)

        self.labels['NObsIDs'] = dict()
        self.labels['NObsIDs']['var'] = tk.StringVar()
        self.labels['NObsIDs']['var'].set("NObs = 0")
        self.labels['NObsIDs']['app'] = tk.Label(self, textvariable=self.labels['NObsIDs']['var'],
                                                 font=medium_font, bd=2, relief='groove')
        self.labels['NObsIDs']['app'].grid(
            row=0, column=0, in_=self.fr1, sticky='w')
        self.labels['NEntries'] = dict()
        self.labels['NEntries']['var'] = tk.StringVar()
        self.labels['NEntries']['var'].set("NEntries = 0")
        self.labels['NEntries']['app'] = tk.Label(self, textvariable=self.labels['NEntries']['var'],
                                                  font=medium_font, bd=2, relief='groove')
        self.labels['NEntries']['app'].grid(
            row=0, column=1, in_=self.fr1, sticky='w')

        self.fr1.grid_columnconfigure(0, weight=1)
        self.fr1.grid_rowconfigure(0, weight=1)

        self.search_EXPLOGs()
        self.get_data()

        self.elementHeader = self.EXPLOG.colnames

        self.build_elementList()
        self.createScrollableTreeview()

        self.fr0.grid_columnconfigure(0, weight=1)
        self.fr0.grid_rowconfigure(0, weight=1)

        self.update()

    def update(self):

        if self.updating:

            self.tree.yview_moveto(1)

            self.search_EXPLOGs()
            self.get_data()

            self.build_elementList()

            self.buildTree()

            self.tree.yview_moveto(1)

        self.after(self.interval, self.update)

    def createScrollableTreeview(self):

        # create a treeview with dual scrollbars
        self.fr2 = tk.Frame(self)
        #self.subframe.pack(fill='both', expand=True)
        self.fr2.grid(row=1, column=0, in_=self.fr0)
        self.tree = ttk.Treeview(
            self, columns=self.elementHeader, show="headings")
        for col in self.elementHeader:
            self.tree.column(col, minwidth=100, width=100, stretch=True)
        vsb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        self.tree.grid(column=0, row=0, sticky='nsew', in_=self.fr2)
        vsb.grid(column=1, row=0, sticky='ns', in_=self.fr2)
        hsb.grid(column=0, row=1, sticky='ew', in_=self.fr2)

        # self.tree.yview_scroll(-1,'units')
        # self.tree.yview_moveto(0.5)

        self.tree['height'] = 30

        self.fr2.grid_columnconfigure(0, weight=1)
        self.fr2.grid_rowconfigure(0, weight=1)

    def build_elementList(self):
        """ """

        for ix in range(len(self.EXPLOG)):
            row = []
            for jx, colname in enumerate(self.elementHeader):
                row.append(self.EXPLOG[colname][ix])
            self.elementList.append(tuple(row))

    def search_EXPLOGs(self):
        """ """

        struct_date = time.strptime(os.path.split(
            st.replace(self.path, '_', ''))[-1], '%d%b%y')
        date_infile = time.strftime('%d%m%y', struct_date)

        self.date = date_infile

        tmp_EL = 'EXP_LOG_*.txt'
        tmp_EL = os.path.join(self.path, tmp_EL)

        explogfs = glob.glob(tmp_EL)

        arethere = len(explogfs) > 0

        if arethere:
            self.explogfs = explogfs
        else:
            print 'EXPLOGs %s not found' % tmp_EL
            self.explogfs = None


    def loadExplogs(self):
        """ """
        ELList = []
        for item in self.explogfs:
            ELList.append(ELtools.loadExpLog(item, elvis=self.elvis))
        if len(ELList) < 2:
            return ELList[0]
        else:
            return ELtools.mergeExpLogs(ELList, addpedigree=False)

    def get_data(self):
        """ """

        if self.explogfs is None:
            return

        sEL = get_bitsize(self.explogfs)

        if sEL <= self.sEL:
            return

        EXPLOG = self.loadExplogs()
        #EXPLOG = ELtools.loadExpLog(self.explogf,elvis=self.elvis)

        NObsIDs = len(np.unique(EXPLOG['ObsID']))

        nEL = len(EXPLOG)  # commented on TESTS
        # nEL = self.nEL + 5 # TESTS
        self.EXPLOG = EXPLOG[self.nEL:nEL]

        self.nEL = nEL
        self.sEL = sEL

        self.labels['NEntries']['var'].set('Nentries = %i' % len(EXPLOG))
        self.labels['NObsIDs']['var'].set('NObsIDs = %i' % NObsIDs)

    def growTree(self):

        for i, item in enumerate(self.elementList):
            #print i+self.nEL

            if (i+self.nEL) % 2 == 0:
                parity = 'odd'
            elif (i+self.nEL) % 2 != 0:
                parity = 'pair'
            self.tree.insert('', 'end', values=item, tags=(parity,))

            # adjust column's width if necessary to fit each value
#            for ix, val in enumerate(item):
#                col_w = tkFont.Font().measure(val)
#                if self.tree.column(self.elementHeader[ix], width=None) < col_w:
#                    self.tree.column(self.elementHeader[ix], width=col_w)

        self.elementList = []

    def buildTree(self):

        # add possibility to sort by each column
        for col in self.elementHeader:
            self.tree.heading(col, text=col,
                              command=lambda c=col: self.sortBy(self.tree, c, 0))
            # adjust the column's width to the header string
            self.tree.column(col, width=int(tkFont.Font(
                font='Helvetica', size=12).measure(col)*1.5))

        self.growTree()

        self.tree.yview_moveto(1)

        # for i in range(1,2*growth+1):
        #    self.tree.yview_scroll(1,'units')

        self.tree.tag_configure('pair', background='#B6D2D2')
        self.tree.tag_configure('odd', background='#AFE0B5')

        self.tree.bind("<Double-1>", self.OnDoubleClick)

    def OnDoubleClick(self, event):
        #item = self.tree.selection()[0]
        item = self.tree.identify('item', event.x, event.y)
        values = self.tree.item(item, "value")
        ObsID = int(values[0])
        try:
            d = pyds9.DS9()
        except NameError:
            print "pyds9 not installed, can't open DS9!"
            return

        for CCD in [1, 2, 3]:
            tmpfits = os.path.join(
                self.path, 'EUC_%i_*D*T_ROE1_CCD%i.fits' % (ObsID, CCD))
            try:
                iFITSf = glob.glob(tmpfits)[0]
            except IndexError:
                #print tmpfits
                iFITSf = None

            if iFITSf is not None:
                d.set("frame %i" % CCD)
                d.set("file %s" % iFITSf)
                d.set("zoom to fit")
                d.set("scale mode zscale")


    def sortBy(self, tree, col, descending):
        """
        sort tree contents when a column header is clicked
        """
        # grab values to sort
        data = [(tree.set(child, col), child)
                for child in tree.get_children('')]
        # if the data to be sorted is numeric change to float
        data = changeNumeric(data)
        # now sort the data in place
        data.sort(reverse=descending)
        for ix, item in enumerate(data):
            tree.move(item[1], '', ix)
        # switch the heading so that it will sort in the opposite direction
        tree.heading(col,
                     command=lambda col=col: self.sortBy(tree, col, int(not descending)))
