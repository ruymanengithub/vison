#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:22:36 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import matplotlib
matplotlib.use("TkAgg")
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#import matplotlib.animation as animation
#from matplotlib import pyplot as plt
import glob


#from pdb import set_trace as stop
#from optparse import OptionParser
#import sys 
import os
#import numpy as np
import time
import string as st
#import datetime

#from vison.datamodel import HKtools 
from vison.datamodel import EXPLOGtools as ELtools
#from vison import data as vdata
#from vison.pipe import lib as pilib
from vison.support import context
#import loadHK_QFM,allHK_keys

#from multiprocessing.dummy import Pool

import Tkinter as tk
import ttk
import tkFont as tkFont

import pyds9


# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)


class ExpLogDisplay(tk.Toplevel):
    """ """
    
    def __init__(self,parent,path,interval,elvis=context.elvis):
        
        self.path = path
        self.elvis = elvis
        self.date = '21-02-80'
        self.interval = interval
        self.explogf = None
        self.EXPLOG = dict()
        self.nEL = 0 # Nr. lines in EXPLOG
        self.sEL = 0 # size of EXPLOG, bytes
        
        self.labels = {}
        self.elementHeader = []
        self.elementList = []
        self.info = ""
        self.tree = None
        
        tk.Toplevel.__init__(self,parent)
        
        #self.minsize(width=850,height=500)
        
        self.info = """\
        Click on header to sort by that column. To change width of column drag boundary.
        """
        
        self.wm_title('EXP-LOG')
        
        self.fr0 = tk.Frame(self)
        self.fr0.pack(fill='both', expand=False)
        #self.fr0.grid()
        
        #msg = ttk.Label(self,wraplength="4i", justify="left", anchor="n",
        #    padding=(10, 2, 10, 6), text=self.info)
        #msg.grid(row=0,column=0,columnspan=2,in_=self.frame)        
        
        self.fr1 = tk.Frame(self)
        self.fr1.grid(row=0,in_=self.fr0)
                
        self.labels['NObsID'] = tk.Label(self, text="NObs = 0", font=small_font)
        self.labels['NObsID'].grid(row=0,column=0,in_=self.fr1,sticky='w')
        self.labels['NEntries'] = tk.Label(self, text="NEntries = 0", font=small_font)
        self.labels['NEntries'].grid(row=0,column=1,in_=self.fr1,sticky='w')
        
        self.fr1.grid_columnconfigure(0, weight=1)
        self.fr1.grid_rowconfigure(0, weight=1)

        
        self.search_EXPLOG()
        self.get_data()
        
        self.elementHeader = self.EXPLOG.colnames
        
        self.build_elementList()
        self.createScrollableTreeview()

        self.fr0.grid_columnconfigure(0, weight=1)
        self.fr0.grid_rowconfigure(0, weight=1)
        
        self.update()
    
    def update(self):
                
        self.tree.yview_moveto(1)
    
        self.get_data()
        self.build_elementList()
                
        self.buildTree()
        
        self.tree.yview_moveto(1)
        
        self.after(self.interval,self.update)


    def createScrollableTreeview(self):

        # create a treeview with dual scrollbars
        self.fr2 = tk.Frame(self) 
        #self.subframe.pack(fill='both', expand=True)
        self.fr2.grid(row=1,column=0, in_=self.fr0)
        self.tree = ttk.Treeview(self,columns=self.elementHeader, show="headings")
        vsb = ttk.Scrollbar(self,orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self,orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        self.tree.grid(column=0, row=0, sticky='nsew', in_=self.fr2)
        vsb.grid(column=1, row=0, sticky='ns', in_=self.fr2)
        hsb.grid(column=0, row=1, sticky='ew', in_=self.fr2)

        #self.tree.yview_scroll(-1,'units')
        #self.tree.yview_moveto(0.5)
        
        self.tree['height'] = 20
        
        self.fr2.grid_columnconfigure(0, weight=1)
        self.fr2.grid_rowconfigure(0, weight=1)

    
    def build_elementList(self):
        """ """

        for ix in range(len(self.EXPLOG)):            
            row = []
            for jx,colname in enumerate(self.elementHeader):
                row.append(self.EXPLOG[colname][ix])
            self.elementList.append(tuple(row))
        
    def search_EXPLOG(self):
        """ """
        
        struct_date = time.strptime(os.path.split(st.replace(self.path,'_',''))[-1],'%d%b%y')
        date_infile = time.strftime('%d%m%y',struct_date)
        
        self.date = date_infile
        
        tmp_EL = 'EXP_LOG_%s.txt' % date_infile
        tmp_EL = os.path.join(self.path,tmp_EL)
        
        isthere = os.path.exists(tmp_EL)
        if isthere:
            self.explogf = tmp_EL
        else:
            print 'EXPLOG %s not found' % tmp_EL
            self.explogf = None
        
        
    def get_data(self):
        """ """
        
        if self.explogf is None:
            return
        
        statinfo = os.stat(self.explogf)
        sEL = statinfo.st_size # commented on TESTS
        #sEL = self.sEL + 1 # TESTS
        
        if sEL <= self.sEL:
            return
    
        EXPLOG = ELtools.loadExpLog(self.explogf,elvis=self.elvis)
        nEL = len(EXPLOG) # commented on TESTS
        #nEL = self.nEL + 5 # TESTS
        self.EXPLOG = EXPLOG[self.nEL:nEL]
        
        self.nEL = nEL
        self.sEL = self.sEL
        
    
    def growTree(self):
        

        for i,item in enumerate(self.elementList):
            if (i+self.nEL) % 2 == 0: parity = 'pair'
            elif (i+self.nEL) % 2 !=0: parity = 'odd'
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
            self.tree.column(col, width=int(tkFont.Font(font='Helvetica',size=12).measure(col)*1.5))
        
        self.growTree()
        
        self.tree.yview_moveto(1)
        
        #for i in range(1,2*growth+1):
        #    self.tree.yview_scroll(1,'units')

        self.tree.tag_configure('pair',background='#B6D2D2')
        self.tree.tag_configure('odd',background='#AFE0B5')
        
        self.tree.bind("<Double-1>",self.OnDoubleClick)
        
        
    
    def OnDoubleClick(self,event):
        #item = self.tree.selection()[0]
        item = self.tree.identify('item',event.x,event.y)
        values = self.tree.item(item,"value")
        ObsID = int(values[0])
        
        d = pyds9.DS9()
        
        for CCD in [1,2,3]:
            tmpfits = os.path.join(self.path,'EUC_%i_*D*T_ROE1_CCD%i.fits' % (ObsID,CCD))
            try: iFITSf = glob.glob(tmpfits)[0]
            except IndexError: 
                print tmpfits
                iFITSf = None
            
            if iFITSf is not None:
                d.set("frame %i" % CCD)
                d.set("file %s" % iFITSf)
                d.set("zoom to fit")
                d.set("scale mode zscale")
        
    
    def isNumeric(self, s):
        """
        test if a string s is numeric
        """
        for c in s:
            if c in "1234567890-.+":
                numeric = True
            else:
                return False
        return numeric

    def changeNumeric(self, data):
        """
        if the data to be sorted is numeric change to float
        """
        new_data = []
        if self.isNumeric(data[0][0]):
            # change child to a float
            for child, col in data:
                new_data.append((float(child), col))
            return new_data
        return data

    def sortBy(self, tree, col, descending):
        """
        sort tree contents when a column header is clicked
        """
        # grab values to sort
        data = [(tree.set(child, col), child) for child in tree.get_children('')]
        # if the data to be sorted is numeric change to float
        data =  self.changeNumeric(data)
        # now sort the data in place
        data.sort(reverse=descending)
        for ix, item in enumerate(data):
            tree.move(item[1], '', ix)
        # switch the heading so that it will sort in the opposite direction
        tree.heading(col,
            command=lambda col=col: self.sortBy(tree, col, int(not descending)))



    