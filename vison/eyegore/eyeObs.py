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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.animation as animation
from matplotlib import pyplot as plt


from pdb import set_trace as stop
from optparse import OptionParser
import sys, os
import numpy as np
import time
import string as st
import datetime

from vison.datamodel import HKtools 
from vison.datamodel import EXPLOGtools as ELtools
from vison import data as vdata
#import loadHK_QFM,allHK_keys

#from multiprocessing.dummy import Pool

import Tkinter as tk
import ttk
import tkFont as tkFont
from PIL import Image, ImageTk

from eyeHK import HKDisplay
from eyeCCDs import ImageDisplay
# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)




class ExpLogDisplay(tk.Toplevel):
    """ """
    
    def __init__(self,parent,path,interval,elvis='6.3.0'):
        
        self.path = path
        self.elvis = elvis
        self.date = '21-02-80'
        self.interval = interval
        self.explogf = None
        self.EXPLOG = dict()
        self.ix0 = 1
        self.nEL = 1
        
        
        self.labels = {}
        self.elementHeader = []
        self.elementList = []
        self.info = ""
        self.tree = None
        
        tk.Toplevel.__init__(self,parent)
        
        self.minsize(width=850,height=500)
        
        #self.wm_title('EXP-LOG')
        
        self.info = """\
        Click on header to sort by that column. To change width of column drag boundary.
        """
        
        self.wm_title('EXP-LOG')
        
        self.frame = tk.Frame(self)
        self.frame.pack(fill='both', expand=False)
        
        #msg = ttk.Label(self,wraplength="4i", justify="left", anchor="n",
        #    padding=(10, 2, 10, 6), text=self.info)
        #msg.grid(row=0,column=0,columnspan=2,in_=self.frame)        
        
        self.labels['NObsID'] = tk.Label(self, text="NObs = 0", font=LARGE_FONT)
        self.labels['NObsID'].grid(row=0,column=0,in_=self.frame,sticky='nsew')
        self.labels['NEntries'] = tk.Label(self, text="NEntries = 0", font=LARGE_FONT)
        self.labels['NEntries'].grid(row=0,column=1,in_=self.frame)
        
        self.search_EXPLOG()
        self.get_data()        
        
        self.elementHeader = self.EXPLOG.colnames
        
        self.build_elementList()
        self.createScrollableTreeview()

        self.frame.grid_columnconfigure(0, weight=1)
        self.frame.grid_rowconfigure(0, weight=1)
        
        self.update()
    
    def update(self):
        
        self.search_EXPLOG()
        self.get_data()      
        self.build_elementList()
        self.buildTree()
        
        self.ix0 += 1
        #print self.ix0
        self.tree.yview_scroll(1,'units')
        
        self.after(self.interval,self.update)


#    def get_updater(self,interval):
#        
#        def update(self):
#            self.search_EXPLOG()
#            self.get_data()      
#            self.build_elementList()
#            self.buildTree()
#            
#            self.ix0 += 1
#            #print self.ix0
#            self.tree.yview_scroll(1,'units')
#            
#            self.after(5000,update)
#
#        return update
    
    def build_elementList(self):
        """ """
    
        elementList = []
        
        for ix in range(len(self.EXPLOG)):            
            row = []
            for jx,colname in enumerate(self.elementHeader):
                row.append(self.EXPLOG[colname][ix])
            elementList.append(tuple(row))
        
        
        self.elementList = elementList    
        
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
        
        with open(self.explogf) as f:
            nEL = len(f.readlines())
    
        #if nEL <= self.nEL:   # COMMENTED ON TESTS
        #    return
    
        self.nEL = nEL
    
        EXPLOG = ELtools.loadExpLog(self.explogf,elvis=self.elvis)
        
        self.EXPLOG = EXPLOG[self.ix0:self.ix0+100]


    def createScrollableTreeview(self):

        # create a treeview with dual scrollbars
        self.subframe = tk.Frame(self.frame) 
        #self.subframe.pack(fill='both', expand=True)
        self.subframe.grid(row=1,column=0, columnspan=2)
        self.tree = ttk.Treeview(self,columns=self.elementHeader, show="headings")
        vsb = ttk.Scrollbar(self,orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self,orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(column=0, row=0, sticky='nsew', in_=self.subframe)
        vsb.grid(column=1, row=0, sticky='ns', in_=self.subframe)
        hsb.grid(column=0, row=1, sticky='ew', in_=self.subframe)
        self.subframe.grid_columnconfigure(0, weight=1)
        self.subframe.grid_rowconfigure(0, weight=1)


    def buildTree(self):
        
        # add possibility to sort by each column
        for col in self.elementHeader:
            self.tree.heading(col, text=col,
                command=lambda c=col: self.sortBy(self.tree, c, 0))
            # adjust the column's width to the header string
            self.tree.column(col, width=int(tkFont.Font(font='Helvetica',size=12).measure(col)*1.5))
            

        for i,item in enumerate(self.elementList):
            if i % 2 == 0: parity = 'pair'
            elif i % 2 !=0: parity = 'odd'
            self.tree.insert('', 'end', values=item, tags=(parity,))

            # adjust column's width if necessary to fit each value
#            for ix, val in enumerate(item):
#                col_w = tkFont.Font().measure(val)
#                if self.tree.column(self.elementHeader[ix], width=None) < col_w:
#                    self.tree.column(self.elementHeader[ix], width=col_w)

        self.tree.tag_configure('pair',background='#B6D2D2')
        self.tree.tag_configure('odd',background='#AFE0B5')
        
        self.tree.bind("<Double-1>",self.OnDoubleClick)
        
        
    
    def OnDoubleClick(self,event):
        #item = self.tree.selection()[0]
        item = self.tree.identify('item',event.x,event.y)
        print "you clicked on", self.tree.item(item,"value")
    
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



    