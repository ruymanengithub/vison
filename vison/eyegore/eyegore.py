# -*- coding: utf-8 -*-
"""

eyegore

data acquisition monitoring script for vison package.


'- You must be Igor...
- No, it's pronounced "Eye-gore".'

Created on Thu Feb  2 15:27:39 2017

:Author: Ruyman Azzollini
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
#import loadHK_QFM,allHK_keys

#from multiprocessing.dummy import Pool

import Tkinter as tk
import ttk
import tkFont as tkFont

# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)


class ImageDisplay(tk.Toplevel):
    """ """
    
    def __init__(self,parent,path):
        """ """
        tk.Toplevel.__init__(self,parent)
        self.parent = parent
        self.path = path
        self.wm_title('Image Display')
        
        self.minsize(width=850,height=400)
        
        
        f1 = plt.figure(figsize=(8,4),dpi=100)
        ax1 = f1.add_subplot(111)
        
        self.f = f1
        self.axs = [ax1]
        
        canvas = FigureCanvasTkAgg(f1, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        #self.gen_render()(0)
        
        #self.start_updating()
        #self.mainloop()

    
    def get_data(self):
        """PENDING"""
        t=datetime.datetime.now().second
        #print t
        yield t
    
    
    def gen_render(self):
        """ """
        def render(image):
            
            t = datetime.datetime.now()
            s = t.strftime('%H:%M:%S')
            
            #ax = f.add_subplot(111)
            x = np.arange(10)
            y = x**2.+image

            self.axs[0].clear()
            self.axs[0].plot(x,y)
            self.axs[0].set_title('%s' % s)
            self.axs[0].set_ylim([0,140])

        return render
    
    def start_updating(self,interval=5000):
        #self.do_update = True
        
        f = self.f
        render = self.gen_render()
        #while self.do_update:
        return animation.FuncAnimation(f, render,self.get_data,interval=interval)
        
    
    def stop_updating(self):
        self.do_update = False
    
    
    def get_animator(self):
        f = self.f
        render = self.gen_render()
        

class HKDisplay(tk.Toplevel):
    """ """
    def __init__(self,parent,path,elvis='6.3.0'):
        
        self.elvis = elvis
        self.date = '21-02-80'
        self.HKfile = None
        self.HK = dict()
        self.page = 1
        self.path = path
        self.nHKlines = 1
        
        self.HKkeys = HKtools.allHK_keys[elvis]
        
        self.search_HKfile()
        
        tk.Toplevel.__init__(self,parent)
        
        self.minsize(width=700,height=925)
        
        self.wm_title('HK Display')
        
        #label = tk.Label(self, text="HK Display", font=LARGE_FONT)
        #label.pack(side=tk.TOP,pady=10,padx=10)
        #label.grid(row=0,pady=10,padx=10,columnspan=2)
        
        frame = tk.Frame(self)
        l1 = tk.Label(frame,text="Page: ",font=LARGE_FONT)
        l1.pack(side="left")
        
        pvar = tk.StringVar() 
        pvar.set('1')
        pentry = tk.Entry(frame,textvariable=pvar)
        
        pentry.pack(side="left")
        frame.pack(side="top")
        
        def update_page(event):
            try:
                self.page = int(pvar.get())
            except ValueError:
                self.page = 1
        
        pentry.bind('<Return>',update_page)
        
        f = plt.figure(figsize=(4,8),dpi=100)
        
        axs = []
        
        for i in range(9):
            
            ax = f.add_subplot(3,3,i+1)            
            axs.append(ax)
        
        self.f = f
        self.axs = axs
        #self.f.autofmt_xdate()
        
        canvas = FigureCanvasTkAgg(self.f, self)
        plt.tight_layout(rect=[0, 0, 1, 1])
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        #cw = canvas.get_tk_widget()
        #cw.grid(row=3,sticky=tk.S,columnspan=2)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        #canvas._tkcanvas.grid(row=3,sticky=tk.S,columnspan=2)
    
    
    def search_HKfile(self):
        """ """
        
        struct_date = time.strptime(os.path.split(st.replace(self.path,'_',''))[-1],'%d%b%y')
        date_infile = time.strftime('%d-%m-%y',struct_date)
        
        self.date = date_infile

        tmp_HK = 'HK_%s_ROE1.txt' % date_infile
        tmp_HK = os.path.join(self.path,tmp_HK)
        
        isthere = os.path.exists(tmp_HK)
        if isthere:
            self.HKfile = tmp_HK
        else:
            print 'HKfile %s not found' % tmp_HK
            self.HKfile = None
    
    def select_HKkeys(self):
        """ """
        
        allHKkeys = self.HKkeys
        nkeys = len(allHKkeys)-1
        page = self.page
        #print page
        
        if page * 9 > nkeys-1:
            return
        
        ix0 = (page-1)*9
        ix1 = ix0 + 9 if ix0+9<=nkeys else None
        
        HKkeys_to_plot = allHKkeys[1:][ix0:ix1]
        
        self.HKkeys_to_plot = HKkeys_to_plot
        
    
    def get_data(self):
        """ """
        
        self.select_HKkeys()
        
        if self.HKfile is None:
            yield self.HK
        
        with open(self.HKfile) as f:
            nHKlines = len(f.readlines())
        
        if nHKlines <= self.nHKlines:
            yield self.HK
        
        self.nHKlines = nHKlines
        
        HK = HKtools.loadHK_QFM(self.HKfile,elvis=self.elvis)
        
        dtobjarr = np.array([datetime.datetime.strptime('%s_%s'  % (self.date,item),'%d-%m-%y_%H:%M:%S') \
                             for item in HK['TimeStamp']])
        
        pHK = dict(time = dtobjarr)
        
        subKeys = [Key for Key in HK.keys() if Key != 'TimeStamp']
        
        for key in subKeys:
            pHK[key] = HK[key].copy()
        
        self.HK = pHK
        
        yield pHK
        
    
    def gen_render(self):
        
        def render(pHK):
            
            nHK = len(self.HKkeys_to_plot)
            
            for i in range(nHK):
                
                ax = self.axs[i]
                
                HKname = self.HKkeys_to_plot[i]
                
                try:
                    x = pHK['time'].copy()
                    y = pHK[HKname].copy()
                except KeyError:
                    x = np.array([0,1])
                    y = np.array([0,1])
                                
                ax.clear()
                
                if np.any(np.isnan(y)):
                    yp = y.copy()
                    yp[np.isnan(y)] = 0
                    ax.plot(x,yp)
                    ax.plot(x[np.isnan(y)],yp[np.isnan(y)],'ro-')                    
                else:
                    ax.plot(x,y)
                
                ax.set_title(HKname)
                
                try:
                    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +\
                        ax.get_xticklabels() + ax.get_yticklabels()):
                        item.set_fontsize(10)
                except:
                    pass
            
            #plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(HKtools.format_date)) 
            try: self.f.autofmt_xdate()
            except:
                pass
            #plt.tight_layout(rect=[0, 0, 1, 1])

        return render
        

    def start_updating(self,interval):
        #self.do_update = True
        
        f = self.f
        render = self.gen_render()
        #while self.do_update:
        return animation.FuncAnimation(f, render,self.get_data,interval=interval)



class ExpLogDisplay(tk.Toplevel):
    """ """
    
    def __init__(self,parent,path,elvis='6.3.0'):
        
        self.elvis = elvis
        self.date = '21-02-80'
        self.path = path
        self.explogf = None
        self.EXPLOG = dict()
        self.nEL = 1

        self.elementHeader = []
        self.elementList = []
        self.info = ""
        self.tree = None

        self.search_EXPLOG()
        self.get_data()
        
        tk.Toplevel.__init__(self,parent)
        
        self.minsize(width=850,height=500)
        
        self.wm_title('EXP-LOG')
        
        self.info = """\
        Click on header to sort by that column.
        To change width of column drag boundary.
        """
        
        elementHeader = self.EXPLOG.colnames
        
        elementList = []
        
        for ix in range(len(self.EXPLOG)):            
            row = []
            for jx,colname in enumerate(elementHeader):
                row.append(self.EXPLOG[colname][ix])
            elementList.append(tuple(row))
        
        self.elementList = elementList
        self.elementHeader = elementHeader
        
        self.setupWidgets()
        #stop()
        self.buildTree()
        
        
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
            return self.EXPLOG
        
        with open(self.explogf) as f:
            nEL = len(f.readlines())
    
        if nEL <= self.nEL:
            return  self.EXPLOG
    
        self.nEL = nEL
    
        EXPLOG = ELtools.loadExpLog(self.explogf,elvis=self.elvis)
        
        self.EXPLOG = EXPLOG
        
        return self.EXPLOG


    def setupWidgets(self):
        
        msg = ttk.Label(self,wraplength="4i", justify="left", anchor="n",
            padding=(10, 2, 10, 6), text=self.info)
        msg.pack(fill='x')
        
        frame = ttk.Frame(self)
        frame.pack(fill='both', expand=True)

        # create a treeview with dual scrollbars
        self.tree = ttk.Treeview(self,columns=self.elementHeader, show="headings")
        vsb = ttk.Scrollbar(self,orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self,orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.grid(column=0, row=0, sticky='nsew', in_=frame)
        vsb.grid(column=1, row=0, sticky='ns', in_=frame)
        hsb.grid(column=0, row=1, sticky='ew', in_=frame)

        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)

    def buildTree(self):

        for col in self.elementHeader:
            self.tree.heading(col, text=col.title(),
                command=lambda c=col: self.sortBy(self.tree, c, 0))
            # adjust the column's width to the header string
            self.tree.column(col, width=tkFont.Font().measure(col.title()))

        for i,item in enumerate(self.elementList):
            if i % 2 == 0: parity = 'pair'
            elif i % 2 !=0: parity = 'odd'
            self.tree.insert('', 'end', values=item, tags=(parity,))

            # adjust column's width if necessary to fit each value
            for ix, val in enumerate(item):
                col_w = tkFont.Font().measure(val)
                if self.tree.column(self.elementHeader[ix], width=None) < col_w:
                    self.tree.column(self.elementHeader[ix], width=col_w)

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
   

class Eyegore(tk.Tk):
    """ """
    
    def __init__(self,path,broadcast,interval=3000):
        """ """
        tk.Tk.__init__(self)
        
        if path[-1] == os.path.sep:
            path = path[:-1]
        self.path = path
        self.interval = interval
        self.broadcast = broadcast
        
        self.run()
        
                
    def run(self):
        

        self.withdraw()
        
        Ds = dict(image=ImageDisplay,hk=HKDisplay,explog=ExpLogDisplay)
        dkeys = ['image','hk','explog']
        
        display1 = Ds[dkeys[0]](self,self.path)        
        ani = display1.start_updating(self.interval)
       
        display2 = Ds[dkeys[1]](self,self.path)
        ani = display2.start_updating(self.interval)
            
        display3 = Ds[dkeys[2]](self,self.path)
        ani = display3.start_updating(self.interval)
        
        self.mainloop()


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-p","--path",dest="path",default='',help="day-path to be monitored.")
    parser.add_option("-B","--broadcast",dest="broadcast",action='store_true',help="")
    
    (options, args) = parser.parse_args()
    
    if options.path == '':
        parser.print_help()
        sys.exit()
        
    path = options.path
    broadcast = options.broadcast
    
    if not os.path.exists(path):
        sys.exit('HKmonitory.py: %s does not exist' % path)
    else:
        
        header = '\n\n#########################################################\n'+\
                 '#                                                       #\n'+\
                 '#                                                       #\n'+\
                 '#       starting EYEGORE on path:                       #\n'+\
                 '#            %s                                         \n'+\
                 '#                                                       #\n'+\
                 '#                                                       #\n'+\
                 '#########################################################\n'
        
        print header % path
    
    
    app = Eyegore(path,broadcast=broadcast)

    
    
