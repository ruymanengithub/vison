#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""


Created on Fri Oct 13 14:11:41 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.animation as animation
from matplotlib import pyplot as plt

import os
from pdb import set_trace as stop
import numpy as np
import time
import string as st
import datetime

import Tkinter as tk
import ttk
import tkFont as tkFont

from vison.support import context
#from vison.pipe import lib as pilib
from vison.datamodel import HKtools
# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)


def _ax_render_HK(ax,x,y,HKlims,HKkey):
    """ """
    
        
    max_xticks = 6
    xloc = plt.MaxNLocator(max_xticks)            
        
    ax.clear()
        
                                
    if np.any(np.isnan(y)):
        yp = y.copy()
        yp[np.isnan(y)] = 0
        ax.plot(x,yp)
        ax.plot(x[np.isnan(y)],yp[np.isnan(y)],'ro-')                    
    else:
        ax.plot(x,y)
        
    if len(HKlims) == 0:
        ax.axhline(y=HKlims[0],ls='--',lw=2,color='r')
    elif len(HKlims) == 2:
        for ik in range(2):
            ax.axhline(y=HKlims[ik],ls='--',lw=2,color='r')
        
    HKtitle = '$%s$' % HKkey.replace('_','\_')
    ax.set_title(HKtitle)
    
    try:
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +\
            ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(10)
    except:
        pass
    
    ax.xaxis.set_major_locator(xloc)
    
    return ax
        



class SingleHKplot(tk.Toplevel):
    """ """
    
    def __init__(self,root):
        """ """
        
        tk.Toplevel.__init__(self,root)
        
        self.wm_title('HK Parameter')
        
        self.minsize(width=450,height=400)
        
        
        #frame = tk.Frame(self)
        
        self.f = plt.figure(figsize=(4,4),dpi=100)
        self.ax = self.f.add_subplot(111)
        
        canvas = FigureCanvasTkAgg(self.f, self)
        plt.tight_layout(rect=[0, 0, 1, 1])
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def render(self,HKkey,HKlims,x,y):
        
        ax = _ax_render_HK(self.ax,x,y,HKlims,HKkey)
        
        try: self.f.autofmt_xdate()
        except:
            pass
        plt.tight_layout(rect=[0, 0, 1, 1])
        
          



class HKButton(tk.Button,object):
    
    def __init__(self,master=None,**options):
        """ """
        
        ix = options.pop('ix') 
        ic = options.pop('ic') 
        ir = options.pop('ir')
        
        super(HKButton,self).__init__(master,options)
        
        self.ix = ix
        self.ir = ir
        self.ic = ic
        self.text = options['text']


class HKFlags(tk.Toplevel):
    """ """
    
    def __init__(self,root,parent,interval=5000,elvis=context.elvis):
        """ """
        
        tk.Toplevel.__init__(self,root)
        self.ncols = 5
        self.elvis = elvis
        
        self.parent = parent
        self.interval = interval
        self.HKkeys = self.parent.HKkeys[1:]
        
        self.wm_title('HK Flags')
        
        self.minsize(width=400,height=400)
        
        self.HKframe = tk.Frame(self)
        self.HKframe.pack(fill='both')
        
        self.setup_flagsgrid()
        
        self.update()
        
    def setup_flagsgrid(self):
        
        HKkeys = self.HKkeys
        
        ncols = self.ncols
        nrows = int(np.ceil(len(HKkeys)/(ncols*1.)))
        
        self.HKflags = []
        
        for ir in range(nrows):
            for ic in range(ncols):
                ix = ir * ncols + ic
                
                try: HKkey = HKkeys[ix]
                except:
                    break
                
                #self.setup_Flag(HKkey,'green',ix,ic,ir)
                
                self.HKflags.append(HKButton(self, text=HKkey,font=small_font,bg='green',
                                             ix=ix,ic=ic,ir=ir))
                self.HKflags[-1].grid(column=ic,row=ir,sticky='nsew',\
                     in_=self.HKframe)
                
                self.HKflags[-1].bind("<Button 1>",self.ResetFlag)
                self.HKflags[-1].bind("<Button 3>",self.showHK)
        
        if ic < ncols-1:
            icr = ic+1
            irr = ir
        else:
            icr = 0
            irr = ir+1
        
        self.resetButton = HKButton(self, text='RESET ALL',font=small_font,bg='blue',
                                             ix=-1,ic=icr,ir=irr)
        self.resetButton.grid(column=icr,row=irr,sticky='nsew',\
                     in_=self.HKframe)
                
        self.resetButton.bind("<Button 1>",self.ResetAllFlags)
    
    def setup_Flag(self,text,color,ix,ic,ir):
        try:
            HKflag = self.HKflags[ix]
        except IndexError:
            self.HKflags.append(HKButton(self,text=text,font=small_font,bg=color,
                                         ix=ix,ic=ic,ir=ir))
            HKflag = self.HKflags[ix]
        try: HKflag.grid(column=ic,row=ir,sticky='nsew',in_=self.HKframe)
        except: stop()
        HKflag.bind("<Button 1>", self.ResetFlag)
        HKflag.bind("<Button 3>)", self.showHK)
    
    def ResetFlag(self,event):
        """ """

        button = event.widget
        ix = button.ix
        
        self.changeColor(ix,'green')
    
    def ResetAllFlags(self,event):
        
        for HKflag in self.HKflags:
            ix = HKflag.ix
            self.changeColor(ix,'green')
            
    def showHK(self,event):
        
        HKflag = event.widget
        ix = HKflag.ix
        HKkey = self.HKkeys[ix]
        HKlim = self.parent.HKlims[HKkey]
        t = self.parent.HK['time'].copy()
        y = self.parent.HK[HKkey].copy()
        #print 'Im here!'
        
        window = SingleHKplot(self.parent.root)
        window.render(HKkey,HKlim,t,y)
        
    
    def update(self):
        
        try: self.find_offlims()
        except: pass    
        #self.find_offlims()
    
        #print self.interval
    
        self.after(self.interval,self.update)
            
    
    def find_offlims(self):
        
        HKkeys = self.HKkeys
        #print self.parent.HK
        
        for ix in range(len(HKkeys)):
            HKlim = self.parent.HKlims[HKkeys[ix]]
            lastval = self.parent.HK[HKkeys[ix]][-1]
            
            isWithin = self.validate(lastval,HKlim)
            if isWithin: continue
            
            self.changeColor(ix,'red')
    

    def validate(self,val,HKlim):
        
        if len(HKlim) == 2:
            if (HKlim[0]<= val) and (HKlim[1]>=val):
                return True
        elif len(HKlim) == 1:
            if HKlim == val:
                return True
        return False
    
    def changeColor(self,ix,color):
        
        ic = self.HKflags[ix].ic
        ir = self.HKflags[ix].ir
        text = self.HKflags[ix].text
        self.HKflags[ix].destroy()
        self.HKflags[ix] = HKButton(self, text=text,font=small_font,bg=color,
                                             ix=ix,ic=ic,ir=ir)
        self.HKflags[ix].grid(column=ic,row=ir,sticky='nsew',\
                in_=self.HKframe)
                
        self.HKflags[ix].bind("<Button 1>",self.ResetFlag)
        self.HKflags[ix].bind("<Button 3>",self.showHK)
        
        #self.setup_Flag(text,color,ix,ic,ir)
        

class HKDisplay(tk.Toplevel):
    """ """
    def __init__(self,root,path,interval,elvis=context.elvis):
        
        self.path = path
        self.interval = interval
        self.elvis = elvis
        self.date = '21-02-80'
        self.HKfile = None
        self.HK = dict()
        self.sizeHK = 0 # size in bytes of HK file
        self.page = 1
        
        #self.nHKlines = 1
        self.root = root
        
        self.HKkeys = HKtools.allHK_keys[elvis]
        self.HKlims = HKtools.HKlims[elvis]['P']
        
        self.search_HKfile()
        
        tk.Toplevel.__init__(self,self.root)
        
        self.minsize(width=700,height=925)
        
        self.wm_title('HK Display')
        
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
        
        self.setup_plotgrid()
        
    
    def setup_plotgrid(self):
        
        f = plt.figure(figsize=(4,8),dpi=100)
        
        axs = []
        
        for i in range(9):
            
            ax = f.add_subplot(3,3,i+1)            
            axs.append(ax)
        
        self.f = f
        self.axs = axs
        
        canvas = FigureCanvasTkAgg(self.f, self)
        plt.tight_layout(rect=[0, 0, 1, 1])
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    
    def search_HKfile(self):
        """ """
        
        struct_date = time.strptime(os.path.split(st.replace(self.path,'_',''))[-1],'%d%b%y')
        date_infile = time.strftime('%d-%m-%y',struct_date)
        
        self.date = date_infile

        tmp_HK = 'HK_%s_ROE1_01.txt' % date_infile
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
        
        if page * 9 > nkeys+9-1:
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
        
        #with open(self.HKfile) as f:
        #    nHKlines = len(f.readlines())
        
        #if nHKlines <= self.nHKlines:
        #    yield self.HK
        
        statinfo = os.stat(self.HKfile)
        sizeHK = statinfo.st_size
        
        if sizeHK <= self.sizeHK:
            yield self.HK
        
        self.sizeHK = sizeHK
        #self.nHKlines = nHKlines
        
        HK = HKtools.loadHK_QFM(self.HKfile,elvis=self.elvis)
        
        
        dtobjarr = np.array([datetime.datetime.strptime(item,'%d-%m-%y_%H:%M:%S') \
                             for item in HK['TimeStamp']])        
    
        pHK = dict(time = dtobjarr)
        
        subKeys = [Key for Key in HK.keys() if Key != 'TimeStamp']
        
        for key in subKeys:
            pHK[key] = HK[key].copy()
        
        pHK['CCD1_OD_T'] += datetime.datetime.now().second / 30. # TEST
        
        self.HK = pHK
        
        yield pHK
        
    def gen_render(self):
        
        HKlims = self.HKlims
        
        def render(pHK):
            
            
            max_xticks = 6
            xloc = plt.MaxNLocator(max_xticks)            
            
            nHK = len(self.HKkeys_to_plot)
            
            for item in self.axs:
                item.clear()
            
            if nHK <9:
                for item in self.axs[-(9-nHK)]:
                    item.get_xaxis().set_visible(False)
                    item.get_yaxis().set_visible(False)
            
            
            for i in range(nHK):
                
                ax = self.axs[i]
                
                HKname = self.HKkeys_to_plot[i]
                _HKlims = HKlims[HKname]
                                
                try:
                    x = pHK['time'].copy()
                    y = pHK[HKname].copy()
                except KeyError:
                    x = np.array([0,1])
                    y = np.array([0,1])
                
                ax = _ax_render_HK(ax,x,y,_HKlims,HKname)
                
            try: self.f.autofmt_xdate()
            except:
                pass
            plt.tight_layout(rect=[0, 0, 1, 1])
            

        return render
        

    def start_updating(self,interval):
        #self.do_update = True
        
        f = self.f
        render = self.gen_render()
        #while self.do_update:
        return animation.FuncAnimation(f, render,self.get_data,interval=interval)

