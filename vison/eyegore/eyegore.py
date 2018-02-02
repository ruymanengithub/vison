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
from vison import data as vdata
#import loadHK_QFM,allHK_keys

#from multiprocessing.dummy import Pool

import time
import Tkinter as tk
import ttk
import tkFont as tkFont
from PIL import Image, ImageTk

from eyeHK import HKDisplay,HKFlags
from eyeCCDs import ImageDisplay
from eyeObs import ExpLogDisplay
from vison.support import context
# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)

_extpath = os.path.join(os.sep,'data2','gaia','usdownloads','EuclidCaldata','Quarantine')

def rsync_to_remote(path):
    """ """
    extpath = os.path.join(_extpath,path)
    command = "rsync -avzq %s raf@msslus.ac.uk:%s" % (path+os.sep,extpath)
    #command = "rsync -avqz TEST_DATA/24_Feb_80 /home/raf/Desktop/24_Feb_80"
    os.system(command)

def rsync_to_altlocalpath(path,altpath):
    """ """
    _altpath = os.path.join(altpath,path)
    command = "rsync -avzq %s %s" % (path+os.sep,_altpath)
    print command
    os.system(command)

class Eyegore(tk.Tk):
    """ """
    
    def __init__(self,path,broadcast,intervals=[20000,20000,1000,20000,20000,20000],
                 elvis=context.elvis,dolite=False,altpath=''):
        """ """
        tk.Tk.__init__(self)
        
        if path[-1] == os.path.sep:
            path = path[:-1]
        self.path = path
        self.intervals = intervals
        self.broadcast = broadcast
        self.elvis = elvis
        self.dolite = dolite
        self.altpath = altpath
        
        #self.withdraw()
        
        self.setup_MasterWG()
        
        self.run()
        #self.mainloop()
    
    def setup_MasterWG(self):
        """ """
        
        self.wm_title('EYEGORE')
        
        fr = tk.Frame(self)
        fr.pack(fill='both',expand=True)
        
        eyegoregif = os.path.join(vdata.__path__[0],'Eyegore.gif')
        im = Image.open(eyegoregif)
        self.tkimg = ImageTk.PhotoImage(im)
        
        l1 = tk.Label(self,text="Path: ",font=LARGE_FONT)
        l1.grid(column=0,row=0,sticky='w',in_=fr)
        v1 = tk.Label(self,text=self.path,font=LARGE_FONT,bg='white',fg='black')
        v1.grid(column=1,row=0,sticky='w',in_=fr)
        
        l2 = tk.Label(self,text="Broadcasting: ",font=LARGE_FONT)
        l2.grid(column=0,row=1,sticky='w',in_=fr)
        v2 = tk.Label(self,text='%s' % self.broadcast,font=LARGE_FONT,bg='white',fg='black')
        v2.grid(column=1,row=1,sticky='w',in_=fr)

        end_button = tk.Button(self, text="EXIT", 
                              command=self.kill)
        end_button.grid(column=0,row=2,sticky='w',padx=0,in_=fr)
        
        #updates_button = tk.Button(self, text="NO-UPDATES", 
        #                      command=self.stop_updates)
        #updates_button.grid(column=1,row=2,sticky='w',padx=0,in_=fr)

        
        self.label = tk.Label(self,image=self.tkimg)
        self.label.grid(row=3,column=0,columnspan=2,sticky='nsew',in_=fr)
        
        fr.grid_columnconfigure(0, weight=1)
        fr.grid_rowconfigure(0, weight=1)
        
                
    def run(self):
        
        Ds = dict(image=ImageDisplay,hk=HKDisplay,
                  hkflags=HKFlags,explog=ExpLogDisplay)
        #dkeys = ['image','hk','hkflags','explog']
        
        if not self.dolite:
            display1 = Ds['image'](self,self.path,elvis=self.elvis)        
            ani1 = display1.start_updating(self.intervals[1])
        
        display2 = Ds['hk'](self,self.path,self.intervals[2],elvis=self.elvis)
        
        ani2 = display2.start_updating(self.intervals[3])
        
        display2b = Ds['hkflags'](self,display2,self.intervals[4],elvis=self.elvis)

        display4 = Ds['explog'](self,self.path,self.intervals[5],elvis=self.elvis)
        
        self.update()
        
        self.mainloop()
        
        
    def kill(self):
        self.destroy()
        sys.exit()
        
    def update(self):

        if self.broadcast:
            rsync_to_remote(self.path)
            
        if self.altpath != '':
            rsync_to_altlocalpath(self.path,self.altpath)

        self.after(self.intervals[0],self.update)

        


if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("-p","--path",dest="path",default='',help="day-path to be monitored.")
    parser.add_option("-B","--broadcast",dest="broadcast",action='store_true',default=False,help="")
    parser.add_option("-E","--elvis",dest="elvis",default=context.elvis,help="ELVIS version.")
    parser.add_option("-L","--lite",dest="lite",action="store_true",default=False,help="Run a lighter version of the program (no autom. HK-plots or image displays).")
    parser.add_option("-r","--rsync",dest="altpath",default='',help="rsync to an alternative local path.")
    
    (options, args) = parser.parse_args()
    
    if options.path == '':
        parser.print_help()
        sys.exit()
        
    path = options.path
    broadcast = bool(options.broadcast)
    elvis = options.elvis
    altpath = options.altpath
    dolite = bool(options.lite)
    
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
    
    
    app = Eyegore(path,broadcast=broadcast,elvis=elvis,dolite=dolite,altpath=altpath)

    
    
