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


from pdb import set_trace as stop
from optparse import OptionParser
import sys
import os

from vison import data as vdata
#import loadHK_QFM,allHK_keys

#from multiprocessing.dummy import Pool

import Tkinter as tk
from PIL import Image, ImageTk

from vison.eyegore.eyeHK import HKDisplay, HKFlags
from vison.eyegore.eyeCCDs import ImageDisplay
from vison.eyegore.eyeObs import ExpLogDisplay
from vison.eyegore import eyeWarnings
from vison.support import context
from vison.support import logger as lg
from vison import __version__
from vison.support import vistime
# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)

_extpath = os.path.join(os.sep, 'data2', 'gaia',
                        'usdownloads', 'EuclidCaldata', 'Quarantine')


def rsync_to_remote(path, broadcast):
    """ """
    #extpath = os.path.join(_extpath,path)
    command = "rsync -e 'ssh' -L -avzq %s raf@msslus:%s%s" % (
        path, broadcast, os.sep)
    print 'syncing to MSSLUS: %s' % command
    os.system(command)


def rsync_to_altlocalpath(path, altpath):
    """ """
    #_altpath = os.path.join(altpath, path)
    command = "rsync -avzq %s %s" % (path, altpath+os.sep)
    print 'SYNCING TO ALTLOCAL: %s' % command
    os.system(command)


class Eyegore(tk.Tk):
    """ """

    def __init__(self, path, broadcast, intervals=None,
                 elvis=context.elvis, blind=False, dolite=False, altpath='',
                 doWarnings=False, dolog=False, ds9target='*', tag=''):
        """ """
        tk.Tk.__init__(self)

        if intervals is None:
            intervals = dict(EYE=20000, 
         image=20000,
         HK=5000, 
         HKplots=500, 
         HKflags=500, 
         EXP=20000)

        if path is not None:
            if path[-1] == os.path.sep:
                path = path[:-1]
        self.path = path
        self.intervals = intervals
        self.broadcast = broadcast
        self.elvis = elvis
        self.blind = blind
        self.dolite = dolite
        self.altpath = altpath
        self.ds9target = ds9target
        self.tag = tag

        if dolog:
            datestamp = vistime.get_time_tag()
            
            if self.tag != '':
                ntag = '_%s' % self.tag
            else:
                ntag = ''
            self.logf = 'Eyegore_%s%s.log' % (datestamp,ntag)

            if os.path.exists(self.logf):
                os.system('rm %s' % self.logf)

            self.log = lg.setUpLogger(self.logf)
            self.log.info(['\n\nStarting Eyegore',
                           'DATE: %s' % datestamp,
                           'vison version: %s\n' % __version__])
        else:
            self.log = None

        if doWarnings:
            self.Warnings = eyeWarnings.EyeWarnings()
        else:
            self.Warnings = None

        self.setup_MasterWG()

        self.run()

    def setup_MasterWG(self):
        """ """
        
        title = 'EYEGORE'
        if self.tag != '':
            title = '%s: %s' % (title,self.tag)

        self.wm_title(title)
        

        fr = tk.Frame(self)
        fr.pack(fill='both', expand=True)
        
        eyegoregif = os.path.join(vdata.__path__[0], 'Eyegore.gif')
        im = Image.open(eyegoregif)
        self.tkimg = ImageTk.PhotoImage(im)

        l1 = tk.Label(self, text="Path: ", font=LARGE_FONT)
        l1.grid(column=0, row=0, sticky='w', in_=fr)
        v1 = tk.Label(self, text=self.path, font=LARGE_FONT,
                      bg='white', fg='black')
        v1.grid(column=1, row=0, sticky='w', in_=fr)

        l2 = tk.Label(self, text="Broadcasting: ", font=LARGE_FONT)
        l2.grid(column=0, row=1, sticky='w', in_=fr)
        v2 = tk.Label(self, text='%s' % self.broadcast,
                      font=LARGE_FONT, bg='white', fg='black')
        v2.grid(column=1, row=1, sticky='w', in_=fr)

        end_button = tk.Button(self, text="EXIT",
                               command=self.kill)
        end_button.grid(column=0, row=2, sticky='w', padx=0, in_=fr)

        # updates_button = tk.Button(self, text="NO-UPDATES",
        #                      command=self.stop_updates)
        # updates_button.grid(column=1,row=2,sticky='w',padx=0,in_=fr)

        self.label = tk.Label(self, image=self.tkimg)
        self.label.grid(row=3, column=0, columnspan=2, sticky='nsew', in_=fr)

        fr.grid_columnconfigure(0, weight=1)
        fr.grid_rowconfigure(0, weight=1)
        

    def run(self):

        Ds = dict(image=ImageDisplay, hk=HKDisplay,
                  hkflags=HKFlags, explog=ExpLogDisplay)
        #dkeys = ['image','hk','hkflags','explog']

        if not (self.dolite or self.blind):
            display1 = Ds['image'](self, self.path, elvis=self.elvis, tag=self.tag)
            ani1 = display1.start_updating(self.intervals['image'])

        display2 = Ds['hk'](
            self, self.path, self.intervals['HK'], elvis=self.elvis, 
                                dolite=self.dolite or self.blind,
                                tag=self.tag)
        
        if not (self.dolite or self.blind):
            ani2 = display2.start_updating_display(self.intervals['HKplots'])

        display2b = Ds['hkflags'](
            self, display2, self.intervals['HKflags'], elvis=self.elvis,tag=self.tag)
        
        if not self.dolite:

            display4 = Ds['explog'](
                self, self.path, self.intervals['EXP'], elvis=self.elvis, ds9target=self.ds9target,
                                               tag=self.tag)
        self.update()

        self.mainloop()

    def kill(self):
        self.destroy()
        sys.exit()

    def update(self):

        if self.broadcast is not None:
            rsync_to_remote(self.path, self.broadcast)

        if self.altpath != '':
            rsync_to_altlocalpath(self.path, self.altpath)

        self.after(self.intervals['EYE'], self.update)


def Eexecuter():

    parser = OptionParser()
    parser.add_option("-p", "--path", dest="path",
                      default=None, help="day-path to be monitored.")
    parser.add_option("-B", "--broadcast", dest="broadcast", default='None',
                      help="Synchronize data to gateway folder at msslus")
    parser.add_option("-E", "--elvis", dest="elvis",
                      default=context.elvis, help="ELVIS version.")
    parser.add_option("-b", "--blind", dest="blind", action="store_true", default=False,
                      help="Run without image or HK displays.")
    parser.add_option("-L", "--lite", dest="lite", action="store_true", default=False,
                      help="Run a lighter version of the program (no image/HK displays and no ExpLog).")
    parser.add_option("-r", "--rsync", dest="altpath", default='',
                      help="rsync to an alternative local path.")
    parser.add_option("-g", "--log", dest="dolog", action="store_true", default=False,
                      help="keep a log")
    parser.add_option("-W", "--Warnings", dest="doWarnings", action="store_true", default=False,
                      help="Raise warnings (via email and/or phone) if critical HK is OOL.")
    parser.add_option("-d", "--DS9", dest="ds9target", default='*',
                      help="Specify DS9 target (pyds9)?")
    parser.add_option("-t","--tag", dest="tag", default="",
                      help="add a tag to the name of windows")
    (options, args) = parser.parse_args()

    if options.path is None:
        parser.print_help()
        ans = raw_input(
            '\nNo path provided... you sure? [press any key to continue] \n')

        # sys.exit()

    path = options.path
    broadcast = options.broadcast
    elvis = options.elvis
    altpath = options.altpath
    blind = bool(options.blind)
    dolite = bool(options.lite)
    dolog = bool(options.dolog)
    doWarnings = bool(options.doWarnings)
    ds9target = options.ds9target
    tag = options.tag
    

    if broadcast != 'None':
        broadcast = os.path.join(_extpath, broadcast)
    else:
        broadcast = None

    # if not os.path.exists(path):
    #    sys.exit('path "%s" does not exist' % path)
    # else:
    header = '\n\n#########################################################\n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#       starting EYEGORE on path:                       #\n' +\
             '            "%s"                                         \n' +\
             '#                                                       #\n' +\
             '#                                                       #\n' +\
             '#########################################################\n'

    print header % path

    app = Eyegore(path, broadcast=broadcast, elvis=elvis, blind=blind,
                  dolite=dolite,altpath=altpath, doWarnings=doWarnings, 
                  dolog=dolog, ds9target=ds9target, tag=tag)


if __name__ == '__main__':

    Eexecuter()
