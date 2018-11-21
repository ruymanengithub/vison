#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Eyegore: House Keeping Monitoring.

Created on Fri Oct 13 14:11:41 2017

:author: raf

"""

# IMPORT STUFF
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.animation as animation
from matplotlib import pyplot as plt
import glob

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
from vison.datamodel.HKtools import format_date, _ax_render_HK
from vison.eyegore.eyelib import get_bitsize, HelpButton
from vison.support import vistime
# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)


class SingleHKplot(tk.Toplevel):
    """ """

    def __init__(self, root):
        """ """

        tk.Toplevel.__init__(self, root)

        self.wm_title('HK Parameter')

        self.minsize(width=450, height=400)

        #frame = tk.Frame(self)

        self.f = plt.figure(figsize=(4, 4), dpi=100)
        self.ax = self.f.add_subplot(111)

        canvas = FigureCanvasTkAgg(self.f, self)
        plt.tight_layout(rect=[0, 0, 1, 1])
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2TkAgg(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def render(self, HKkey, HKlims, x, y):

        ax = _ax_render_HK(self.ax, x, y, HKlims, HKkey)
        try:
            plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
            self.f.autofmt_xdate()
        except:
            pass
        plt.tight_layout(rect=[0, 0, 1, 1])


class HKButton(tk.Button, object):

    def __init__(self, master=None, **kwargs):
        """ """

        options = dict(status=False)
        options.update(kwargs)

        ix = options.pop('ix')
        ic = options.pop('ic')
        ir = options.pop('ir')
        status = options.pop('status')

        super(HKButton, self).__init__(master, options)

        self.ix = ix
        self.ir = ir
        self.ic = ic
        self.text = options['text']
        self.status = status  # True/False/None


def validate_within_HKlim(val, HKlim):
    """ 
    violation:
        0: None
        -1: below lower limit
        1: above upper limit
        2: different from limit, if limit is a single value

    """

    if len(HKlim) == 2:
        if HKlim[0] > val:
            return False, -1
        elif HKlim[1] < val:
            return False, 1
    elif len(HKlim) == 1:
        if HKlim != val:
            return False, 2
    return True, 0


def sort_HKfiles(HKfiles):
    """ """

    # First we order by dates (ELVIS generates a
    # new master HK file after midnight in the same day-folder)

    alldays = []
    for iHKfile in HKfiles:
        itstr = os.path.split(iHKfile)[-1][3:11]
        it = datetime.datetime.strptime(itstr, '%d-%m-%y')
        alldays.append(it)
    dorder = np.argsort(alldays)

    oHKfiles = []

    # Then we sort by the appendix number among the master HK files of a
    # given date (there may be more than one, I've been told...)

    ordered_dates = np.unique(np.array(alldays)[dorder])

    for iod, idate in enumerate(ordered_dates):
        idate = ordered_dates[iod]
        idatestr = idate.strftime('%d-%m-%y')
        subfiles = [item for item in HKfiles if 'HK_%s' % idatestr in item]
        oHKfiles += np.sort(subfiles).tolist()

    return oHKfiles


flagshelpList = ['HK Flags:  HK parameters monitoring',
                 'Green means HK flag is within limits ("lowered")',
                 'Red means HK flag is outside limits ("raised")',
                 'Gray means HK flag is "muted" (not being monitorized)',
                 'Left-click: lower HK flag (reset to "green")',
                 'Right-click: show plot of HK parameter vs. time',
                 'Middle-click: mute/unmute HK parameter']


class HKFlags(tk.Toplevel):
    """ """

    def __init__(self, root, parent, interval=5000, elvis=context.elvis):
        """ """

        tk.Toplevel.__init__(self, root)
        self.ncols = 5
        self.elvis = elvis
        self.log = root.log
        self.Warnings = root.Warnings
        if self.Warnings is not None:
            self.Warnings.parent = self
        self.parent = parent
        self.interval = interval
        self.HKkeys = self.parent.HKkeys[1:]
        self.HK = self.parent.HK

        self.wm_title('HK Flags')

        self.minsize(width=400, height=400)

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

                try:
                    HKkey = HKkeys[ix]
                except IndexError:
                    break

                # self.setup_Flag(HKkey,'green',ix,ic,ir)

                self.HKflags.append(HKButton(self, text=HKkey, font=small_font, bg='green',
                                             ix=ix, ic=ic, ir=ir, status=False))
                self.HKflags[-1].grid(column=ic, row=ir, sticky='nsew',
                                      in_=self.HKframe)

                self.bind_buttons_to_methods(-1)

        if ic < ncols-1:
            icr = ic+1
            irr = ir
        else:
            icr = 0
            irr = ir+1

        self.resetButton = HKButton(self, text='RESET ALL', font=small_font, bg='blue',
                                    ix=-1, ic=icr, ir=irr)
        self.resetButton.grid(column=icr, row=irr, sticky='nsew',
                              in_=self.HKframe)

        self.resetButton.bind("<Button 1>", self.ResetAllFlags)

        self.HelpButton = HelpButton(self, helplist=flagshelpList, text='?',
                                     font=small_font, bg='pink',
                                     label='HK Flags: HELP')

        self.HelpButton.grid(column=icr+1, row=irr, sticky='nsew',
                             in_=self.HKframe)

    def ToggleMute(self, event):
        #ix = event.widget
        status = event.widget.status
        if status is None:
            self.UnmuteFlag(event)
        else:
            self.MuteFlag(event)

    def MuteFlag(self, event):
        """ """
        ix = event.widget.ix
        self.changeColor(ix, 'gray')
        self.HKflags[ix].status = None
        if self.log is not None:
            HKkey = self.HKflags[ix].text
            self.log.info('MUTING %s' % HKkey)

    def UnmuteFlag(self, event):
        """ """
        ix = event.widget.ix
        self.lowerflag(ix)
        if self.log is not None:
            HKkey = self.HKflags[ix].text
            self.log.info('UNMUTING %s' % HKkey)

    def ResetFlag(self, event):
        """ """
        button = event.widget
        ix = button.ix
        self.lowerflag(ix)
        if self.log is not None:
            HKkey = self.HKflags[ix].text
            self.log.info('LOWERED %s' % HKkey)

    def ResetAllFlags(self, event):
        for ix, HKflag in enumerate(self.HKflags):
            self.lowerflag(ix)
        if self.log is not None:
            self.log.info('LOWERED ALL FLAGS')

    def showHK(self, event):

        HKflag = event.widget
        ix = HKflag.ix
        HKkey = self.HKkeys[ix]
        HKlim = self.parent.HKlims[HKkey][1:]
        try:
            t = self.parent.HK['time'].copy()
            y = self.parent.HK[HKkey].copy()
            window = SingleHKplot(self.parent.root)
            window.render(HKkey, HKlim, t, y)
        except KeyError:
            print 'Cant find %s in HK!\n' % HKkey

    def update(self):

        try:
            self.find_offlims()
        except:
            pass

        self.after(self.interval, self.update)

    def find_offlims(self):

        HKkeys = self.HKkeys
        #print self.parent.HK

        for ix in range(len(HKkeys)):
            HKlim = self.parent.HKlims[HKkeys[ix]][1:]
            lastval = self.parent.HK[HKkeys[ix]][-1]

            isWithin, violation_type = validate_within_HKlim(lastval, HKlim)
            if isWithin or self.isflagraised(ix):
                continue
            time_stamp = vistime.get_time_tag()
            self.raiseflag(ix)
            if self.log is not None:
                self.log.info('RAISED %s: %s, lims=%s (timestamp=%s)' %
                              (HKkeys[ix], lastval.__str__(), HKlim.__str__(), time_stamp))
            if self.Warnings is not None and self.HKflags[ix].status is not None:
                self.Warnings.process_event(
                    HKkeys[ix], violation_type, lastval, HKlim, time_stamp)

    def isflagraised(self, ix):
        """ """
        return self.HKflags[ix].status

    def lowerflag(self, ix):
        """ """
        self.HKflags[ix].status = False
        self.changeColor(ix, 'green')

    def raiseflag(self, ix):
        """ """
        self.HKflags[ix].status = True
        self.changeColor(ix, 'red')

    def changeColor(self, ix, color):
        """ """
        ic = self.HKflags[ix].ic
        ir = self.HKflags[ix].ir
        text = self.HKflags[ix].text
        status = self.HKflags[ix].status
        self.HKflags[ix].destroy()
        self.HKflags[ix] = HKButton(self, text=text, font=small_font, bg=color,
                                    ix=ix, ic=ic, ir=ir, status=status)
        self.HKflags[ix].grid(column=ic, row=ir, sticky='nsew',
                              in_=self.HKframe)

        self.bind_buttons_to_methods(ix)

    def bind_buttons_to_methods(self, ix):
        """ """
        self.HKflags[ix].bind("<Button 1>", self.ResetFlag)
        self.HKflags[ix].bind("<Button 3>", self.showHK)
        self.HKflags[ix].bind("<Button 2>", self.ToggleMute)


class HKDisplay(tk.Toplevel):
    """ """

    def __init__(self, root, path, interval, elvis=context.elvis):

        self.path = path
        self.interval = interval
        self.elvis = elvis
        self.date = '21-02-80'
        self.HKfiles = None
        self.HK = dict()
        self.sizeHK = 0  # size in bytes of HK file
        self.page = 1
        self.log = root.log
        self.HKkeys_to_plot = []

        #self.nHKlines = 1
        self.root = root

        self.HKkeys = HKtools.allHK_keys[elvis]
        self.HKlims = HKtools.HKlims[elvis]['S']

        if self.path is not None:
            self.search_HKfiles()

        tk.Toplevel.__init__(self, self.root)

        self.minsize(width=700, height=925)

        self.wm_title('HK Display')

        frame = tk.Frame(self)
        l1 = tk.Label(frame, text="Page: ", font=LARGE_FONT)
        l1.pack(side="left")

        pvar = tk.StringVar()
        pvar.set('1')
        pentry = tk.Entry(frame, textvariable=pvar)

        pentry.pack(side="left")
        frame.pack(side="top")

        def update_page(event):
            try:
                self.page = int(pvar.get())
            except ValueError:
                self.page = 1

        pentry.bind('<Return>', update_page)

        self.setup_plotgrid()

    def setup_plotgrid(self):

        f = plt.figure(figsize=(4, 8), dpi=100)

        axs = []

        for i in range(9):

            ax = f.add_subplot(3, 3, i+1)
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

    def search_HKfiles(self):
        """ """

        if self.path is None:
            self.HKfiles = None
            return

        struct_date = time.strptime(os.path.split(
            st.replace(self.path, '_', ''))[-1], '%d%b%y')
        date_infile = time.strftime('%d-%m-%y', struct_date)

        self.date = date_infile

        #tmp_HK = 'HK_%s_ROE1_??.txt' % date_infile
        tmp_HK = 'HK_??-??-??_ROE1_??.txt'
        tmp_HK = os.path.join(self.path, tmp_HK)

        HKfiles = glob.glob(tmp_HK)

        HKfiles = sort_HKfiles(HKfiles)

        #isthere = os.path.exists(tmp_HK)
        arethere = len(HKfiles) > 0

        if arethere:
            self.HKfiles = HKfiles
        else:
            print 'HKfiles %s not found' % tmp_HK
            self.HKfiles = None

    def select_HKkeys(self):
        """ """

        allHKkeys = self.HKkeys
        nkeys = len(allHKkeys)-1
        page = self.page
        #print page

        if page * 9 > nkeys+9-1:
            return

        ix0 = (page-1)*9
        ix1 = ix0 + 9 if ix0+9 <= nkeys else None

        HKkeys_to_plot = allHKkeys[1:][ix0:ix1]

        self.HKkeys_to_plot = HKkeys_to_plot

    def get_data(self):
        """ """

        self.select_HKkeys()

        if self.HKfiles is None:
            yield self.HK
            return

        sizeHK = get_bitsize(self.HKfiles)

        if sizeHK <= self.sizeHK:
            yield self.HK
            return

        self.sizeHK = sizeHK

        print 'loading HK...'

        try: HK = HKtools.loadHK_QFM(self.HKfiles, elvis=self.elvis)
        except IOError:
            yield self.HK
            return

        print 'done loading HK!'

        dtobjarr = np.array([datetime.datetime.strptime(item, '%d-%m-%y_%H:%M:%S')
                             for item in HK['TimeStamp']])

        pHK = dict(time=dtobjarr)

        subKeys = [Key for Key in HK.keys() if Key != 'TimeStamp']

        for key in subKeys:
            pHK[key] = HK[key].data.copy()

        self.HK = pHK.copy()

        yield pHK

    def gen_render(self):

        HKlims = self.HKlims

        def render(pHK):

            nHK = len(self.HKkeys_to_plot)

            for item in self.axs:
                item.clear()

            if nHK < 9:
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
                    x = np.array([0, 1])
                    y = np.array([0, 1])

                self.axs[i] = _ax_render_HK(ax, x, y, _HKlims, HKname)

            try:
                plt.gca().xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(format_date))
                self.f.autofmt_xdate()
            except:
                pass
            plt.tight_layout(rect=[0, 0, 1, 1])

        return render

    def start_updating(self, interval):
        f = self.f
        render = self.gen_render()
        return animation.FuncAnimation(f, render, self.get_data, interval=interval)
