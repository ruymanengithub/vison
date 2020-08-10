#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Eyegore: CCDs display.

Created on Fri Oct 13 16:16:08 2017

:author: raf

"""

# IMPORT STUFF
#import matplotlib
# matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.animation as animation
from matplotlib import pyplot as plt


from pdb import set_trace as stop
from optparse import OptionParser
import sys
import os
import numpy as np
import time
import string as st
import datetime
import glob
from skimage import exposure

#from multiprocessing.dummy import Pool

import tkinter as tk
import tkinter.ttk
import tkinter.font as tkFont
#from PIL import Image, ImageTk

from vison.datamodel import ccd as ccdmod
from vison.support import context

# END IMPORT


LARGE_FONT = ("Helvetica", 12)
small_font = ("Verdana", 8)


class ImageDisplay(tk.Toplevel):
    """ """

    def __init__(self, parent, path, elvis=context.elvis, tag=''):
        """ """

        tk.Toplevel.__init__(self, parent)
        self.parent = parent
        self.path = path
        self.tag = tag

        title = 'Image Display'
        if self.tag != '':
            title = '%s: %s' % (title, self.tag)

        self.wm_title(title)

        self.log = parent.log

        self.minsize(width=850, height=400)

        f1 = self.setup_fig()

        canvas = FigureCanvasTkAgg(self.f, self)
        canvas.show()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def setup_fig(self):
        """ """
        f1 = plt.figure(figsize=(8, 3), dpi=100)
        ax1 = f1.add_subplot(131)
        ax1.set_title('CCD1')
        ax2 = f1.add_subplot(132)
        ax2.set_title('CCD2')
        ax3 = f1.add_subplot(133)
        ax3.set_title('CCD3')

        self.f = f1
        self.axs = [ax1, ax2, ax3]

        return

    def get_data(self):
        """"""
        IMGs = dict(ObsID=-1)
        for CCD in [1, 2, 3]:
            IMGs['CCD%i' % CCD] = None

        if self.path is None:
            yield IMGs
            return

        FITSs = sorted(glob.glob(os.path.join(
            self.path, '*.fits')), key=os.path.getmtime)
        lastFITS = FITSs[-1]
        lastObs = int(os.path.split(lastFITS)[-1].split('_')[1])

        for CCD in [1, 2, 3]:
            try:
                CCDF = glob.glob(os.path.join(self.path, 'EUC_%i_*D*T_ROE1_CCD%i.fits' %
                                              (lastObs, CCD)))[0]
                img = ccdmod.CCD(
                    CCDF, extensions=[-1]).extensions[0].data.copy()
            except BaseException:
                img = None
            IMGs['CCD%i' % CCD] = img
        IMGs['ObsID'] = lastObs

        yield IMGs

    def gen_render(self):
        """ """

        def render(IMGs):

            ObsID = IMGs['ObsID']
            t = datetime.datetime.now()
            s = t.strftime('%H:%M:%S')

            for CCD in [1, 2, 3]:
                self.axs[CCD - 1].clear()
                img = IMGs['CCD%i' % CCD]
                if img is None:
                    img = np.zeros((2119 * 2, 2086 * 2), dtype='int32')
                eimg = exposure.equalize_hist(img, nbins=2**16)
                self.axs[CCD - 1].imshow(eimg.T, cmap=plt.cm.gray, origin='lower left')
                self.axs[CCD - 1].get_xaxis().set_visible(False)
                self.axs[CCD - 1].get_yaxis().set_visible(False)
                self.axs[CCD - 1].set_title('CCD%i' % CCD)
            self.f.suptitle('OBS %i : %s' % (ObsID, s))
            self.f.tight_layout()

        return render

    def start_updating(self, interval=5000):
        return animation.FuncAnimation(self.f, self.gen_render(), self.get_data, interval=interval)
