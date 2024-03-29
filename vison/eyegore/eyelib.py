#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""


Eyegore: common use library.

Created on Thu Apr 19 15:55:48 2018

:author: Ruyman Azzollini

"""
# IMPORT STUFF
from pdb import set_trace as stop
import os
import tkinter as tk
import pprint
import string as st
# END IMPORT


def get_bitsize_file(item):
    """ """
    return os.stat(item).st_size


def get_bitsize(filelist):
    """ """
    bitsize = 0
    for item in filelist:
        try:
            nbitsize = get_bitsize_file(item)
        except BaseException:
            nbitsize = 0
        bitsize += nbitsize
    return bitsize


class HelpButton(tk.Button, object):

    def __init__(self, master=None, helplist=None, **kwargs):
        """ """

        if helplist is None:
            helplist = ['There is no help...', 'Who needs help anyway?',
                        'You are better off without help!']

        options = dict(text='?', bg='pink', label='HELP')
        options.update(kwargs)

        label = options.pop('label')

        super(HelpButton, self).__init__(master, options)

        self.helplist = helplist
        self.bind("<Button 1>", self.showHelp)
        self.master = master
        self.label = label

    def showHelp(self, event):
        """ """
        #pp = pprint.PrettyPrinter()
        # pp.pprint(['','','']+self.helplist)
        newwin = tk.Toplevel(self.master)
        #newwin = tk.Toplevel()
        newwin.wm_title(self.label)
        display = tk.Label(newwin, text='\n'.join(self.helplist))
        display.pack()
