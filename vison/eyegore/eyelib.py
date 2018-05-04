#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 19 15:55:48 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""
#IMPORT STUFF
import os
import Tkinter as tk
import pprint
#END IMPORT

def get_bitsize_file(item):
    """ """
    return os.stat(item).st_size

def get_bitsize(filelist):
    """ """
    bitsize = 0
    for item in filelist:
        bitsize += get_bitsize_file(item)
    return bitsize


class HelpButton(tk.Button, object):

    def __init__(self, master=None, helplist=None, **kwargs):
        """ """
        
        if helplist is None:
            helplist = ['There is no help']
        
        options = dict(text='?',bg='pink')
        options.update(kwargs)
        
        super(HelpButton,self).__init__(master,options)

        self.helplist = helplist
        self.bind("<Button 1>", self.showHelp)
        
        
    def showHelp(self,event):
        """ """
        pp = pprint.PrettyPrinter()
        pp.pprint(['','','']+self.helplist)

