#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 19:11:53 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
import string as st
# END IMPORT

def add_checkHK_report(self,report_HK,tag):
        
    msg_HK = ['$\\bf{HK-%s}$, Listing parameters not passing Tests:\n' % tag]
    msg_HK += [st.join(['\\bf{%s}' % (st.replace(key,'_','\_'),) for key in report_HK if not report_HK[key]],', ')]
    if len(msg_HK[-1]) == 0: msg_HK += ['$\\bf{ALL\ PASSED}$\n']
    
    if self.report is not None:
        self.report.add_text(st.join(msg_HK,'\n'))
        self.report.add_text('\\')
    
    return msg_HK