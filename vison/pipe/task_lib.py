#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 19:11:53 2017

:author: raf
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import string as st
from vison.datamodel import HKtools
# END IMPORT




def check_HK(self,HKKeys,reference='command',limits='P',tag='',doReport=False,doLog=True):
    """ """
    
    checkers = dict(command=HKtools.check_HK_vs_command,
                    abs=HKtools.check_HK_abs)
    
    report_HK = checkers[reference](HKKeys,self.dd,limits=limits,elvis=self.elvis) 
    
    if doReport and self.report is not None:
        nicetag = st.replace(tag,' ','\ ')
        msg_HK = ['$\\bf{HK-%s}$, Parameters not passing Tests:\n' % nicetag]
        msg_HK += [st.join(['$\\bf{\\textcolor{red}{%s}}$' % (st.replace(key,'_','\_'),) for key in report_HK if not report_HK[key]],', ')]
        if len(msg_HK[-1]) == 0: msg_HK += ['$\\bf{ALL\ PASSED}$\n']
        self.report.add_Text(st.join(msg_HK,'\n'))
        self.report.add_Text('\\')
        
    if doLog and self.log is not None:
        msg_HK = ['HK-%s, Parameters not passing Tests:\n' % tag]
        msg_HK += [st.join([key for key in report_HK if not report_HK[key]],', ')]
        if len(msg_HK[-1]) == 0: msg_HK += ['ALL PASSED']
        self.log.info(msg_HK)
    
    return report_HK


def add_checkHK_report(self,report_HK,tag):
        
    msg_HK = ['$\\bf{HK-%s}$, Listing parameters not passing Tests:\n' % tag]
    msg_HK += [st.join(['\\bf{%s}' % (st.replace(key,'_','\_'),) for key in report_HK if not report_HK[key]],', ')]
    if len(msg_HK[-1]) == 0: msg_HK += ['$\\bf{ALL\ PASSED}$\n']
    
    if self.report is not None:
        self.report.add_Text(st.join(msg_HK,'\n'))
        self.report.add_Text('\\')
    
    return msg_HK

