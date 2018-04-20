#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Module to handle HK-OOL Warnings

Created on Thu Apr 19 16:09:02 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT  STUFF
import os
from pdb import set_trace as stop
from vison.support import ET
# END IMPORT

critical_HKkeys = ['CCD3_TEMP_T', 'CCD2_TEMP_T', 'CCD1_TEMP_T', 
                   'VID_PCB_TEMP_T', 'FPGA_PCB_TEMP_T',
                   'CCD1_TEMP_B', 'CCD2_TEMP_B', 'CCD3_TEMP_B',
                   'VID_PCB_TEMP_B', 'FPGA_PCB_TEMP_B']


class Warnings(object):
    
    critical_HKkeys = critical_HKkeys
    severidict = dict()
    
    def __init__(self):
        """ """
    
    def process_event(self,HKKey,value,HKlim,violation_type):
        raise NotImplementedError
    
    def iscritical(self,HKkey):
        return bool(HKkey in self.critical_HKkeys)
    
    def assess_OOL_incident(self,HKkey,data=None):
        """ """
        if self.iscritical(HKkey):
            severity = self.severidict[HKkey]
            self.issue_warning(HKkey,data,severity)
    
    
    def send_email(self,subject,body,recipient):
        os.system('mail -s "%s" < %s %s' % (subject,body,recipient))
        # including an attachment (e.g. a figure) requires using "mutt" or some
        # or some other program that can handle that... example:
        #echo "This is the message body" | mutt -a file.to.attach -s "subject of message" recipient@domain.com
    
    def warn_via_email(self,HKkey,data=None):
        """ """
        
        #os.system('mail -s "ALARM (ALMOST) FULL DISK" < fullalarm.mail r.azzollini@ucl.ac.uk')
        self.send_email(subject,body,recipient,attachment)
        
        
    def warn_via_phone(self,HKkey):
        """ """
        
        
    def issue_warning(self,HKkey,data=None,severity=1):
        """ """
        self.warn_via_email(HKkey,data)
        if severity > 1:
            self.warn_via_phone(HKkey)
    
        