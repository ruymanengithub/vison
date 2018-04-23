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
import string as st
import tempfile
import pandas as pd
from collections import OrderedDict

from vison.support import ET
# END IMPORT

critical_HKkeys = ['CCD3_TEMP_T', 'CCD2_TEMP_T', 'CCD1_TEMP_T', 
                   'VID_PCB_TEMP_T', 'FPGA_PCB_TEMP_T',
                   'CCD1_TEMP_B', 'CCD2_TEMP_B', 'CCD3_TEMP_B',
                   'VID_PCB_TEMP_B', 'FPGA_PCB_TEMP_B']

severitydict = dict(CCD1_TEMP_T=dict(Tm1=2,T1=2),
                    CCD1_TEMP_B=dict(Tm1=2,T1=2),
                    CCD2_TEMP_T=dict(Tm1=2,T1=2),
                    CCD2_TEMP_B=dict(Tm1=2,T1=2),
                    CCD3_TEMP_T=dict(Tm1=2,T1=2),
                    CCD3_TEMP_B=dict(Tm1=2,T1=2),
                    VID_PCB_TEMP_T=dict(Tm1=1,T1=2), 
                    VID_PCB_TEMP_B=dict(Tm1=1,T1=2), 
                    FPGA_PCB_TEMP_T=dict(Tm1=1,T1=2),
                    FPGA_PCB_TEMP_B=dict(Tm1=1,T1=2), 
                                    )

# NEEDS WORK, populate URLdict
URLdict = dict(NonSpecificWarning="https://visonwarningcall.000webhostapp.com/visonwarningcall.xml")
masked_recipient = 'ruyman.azzollini_at_gmail.com' # NEEDS BETTER/SAFER SOLUTION!


class EyeWarnings(object):
    
    critical_HKkeys = critical_HKkeys
    severitydict = severitydict
    
    def __init__(self,parent_flags_obj=None):
        """ """
        self.parent = parent_flags_obj
        self.log = self.parent.log
        
        self.ETavailable = True
        try:
            self.et = ET()
        except IOError:
            print 'vison.support.ET: Phone Calling is limited to personel with access details.'
            self.ETavailable = False
    
    def process_event(self,HKkey,violation_type,value,HKlim,timestamp):
        
        if not self.iscritical(HKkey):
            return None
        self.assess_OOL_incident(HKkey,violation_type,value,HKlim,timestamp)
        
    
    def iscritical(self,HKkey):
        return bool(HKkey in self.critical_HKkeys)
    
    def assess_OOL_incident(self,HKkey,violation_type,value,HKlim,timestamp):
        """ """
        if self.iscritical(HKkey):
            violation_key = st.replace('T%i' % violation_type, '-', 'm')            
            Kseveritydict = self.severitydict[HKkey]
            try: 
                severity = Kseveritydict[violation_key]
            except KeyError:
                severity = 0
            
            self.issue_warning(HKkey,severity,violation_type,value,HKlim,timestamp)
    
    
    def send_email(self,subject,bodyList,recipient):
        
        with tempfile.NamedTemporaryFile(mode='w+a',delete=False) as f:
            for line in bodyList: print >> f, line
            f.close()
            try:
                os.system('mail -s "%s" %s < %s' % (subject,recipient,f.name))
            except:
                if self.log is not None:
                    self.log.info('WARNING email not sent! [subject: %s]' % subject)
        os.unlink(f.name)
        
    
    def warn_via_email(self,HKkey,value,HKlim,timestamp,HKdata=None):
        """ """
        
        recipient = st.replace(masked_recipient,'_at_','@') # NEEDS IMPROVEMENT, UNSAFE
        
        subject = 'Eyegore HK WARNING: %s (DT=%s)' % (HKkey,timestamp)
        bodyList = ['HK OOL WARNING: %s' % HKkey,
            'value = %s, limits = %s' % (value,HKlim),
            'at %s\n\n' % timestamp]
        if HKdata is not None:
            _data = OrderedDict(time=HKdata[HKkey])
            _data[HKkey] = HKdata[HKkey]
            df = pd.DataFrame.from_dict(_data)
            bodyList.append('LATEST VALUES after the jump\n\n')
            bodyList.append(df.to_string(index=False))
        
        self.send_email(subject,bodyList,recipient)
    
    def do_phone_call(self,urlkey):
        """ """
        self.et.dial_numbers(URLdict[urlkey])
    
    def warn_via_phone(self,HKkey,violation_type):
        """ """
        
        if not self.ETavailable: 
            if self.log is not None:
                    self.log.info('VOICE WARNING not sent! ET not available')
            return None
        
        violation_key = st.replace('T%i' % violation_type, '-', 'm')   
        urlkey = '%s_%s' % (HKkey,violation_key)
        if urlkey not in URLdict:
            if self.log is not None:
                    self.log.info('VOICE WARNING not sent! urlkey not known: %s' % urlkey)
            return None

        try:
            self.do_phone_call(urlkey)
        except:
            if self.log is not None:
                    self.log.info('VOICE WARNING not sent! [%s]' % HKkey)
    
    def get_parent_HK_data(self,HKkey):
        
        if self.parent is not None:
            HKdata = self.parent.parent.HK # alias
            data = dict(time=HKdata['time'][-100:],
                        HKkey= HKdata['time'][-100:])
        else:
            data = None
            
        return data
        
    def issue_warning(self,HKkey,severity,violation_type,value,HKlim,timestamp):
        """ """
        
        HKdata = self.get_parent_HK_data(HKkey)
        
        if severity>0:
            self.warn_via_email(HKkey,value,HKlim,timestamp,HKdata)
        if severity > 1:
            self.warn_via_phone(HKkey)
    
        