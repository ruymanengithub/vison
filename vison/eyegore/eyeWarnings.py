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
import re

from vison.support import ET
from vison.support import vjson
from vison.support import utils

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

rootURL = "https://visonwarningcall.000webhostapp.com"
subURLs = dict(
    NonSpecificWarning="visonwarningcall.xml",
    LoCCDTemp="low_CCD_TEMP.xml",
    HiCCDTemp="hi_CCD_TEMP.xml",
    LoVIDTemp="low_VID_PCB_TEMP.xml",
    HiVIDTemp="hi_VID_PCB_TEMP.xml",
    LoFPGATemp="low_FPGA_PCB_TEMP.xml",
    HiFPGATemp="hi_FPGA_PCB_TEMP.xml")
URLs = dict()
for key,value in subURLs.iteritems():
    URLs[key] = '%s/%s' % (rootURL,value)


recipient = vjson.load_jsonfile(os.path.join(utils.credentials_path,'recipients_eyegore'))['main']


class EyeWarnings(object):
    
    critical_HKkeys = critical_HKkeys
    severitydict = severitydict
    
    def __init__(self,parent_flags_obj=None):
        """ """
        self.parent = parent_flags_obj
        if self.parent is not None: 
            self.log = self.parent.log
        self.recipient = recipient
        
        self.ETavailable = True
        try:
            self.et = ET.ET()
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
        
        self.send_email(subject,bodyList,self.recipient)
    
    def do_phone_call(self,url):
        """ """
        self.et.dial_numbers(url)
    
    def get_phone_url(self,HKkey,violation_type):
        """ """        
        HKkeys_patterns = dict(
        CCDTemp='CCD\d_TEMP_[T,B]',
        VIDTemp='VID_PCB_TEMP_[T,B]',
        FPGATemp='FPGA_PCB_TEMP_[T,B]')
        
        matches_expression = lambda pair: re.match(pair[0],pair[1]) is not None
                                                  
        
        def get_matching_HKurl(HKkey):
            for key,value in HKkeys_patterns.iteritems():
                if matches_expression((value,HKkey)): return key
            return None
        
        HKurl = get_matching_HKurl(HKkey)
        
        stop()
        
        violation_key = st.replace('T%i' % violation_type, '-', 'm')
        urlkey = '%s_%s' % (HKkey,violation_key)
        if urlkey not in URLdict:
            if self.log is not None:
                    self.log.info('VOICE WARNING not sent! urlkey not known: %s' % urlkey)
        return None
    
    def warn_via_phone(self,HKkey,violation_type):
        """ """
        
        if not self.ETavailable: 
            if self.log is not None:
                    self.log.info('VOICE WARNING not sent! ET not available')
            return None
        
        url = self.get_phone_url(HKkey,violation_type)

        try:
            self.do_phone_call(url)
        except:
            if self.log is not None:
                    self.log.info('VOICE WARNING not sent! [%s]' % HKkey)
    
    def get_parent_HK_data(self,HKkey):
        
        if self.parent is not None:
            HKdata = self.parent.HK # alias
            data = dict(time=HKdata['time'][-100:],
                        HKkey= HKdata[HKkey][-100:])
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
    

def test_URLs():
    """ """
    import urllib2
    
    for key,value in URLs.iteritems():
        try:
            urllib2.urlopen(value)
            print 'Found %s' % key, value
        except urllib2.HTTPError,e:
            print key, value, e.code
        except urllib2.URLError,e:
            print key, value, e.args

if __name__ == '__main__':
    test_URLs()
        