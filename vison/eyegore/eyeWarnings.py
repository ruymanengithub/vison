#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Module to handle HK-OOL Warnings

Created on Thu Apr 19 16:09:02 2018

:author: Ruyman Azzollini

"""

# IMPORT  STUFF
from __future__ import print_function
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

severitydict = dict(CCD1_TEMP_T=dict(Tm1=2, T1=2),
                    CCD1_TEMP_B=dict(Tm1=2, T1=2),
                    CCD2_TEMP_T=dict(Tm1=2, T1=2),
                    CCD2_TEMP_B=dict(Tm1=2, T1=2),
                    CCD3_TEMP_T=dict(Tm1=2, T1=2),
                    CCD3_TEMP_B=dict(Tm1=2, T1=2),
                    VID_PCB_TEMP_T=dict(Tm1=1, T1=2),
                    VID_PCB_TEMP_B=dict(Tm1=1, T1=2),
                    FPGA_PCB_TEMP_T=dict(Tm1=1, T1=2),
                    FPGA_PCB_TEMP_B=dict(Tm1=1, T1=2),
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
for key, value in subURLs.items():
    URLs[key] = '%s/%s' % (rootURL, value)


try:
    recipients_dict = vjson.load_jsonfile(os.path.join(
        utils.credentials_path, 'recipients_eyegore'))
    recipients = [recipients_dict['main'], recipients_dict['secondary']]
except IOError:
    recipients = [None]


def matches_expression(pair): return re.match(pair[0], pair[1]) is not None


def _get_matching_HKurl(HKkey, patterns):
    for key, value in patterns.items():
        if matches_expression((value, HKkey)):
            return key
    return None


class EyeWarnings(object):

    critical_HKkeys = critical_HKkeys
    severitydict = severitydict

    def __init__(self, parent_flags_obj=None):
        """ """
        self.parent = parent_flags_obj
        self.log = None
        if self.parent is not None:
            self.log = self.parent.log

        self.recipients = recipients

        try:
            self.et = ET.ET()
        except IOError:
            print('vison.support.ET: Phone Calling is limited to personel with access details.')
            self.et = None

    def process_event(self, HKkey, violation_type, value, HKlim, timestamp):

        if not self.iscritical(HKkey):
            return None
        self.assess_OOL_incident(
            HKkey, violation_type, value, HKlim, timestamp)

    def iscritical(self, HKkey):
        return bool(HKkey in self.critical_HKkeys)

    def assess_OOL_incident(self, HKkey, violation_type, value, HKlim, timestamp):
        """ """
        if self.iscritical(HKkey):

            tvar = 'T%i' % violation_type
            violation_key = tvar.replace('-', 'm')
            Kseveritydict = self.severitydict[HKkey]

            try:
                severity = Kseveritydict[violation_key]
            except KeyError:
                severity = 0

            self.issue_warning(HKkey, severity, violation_type,
                               value, HKlim, timestamp)

    def send_email(self, subject, bodyList, recipient):

        with tempfile.NamedTemporaryFile(mode='w+a', delete=False) as f:
            for line in bodyList:
                print(line, file=f)
            f.close()
            try:
                os.system('mail -s "%s" %s < %s' %
                          (subject, recipient, f.name))
            except BaseException:
                print(('WARNING email not sent! [subject: %s]' % subject))
                if self.log is not None:
                    self.log.info(
                        'WARNING email not sent! [subject: %s]' % subject)
        os.unlink(f.name)

    def warn_via_email(self, HKkey, value, HKlim, timestamp, HKdata=None):
        """ """

        if self.recipients[0] is None:
            if self.log is not None:
                self.log.info(
                    "warn_via_email: recipient must be valid (%s)" % self.recipients.__str__())

        subject = 'Eyegore HK WARNING: %s (DT=%s)' % (HKkey, timestamp)
        bodyList = ['HK OOL WARNING: %s' % HKkey,
                    'value = %s, limits = %s' % (value, HKlim),
                    'at %s\n\n' % timestamp]
        if HKdata is not None:
            _data = OrderedDict(time=HKdata['time'])
            _data[HKkey] = HKdata[HKkey]
            df = pd.DataFrame.from_dict(_data)
            bodyList.append('LATEST VALUES after the jump\n\n')
            bodyList.append(df.to_string(index=False))

        for recipient in self.recipients:
            self.send_email(subject, bodyList, recipient)

    def do_phone_call(self, url):
        """Does phone call via self.et object"""
        self.et.dial_numbers(url)

    def get_phone_url(self, HKkey, violation_type):
        """ """
        HKkeys_patterns = dict(
            CCDTemp='CCD\d_TEMP_[T,B]',
            VIDTemp='VID_PCB_TEMP_[T,B]',
            FPGATemp='FPGA_PCB_TEMP_[T,B]')

        HKurl = _get_matching_HKurl(HKkey, HKkeys_patterns)

        if violation_type == -1:
            prefix = 'Lo'
        elif violation_type == 1:
            prefix = 'Hi'
        elif violation_type == 2:
            prefix = 'Ne'
        else:
            prefix = 'Un'

        urlkey = '%s%s' % (prefix, HKurl)

        if urlkey in URLs:
            return URLs[urlkey]

        return None

    def warn_via_phone(self, HKkey, violation_type):
        """ """

        if self.et is None:
            if self.log is not None:
                self.log.info('VOICE WARNING not sent! ET not available')
            return None

        url = self.get_phone_url(HKkey, violation_type)

        try:
            self.do_phone_call(url)
        except BaseException:
            if self.log is not None:
                self.log.info('VOICE WARNING not sent! [%s]' % HKkey)

    def warn_via_sms(self, HKkey, value, HKlim, timestamp):
        """ """

        if self.et is None:
            if self.log is not None:
                self.log.info('SMS WARNING not sent! ET not available.')

        body = ''.join(['HK OOL WARNING: %s. ' % HKkey,
                        'value = %s, limits = %s. ' % (value, HKlim),
                        'at %s' % timestamp])

        try:
            self.send_sms(body)
        except BaseException:
            if self.log is not None:
                self.log.info('SMS WARNING not sent! ET not available.')

    def send_sms(self, body):
        """Sends text message via self.et object"""
        self.et.send_sms(body)

    def get_parent_HK_data(self, HKkey):

        if self.parent is not None:
            HKdata = self.parent.HK  # alias
            data = {'time': HKdata['time'][-100:],
                    HKkey: HKdata[HKkey][-100:]}
        else:
            data = None

        return data

    def issue_warning(self, HKkey, severity, violation_type, value, HKlim, timestamp):
        """ """

        HKdata = self.get_parent_HK_data(HKkey)

        if severity > 0:
            self.warn_via_email(HKkey, value, HKlim, timestamp, HKdata)
        if severity > 1:
            #print 'WARNINGS via phone-calls/sms DISABLED by now in eyeWarnings'
            # self.warn_via_phone(HKkey,violation_type) DISABLED BY NOW
            self.warn_via_sms(HKkey, value, HKlim, timestamp)


def test_URLs():
    """ """
    import urllib.request, urllib.error, urllib.parse

    for key, value in URLs.items():
        try:
            urllib2.urlopen(value)
            print('Found %s' % key, value)
        except urllib2.HTTPError as e:
            print(key, value, e.code)
        except urllib2.URLError as e:
            print(key, value, e.args)


def test_get_URLs():
    ew = EyeWarnings()

    assert ew.get_phone_url('CCD1_TEMP_T', -1) == URLs['LoCCDTemp']
    assert ew.get_phone_url('CCD1_TEMP_T', 1) == URLs['HiCCDTemp']
    assert ew.get_phone_url('VID_PCB_TEMP_T', -1) == URLs['LoVIDTemp']
    assert ew.get_phone_url('VID_PCB_TEMP_T', 1) == URLs['HiVIDTemp']
    assert ew.get_phone_url('FPGA_PCB_TEMP_T', -1) == URLs['LoFPGATemp']
    assert ew.get_phone_url('FPGA_PCB_TEMP_T', 1) == URLs['HiFPGATemp']

    print('test_get_URLs passed!')


def test_do_phone_call(urlkey):
    ew = EyeWarnings()
    ew.do_phone_call(URLs[urlkey])


def test_do_send_sms(urlkey):
    ew = EyeWarnings()
    ew.send_sms(urlkey)


def test_assess_OOL():

    ew = EyeWarnings()
    ew.assess_OOL_incident('CCD1_TEMP_T', -1, -160,
                           [-133, 35], '24042018_1828')
    ew.assess_OOL_incident('VID_PCB_TEMP_T', -1, -15.,
                           [-10, 40], '24042018_1828')


if __name__ == '__main__':
    test_URLs()
    test_get_URLs()
    # test_assess_OOL() # does phone call
    # test_do_phone_call('LoCCDTemp') # does phone call
    # test_do_send_sms('LoCCDTemp')
