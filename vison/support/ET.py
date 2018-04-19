#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Module to issue WARNING / ALERT phone calls to designated phone numbers.

'... E.T. phone home...'

Created on Thu Sep 14 10:13:12 2017

@author: raf
"""

# IMPORT STUFF
from pdb import set_trace as stop
from twilio.rest import TwilioRestClient
import cPickle as pickle
import os
# END IMPORT


URLs = dict(
    NonSpecificWarning="https://visonwarningcall.000webhostapp.com/visonwarningcall.xml")


def grab_numbers_and_codes():
    """Retrieves phone numbers and access codes necessary to make the phone calls."""

    detailsf = os.path.join(os.path.pathsep, 'home', 'raf',
                            'CREDENTIALS', 'credentials_twilio.pick')[1:]
    details = pickle.load(open(detailsf))
    return details


class ET(object):
    """Class to do phone calls."""

    def __init__(self,):

        resources = grab_numbers_and_codes()
        self.TWILIO_PHONE_NUMBER = resources['TWILIO_PHONE_NUMBER']
        self.SIDnToken = resources['SIDnToken']
        self.DIAL_NUMBERS = resources['DIAL_NUMBERS']

    def dial_numbers(self, url):
        """Dials one or more phone numbers from a Twilio phone number.

        :param url: char, URL with the TwiML code that Twilio uses as instructions 
                    on call. Basically, it provides a message to be voiced, 
                    as intended.

        """

        SID, Token = self.SIDnToken
        TWILIO_PHONE_NUMBER = self.TWILIO_PHONE_NUMBER
        DIAL_NUMBERS = self.DIAL_NUMBERS

        client = TwilioRestClient(SID, Token)
        for number in DIAL_NUMBERS:
            print("Dialing " + number)
            # set the method to "GET" from default POST because Amazon S3 only
            # serves GET requests on files. Typically POST would be used for apps
            client.calls.create(to=number, from_=TWILIO_PHONE_NUMBER,
                                url=url, method="GET")


if __name__ == '__main__':

    ETavailable = True
    try:
        et = ET()
    except IOError:
        print 'vison.support.ET: Phone Calling is limited to personel with access details.'
        ETavailable = False

    if ETavailable:
        et.dial_numbers(URLs['NonSpecificWarning'])
