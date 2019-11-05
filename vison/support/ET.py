#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Module to issue WARNING / ALERT phone calls to designated phone numbers.
Uses Twilio.

'... E.T. phone home...'

Created on Thu Sep 14 10:13:12 2017

:author: raf

"""

# IMPORT STUFF
from pdb import set_trace as stop
from twilio.rest import TwilioRestClient
import cPickle as pickle
import os
from vison.support import utils
# END IMPORT

testURL = "https://visonwarningcall.000webhostapp.com/visonwarningcall.xml"


def grab_numbers_and_codes():
    """Retrieves phone numbers and access codes necessary to make the phone calls."""

    detailsf = os.path.join(utils.credentials_path, 'credentials_twilio.pick')
    details = pickle.load(open(detailsf))
    return details


class ET(object):
    """Class to do phone calls."""

    def __init__(self,):

        resources = grab_numbers_and_codes()
        self.TWILIO_PHONE_NUMBER = resources['TWILIO_PHONE_NUMBER']
        self.SIDnToken = resources['SIDnToken']
        self.DIAL_NUMBERS = resources['DIAL_NUMBERS']
        self.client = TwilioRestClient(*self.SIDnToken)

    def dial_numbers(self, url):
        """Dials one or more phone numbers from a Twilio phone number.

        :param url: char, URL with the TwiML code that Twilio uses as instructions
                    on call. Basically, it provides a message to be voiced,
                    as intended.

        """

        for number in self.DIAL_NUMBERS:
            print("Dialing " + number)
            # set the method to "GET" from default POST because Amazon S3 only
            # serves GET requests on files. Typically POST would be used for apps
            self.client.calls.create(to=number, from_=self.TWILIO_PHONE_NUMBER,
                                     url=url, method="GET")

    def send_sms(self, body):
        """ """
        for number in self.DIAL_NUMBERS:
            print("Texting " + number)
            self.client.messages.create(
                to=number,
                from_=self.TWILIO_PHONE_NUMBER,
                body=body)


def test_do_phonecall():

    ETavailable = True
    try:
        et = ET()
    except IOError:
        print 'vison.support.ET: Phone Calling is limited to personel with access details.'
        ETavailable = False

    if ETavailable:
        et.dial_numbers(testURL)


def test_send_sms():
    ETavailable = True
    try:
        et = ET()
    except IOError:
        print 'vison.support.ET: Phone Calling is limited to personel with access details.'
        ETavailable = False

    if ETavailable:
        body = "Hellow World!"
        et.send_sms(body)


if __name__ == '__main__':

    # test_do_phonecall()

    test_send_sms()
