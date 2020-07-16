#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Script to build test SEQUENCES easily and reliably.

:History:
Created on Thu Mar 14 10:30:33 2019

:author: raf

"""

# IMPORT STUFF
from __future__ import print_function
from pdb import set_trace as stop
from optparse import OptionParser
import string as st
import sys
import os
import glob

from vison.support import vjson
# END IMPORT


def session_builder(inputs):
    """ """

    inpath = inputs['inpath']
    outpath = inputs['outpath']
    sess_dict = inputs['sessions'].copy()

    assert os.path.exists(inpath)

    if not os.path.exists(outpath):
        os.system('mkdir %s' % outpath)

    session_names = list(sess_dict.keys())

    for session in session_names:

        print(('session: %s' % session))

        tests = sess_dict[session]
        sesspath = os.path.join(outpath, session)
        if not os.path.exists(sesspath):
            os.system('mkdir %s' % sesspath)
        else:
            os.system('rm %s' % os.path.join(sesspath, '*'))

        scripts = []

        for test in tests:
            tmpscriptL = glob.glob(os.path.join(inpath, '*%s*xlsx' % test))
            if len(tmpscriptL) == 0:
                print(("Did not find test %s in session %s, Quitting!" %
                      (test, session)))
                sys.exit()
            if len(tmpscriptL) > 1:
                print(("Found more than 1 script for test %s in session %s, Quitting!" %
                       (test, session)))
            scripts.append(tmpscriptL[0])

        seq_f = os.path.join(sesspath, 'TEST_SEQUENCE_%s.txt' % session)

        with open(seq_f, 'w') as f:
            for script in scripts:
                os.system('cp %s %s/' % (script, sesspath))
                stripscript = os.path.split(script)[-1]
                print(stripscript, file=f)
                print(stripscript)
            f.close()


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-j", "--json", dest="json",
                      default='', help="json file with inputs")

    (options, args) = parser.parse_args()

    if options.json == '':
        parser.print_help()
        sys.exit()

    inputs = vjson.load_jsonfile(options.json, useyaml=True)

    # MISSING: INPUTS VALIDATION

    session_builder(inputs)
