#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 15:20:01 2018

@author: raf
"""

# IMPORT STUFF
import os
from pdb import set_trace as stop
import string as st
from optparse import OptionParser
import sys

from vison.datamodel.EXPLOGtools import loadExpLog, mergeExpLogs
from vison.support import context
from astropy.table import vstack
# END IMPORT


def _loader(ELf, elvis):
    if ELf[0] == '#':
        return None
    try:
        return loadExpLog(ELf, elvis=elvis)
    except:
        return None


def explog_merger(ELlist, output='EXP_LOG_merged.txt', elvis=context.elvis):
    """ """

    f = open(ELlist)
    ELfs = f.readlines()
    f.close()
    ELfs = [item[0:-1] for item in ELfs]

    ELs = [_loader(ELf, elvis=elvis) for ELf in ELfs]
    ELs = [item for item in ELs if item is not None]

    EL = mergeExpLogs(ELs, addpedigree=False, verbose=True)

    EL.write(output, format='ascii', overwrite=False)


if __name__ == '__main__':
    """ """

    parser = OptionParser()

    parser.add_option("-L", "--list", dest="ELlist", default='None',
                      help="Text file with exposure logs to be merged.")
    parser.add_option("-O", "--output", dest="output",
                      default='EXP_LOG_merged.txt', help="Output Merged EXP LOG.")
    parser.add_option("-E", "--elvis", dest="elvis",
                      default=context.elvis, help="ELVIS version.")

    (options, args) = parser.parse_args()

    if options.ELlist == 'None':
        parser.print_help()
        sys.exit()

    ELlist = options.ELlist
    output = options.output
    elvis = options.elvis

    explog_merger(ELlist, output, elvis)
