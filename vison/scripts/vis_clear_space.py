#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Script to clean-up processing files produced by the pipeline.
Useful to free-up space.

**USE WITH CAUTION**.

:History:
Created on Thu Dec  6 10:19:24 2018

:author: raf

"""

# IMPORT STUFF
from pdb import set_trace as stop
from optparse import OptionParser
import sys
import os
import glob
import string as st
# END IMPORT


def find_and_erase(path, keyword):
    """ """

    execline1 = "find %s -type d -name '%s'" % (path, keyword)
    print('\nDirectories that will be WIPED OUT clear:\n')
    os.system(execline1)
    print('\n')
    ans1 = eval(input('Are you happy with the selection? yes=y/Y '))
    if ans1.lower() != 'y':
        sys.exit()

    selpaths = os.popen(execline1).read().split( '\n')
    selpaths = [item for item in selpaths if len(item) > 0]

    print('\nPreparing to CLEAR OUT:\n')
    totNfiles = 0
    totsize = 0.
    for selpath in selpaths:
        Nselfiles = sum([len(files) for directory, folders, files in os.walk(selpath)])
        selsize = sum([sum([os.path.getsize(os.path.join(directory, fname)) for fname in files])
                       for directory, folders, files in os.walk(selpath)])
        print('%s, %s files, %.1e bytes' % (selpath, Nselfiles, selsize))
        totNfiles += Nselfiles
        totsize += selsize
    print('\nTOTAL: %s files, %.1e bytes\n' % (totNfiles, totsize))

    execline2 = "find %s -type d -name '%s' -exec sh -c '%s' {} \;" % (
        path, keyword, 'rm -r "$0"/*')

    print(('\nTo be executed: "%s"\n' % execline2))

    ans2 = eval(input('Still want to proceed? yes=y/Y '))


    if ans2.lower() != 'y':
        sys.exit()

    #print execline2
    os.system(execline2)


if __name__ == '__main__':

    usage = 'usage: %prog [options]'
    parser = OptionParser(usage)
    parser.add_option("-p", "--path", dest="path", default='',
                      help="path where to look for files. Search IS recursive.")
    parser.add_option("-k", "--keyword", dest="keyword", default='',
                      help="keyword to identify target directories. CANT BE EMPTY.")

    (options, args) = parser.parse_args()

    if options.keyword == '':
        parser.print_help()
        sys.exit()

    path = options.path
    keyword = options.keyword

    find_and_erase(path, keyword)
