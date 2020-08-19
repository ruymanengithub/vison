#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python
"""

Finds DataDicts (picks) in a path and converts all string arrays 
in them to decoded UTF8 arrays. Needed when using DataDicts
generated with vison:py2.7 with vison:py3.6.



:History:
Created on Thu Aug 19 12:51:00 2020

:autor: Ruyman Azzollini

"""

# IMPORT STUFF
from pdb import set_trace as stop
import sys
from optparse import OptionParser
from glob import glob
import os

from vison.support import files


# END IMPORT


if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option("-d", "--directory", dest="directory",
        default = './', help = 'directory where it will search for DataDict pick files.')
    parser.add_option("-C", "--change", dest="doChange", action="store_true",
        default=False, help="effect Changes. If not set, program will just list DataDict pick files.")

    (options, args) = parser.parse_args()

    doChange = options.doChange
    directory = options.directory

    subfolders = '%s**%s' % (os.path.sep,os.path.sep)
    allDataDicts = glob(os.path.join(directory+subfolders,'*_DataDict.pick'), recursive=True)

    print('DataDicts found:\n')
    for ddname in allDataDicts: print(ddname)
    print('\n')

    if not doChange:
        sys.exit()

    ans = input('Want to update these files? y/n ').lower()

    if ans == 'y':
        print('Hey!')


