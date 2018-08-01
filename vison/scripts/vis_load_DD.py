#! $HOME/SOFTWARE/anaconda2/envs/vison/bin/ python
"""

VIS Ground Calibration Campaign

Loading a DataDict object for inspection.

Created on Wed Aug 01 10:00:00 2018

:autor: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
import sys
from optparse import OptionParser

from pylab import plot,show

from vison.support import files


# END IMPORT





if __name__ == '__main__':


    parser = OptionParser()
    parser.add_option("-f", "--file", dest="file",
                      default='', help="Pickle file with DataDict to inspect")

    (options, args) = parser.parse_args()

    if options.file == '':
        parser.print_help()
        sys.exit()
        
    dd = files.cPickleRead(options.file)
    stop()

