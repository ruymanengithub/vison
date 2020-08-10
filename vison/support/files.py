"""
IO related functions.

:requires: PyFITS
:requires: NumPy

:author: Sami-Matias Niemi

"""
import datetime
import pickle
import os
import numpy as np
from pdb import set_trace as stop


def cPickleDumpDictionary(dictionary, output, protocol=2):
    """
    Dumps a dictionary of data to a cPickled file.

    :param dictionary: a Python data container does not have to be a dictionary
    :param output: name of the output file

    :return: None
    """
    out = open(output, 'wb')
    pickle.dump(dictionary, out, protocol=protocol)
    out.close()


def cPickleRead(ffile):
    """
    Loads data from a pickled file.
    """
    with open(ffile, 'rb') as f:
        inpick = pickle.load(f)
    f.close()
    return inpick


def cPickleDump(data, output, protocol=2):
    """
    Dumps data to a cPickled file.

    :param data: a Python data container
    :param output: name of the output file

    :return: None
    """
    out = open(output, 'wb')
    pickle.dump(data, out, protocol=protocol)
    out.close()


def convert_fig_to_eps(figname):
    """Converts a figure to .eps. Returns new file name."""
    root = os.path.splitext(figname)
    epsname = '%s.eps' % root
    os.system('convert %s %s' % (figname, epsname))
    return epsname


def test():

    data = ['a']
    picklef = 'test_cpickle.pick'

    cPickleDump(data, picklef)
    data = cPickleRead(picklef)
    stop()


if __name__ == '__main__':
    test()
