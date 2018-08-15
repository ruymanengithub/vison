#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 10:34:43 2018

:author: Ruyman Azzollini
:contact: r.azzollini_at_ucl.ac.uk

"""

# IMPORT STUFF
from pdb import set_trace as stop
from collections import OrderedDict
import numpy as np
# END IMPORT

CommonTaskInputs = OrderedDict(sorted([
    ('BLOCKID', ([str], 'Unique Detection Block Identifier.')),
    ('CHAMBER', ([str, type(None)], 'Test Chamber.')),
    ('processes', ([int], 'Number of Threads to use.')),
    ('datapath', ([str], 'Path to Data.')),
    ('diffvalues', ([
        dict], 'Dictionary with values of variable entries in acq. script (e.g. CCD S/Ns) [opt.].')),
    ('elvis', ([str], 'ELVIS version.')),
    ('explogf', ([str, list], 'Exposure Log File(s).')),
    ('ID', ([str], 'Pipeline Execution ID.')),
    ('inCDPs', ([
        dict], 'Input Calibration Data Products, used in processing of data.')),
    ('OBSID_lims',
     ([list], 'OBSID Limits to search for frames of Test in Exp-Log.')),
    ('perflimits', ([
        dict], '(Alt.) Performance limits: RON, BIAS levels, gradients, etc.')),
    ('resultsroot',
     ([str], 'Path to results directory (Pipeline Session Level).')),
    ('resultspath',
     ([str], 'Path to results directory (Task Level).')),
    ('subpaths',
     ([dict], 'Paths to results sub-folders. Relative to resultspath.')),
    ('structure', ([dict], 'Test Structure.')),
    ('test', ([str], 'Test/Task Name')),
    ('todo_flags', ([dict], 'Sub-Tasks to-do indicator.')),
    ('preprocessing', ([dict], 'Pre-processing kwargs dictionary.'))
]))


class Inputs(dict):
    """Class to hold, transfer and 'document' Task Inputs."""

    manifesto = OrderedDict([('a', ([int], 'first variable')),
                             ('b', ([int], 'second variable'))])

    def __init__(self, *args, **kwargs):
        
        for key in self.manifesto:
            kmv = self.manifesto[key]
            kmvl = list(kmv)
            if type(None) not in kmvl[0]:
                kmvl[0] += [type(None)]
            self.manifesto[key] = tuple(kmvl)
            
        self._fill_with_Nones()
        
        if 'args' in locals():
            for item in locals()['args']:
                _key, _val = item
                self.assertinputs(_key, _val)
        if 'kwargs' in locals():
            for key in kwargs:
                self.assertinputs(key, kwargs[key])
        super(Inputs, self).__init__(*args, **kwargs)
        
    def _fill_with_Nones(self):
        for key in self.manifesto.keys():
            self[key] = None
    

    def assertinputs(self, key, value):
        # return # TESTS
        assert key in self.manifesto, '"%s" not in %s' % (key, self.__class__)
        # assert isinstance(value,self.manifesto[key][0]), \
        #     '"%s" %s does not match expectation: %s' % (key,type(value),self.manifesto[key][0])
        assert np.any([isinstance(value, item) for item in self.manifesto[key][0]]), \
            '"%s" %s does not match expectation: %s' % (
                key, type(value), self.manifesto[key][0])

    def __setitem__(self, key, value):
        self.assertinputs(key, value)
        super(Inputs, self).__setitem__(key, value)

    def update(self, *args, **kwargs):
        if 'args' in locals():
            for argdict in locals()['args']:
                for key in argdict:
                    self.assertinputs(key, argdict[key])
        if 'kwargs' in locals():
            for key in kwargs:
                self.assertinputs(key, kwargs[key])
        super(Inputs, self).update(*args, **kwargs)

    def show_manifesto(self):
        return self.manifesto.__repr__()


class ChildInputs(Inputs):
    manifesto = dict(a=(float, 'first variable'),
                     b=(float, 'second variable'),
                     c=(int, 'third variable'))


if __name__ == '__main__':

    inputs = ChildInputs([('b', 2.)], a=1., c=2)

    stop()
