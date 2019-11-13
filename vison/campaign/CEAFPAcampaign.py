#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Description of the CEA FPA Campaign.


Created on Wed Oct 02 17:26:00 2019

:author: Ruyman Azzollini

"""

# IMPORT STUFF
from collections import OrderedDict
from pdb import set_trace as stop
import numpy as np
import copy

#from vison.pipe import lib as pilib
#from vison.support import context
from vison.support import utils

from vison.fpatests.cea_dec19 import FWD_WARM
from vison.fpatests.cea_dec19 import FPA_BIAS

# END IMPORT


def generate_test_sequence(toGen, elvis='FPA', FPAdesign='final'):
    """Now supporting test repetitions."""
    taskslist = toGen.keys()
    test_sequence = OrderedDict()
    for taskname in taskslist:
        if not toGen[taskname]:
            continue

        strip_taskname, iteration = utils.remove_iter_tag(taskname, Full=True)
        _toGen = OrderedDict()
        _toGen[strip_taskname] = True
        ans = _generate_test_sequence(_toGen, elvis=elvis, FPAdesign=FPAdesign)

        if iteration is not None:
            for key in ans.keys():
                test_sequence['%s.%i' % (key, iteration)] = copy.deepcopy(ans[key])
        else:
            for key in ans.keys():
                test_sequence[key] = copy.deepcopy(ans[key])

    return test_sequence


def _generate_test_sequence(toGen, elvis='FPA', FPAdesign='final'):
    """ """

    #print 'GENERATING TEST SEQUENCE...'

    test_sequence = OrderedDict()

    _toGen = dict(FWD_WARM=False,
                  BIAS_RWDVS_WARM=False,
                  BIAS_RWDV_WARM=False,
                  BIAS_RWDVS_COLD=False,
                  BIAS_RWDV_COLD=False,
                  BIAS_FWD_COLD=False)

    _toGen.update(toGen)

    commoninputs = dict(elvis=elvis,
                        FPAdesign=FPAdesign)

    # DARK-CURRENT RAMP

    if _toGen['FWD_WARM']:

        fwd_warm_inp = dict(
            test='FWD_WARM')
        fwd_warm_inp.update(commoninputs)

        fwd_warm = FWD_WARM.FWD_WARM(inputs=fwd_warm_inp.copy())

        #structBIAS01 = bias01.build_scriptdict(elvis=elvis)
        test_sequence['FWD_WARM'] = copy.deepcopy(fwd_warm)

    if _toGen['BIAS_RWDVS_WARM']:
        
        rwdvs_warm_inp = dict(
            test='BIAS_RWDVS_WARM',
            temperature='WARM',
            readmode='RWDVS')

        rwdvs_warm_inp.update(commoninputs)

        rwdvs_warm = FPA_BIAS.FPA_BIAS(inputs=rwdvs_warm_inp.copy())     

        test_sequence['BIAS_RWDVS_WARM'] = copy.deepcopy(rwdvs_warm)

    if _toGen['BIAS_RWDV_WARM']:

        
        rwdv_warm = FPA_BIAS.FPA_BIAS(inputs=rwdv_warm_inp.copy())     

        test_sequence['BIAS_RWDV_WARM'] = copy.deepcopy(rwdv_warm)

    if _toGen['BIAS_RWDVS_COLD']:


        rwdvs_cold = FPA_BIAS.FPA_BIAS(inputs=rwdvs_cold_inp.copy())     

        test_sequence['BIAS_RWDVS_COLD'] = copy.deepcopy(rwdvs_cold)

    if _toGen['BIAS_RWDV_COLD']:

        rwdv_cold = FPA_BIAS.FPA_BIAS(inputs=rwdv_cold_inp.copy())     

        test_sequence['BIAS_RWDV_COLD'] = copy.deepcopy(rwdv_cold)

    if _toGen['BIAS_FWD_COLD']:

        fwdv_cold = FPA_BIAS.FPA_BIAS(inputs=fwd_cold_inp.copy())     

        test_sequence['BIAS_FWD_COLD'] = copy.deepcopy(fwd_cold)    
        
    


    return test_sequence
