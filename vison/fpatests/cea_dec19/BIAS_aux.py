"""

Auxiliary Functions and resources to FPA_BIAS.

Created on Thu Nov 14 10:50:00 2019

:author: Ruyman Azzollini

"""


# IMPORT STUFF
from pdb import set_trace as stop
import numpy as np
import os
from collections import OrderedDict

from vison.plot import figclasses
from vison.datamodel import cdp


# END IMPORT

def _get_raw_img_dict(readmode, temperature):

    raw_img_dict = dict(
    figname='BIAS_Raw_Img_%s_%s.png' % (readmode, temperature),
    caption='FPA Bias Image Display: %s - %s' % (readmode, temperature),
    meta=dict(suptitle='%s - %s' % (readmode, temperature)))

    return raw_img_dict

def _get_DiffOffsetsMap_dict(readmode, temperature):
    """ """

    DiffOffsetsMap_dict = dict(
        figname='BIAS_%s_%s_DiffOffsetsMap.png' % (readmode, temperature),
        caption='Difference in Offset with Reference image.',
        meta=dict(suptitle='Diff. Offsets'))

    return DiffOffsetsMap_dict

def _get_RatioRonMap_dict(readmode, temperature):
    """ """

    RatioRonMap_dict = dict(
        figname='BIAS_%s_%s_RatioRonMap.png' % (readmode, temperature),
        caption='Ratio of RON measured / reference.',
        meta=dict(suptitle='Ratio RONs'))

    return RatioRonMap_dict


def _get_VPROFS_dict(readmode,temperature):
    """ """

    VPROFS_dict = dict(
    figname='BIAS_%s_%s_VPROFS.png' % (readmode, temperature),
    caption='Vertical Average Profiles (along columns).',
    meta=dict(
        suptitle='Avg. Column Offset Profiles.',
        doLegend=True,
        corekwargs=dict(PRE=dict(marker='', linestyle='-', color='r'),
            IMG=dict(marker='', linestyle='-', color='g'),
            OVE=dict(marker='', linestyle='-', color='b'))))

    return VPROFS_dict

def get_Bfigs(readmode, temperature):
    """ """

    Bfigs = dict()
    Bfigs['raw_img'] = [figclasses.Fig_Dynamic, \
            _get_raw_img_dict(readmode, temperature)]
    Bfigs['DIFFOFFSETSMAP'] = [figclasses.Fig_Dynamic, 
            _get_DiffOffsetsMap_dict(readmode,temperature)]
    Bfigs['RATIORONSMAP'] = [figclasses.Fig_Dynamic, 
            _get_RatioRonMap_dict(readmode,temperature)]
    Bfigs['VPROFS'] = [figclasses.Fig_Dynamic, 
            _get_VPROFS_dict(readmode,temperature)]
    Bfigs['BlueScreen'] = [figclasses.BlueScreen, dict()]

    return Bfigs


def get_CDP_lib(readmode,temperature):
    """ """

    OFF_cdp = cdp.Json_CDP()
    OFF_cdp.rootname = 'BIAS_%s_%s_OFFSETS' % (readmode, temperature)
    RON_cdp = cdp.Json_CDP()
    RON_cdp.rootname = 'BIAS_%s_%s_RONS' % (readmode, temperature)

    CDP_lib = dict(OFF_JSON=OFF_cdp,
                    RON_JSON=RON_cdp)

    return CDP_lib